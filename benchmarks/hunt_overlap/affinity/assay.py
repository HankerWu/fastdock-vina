"""Run the fixed-input Gal-1 / Gal-3 panel and audit all retained parent poses."""

import argparse
import csv
import hashlib
import itertools
import json
import os
import random
import re
import subprocess
import tarfile
import time
from pathlib import Path

STATES = {"Gal1": ("1GZW", "4Q26"), "Gal3": ("1KJL", "3T1L")}
COUNT_KEYS = (
    "hunt_calls",
    "hunt_evals",
    "second_calls",
    "second_evals",
    "explicit_calls",
    "explicit_evals",
)


def write(path, rows):
    with path.open("w", newline="") as stream:
        out = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        out.writeheader()
        out.writerows(rows)


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(args):
    archive = args.repo / "benchmarks/galectin_affinity/evidence.tar.xz"
    with tarfile.open(archive) as tar:
        members = [
            m
            for m in tar
            if m.name == "panel.csv"
            or (
                m.name.startswith("prepared/")
                and m.name.endswith((".pdbqt", ".anchor", ".box"))
            )
        ]
        tar.extractall(args.work, members=members, filter="data")


def cells(args):
    prepare(args)
    bins = {
        k: getattr(args, k)
        for k in ("upstream", "counts", "parent", "candidate", "control")
        if getattr(args, k)
    }
    cases = list(itertools.product(sum(STATES.values(), ()), range(5, 18)))
    rng = random.Random(337383)
    rng.shuffle(cases)
    records = []
    data = args.work / "prepared"
    raw = args.work / "raw"
    raw.mkdir(exist_ok=True)
    for rec, cid in cases:
        arms = list(bins)
        rng.shuffle(arms)
        for arm in arms:
            key = f"{rec}-{cid}-{args.seed}-{arm}-{args.cap}"
            out = raw / (key + ".pdbqt")
            cmd = [
                str(bins[arm]),
                "--seed",
                str(args.seed),
                "--exhaustiveness",
                "1",
                "--max_evals",
                str(args.cap),
                "--num_modes",
                "9",
                "--energy_range",
                "100",
                "--cpu",
                "1",
                "--verbosity",
                "0",
                "--receptor",
                str(data / (rec + ".pdbqt")),
                "--ligand",
                str(data / f"{cid}.pdbqt"),
                "--config",
                str(data / (rec + ".box")),
                "--out",
                str(out),
            ]
            if arm in ("parent", "candidate", "control"):
                cmd.append("--search_counters")
            start = time.perf_counter()
            log = subprocess.check_output(cmd, text=True, stderr=subprocess.STDOUT)
            seconds = time.perf_counter() - start
            (raw / (key + ".log")).write_text(log)
            scores = re.findall(
                r"^REMARK VINA RESULT:\s+([-\d.]+)", out.read_text(), re.MULTILINE
            )
            counts = dict(re.findall(r"\b(\w+)=(\d+)", log))
            records.append(
                dict(
                    key=key,
                    arm=arm,
                    receptor=rec,
                    compound=cid,
                    seed=args.seed,
                    cap=args.cap,
                    top_score=scores[0],
                    modes=len(scores),
                    seconds=f"{seconds:.7f}",
                    sha256=sha(out),
                    **{
                        k: int(counts[k]) if arm != "upstream" else ""
                        for k in COUNT_KEYS
                    },
                )
            )
        if (
            "upstream" in bins
            and "counts" in bins
            and (
                len(
                    {
                        r["sha256"]
                        for r in records
                        if r["receptor"] == rec
                        and r["compound"] == cid
                        and r["arm"] in ("upstream", "counts")
                    }
                )
                != 1
            )
        ):
            raise ValueError("passive counters changed upstream pool")
    tag = f"{args.seed}-{args.cap}-" + "-".join(bins)
    write(args.work / ("runs-" + tag + ".csv"), records)
    (args.work / ("binaries-" + tag + ".json")).write_text(
        json.dumps({k: sha(p) for k, p in bins.items()}, indent=2) + "\n"
    )


def fixed(args):
    raw = args.work / "raw"
    data = args.work / "prepared"
    frozen = args.work / "frozen.pdbqt"
    rows = []
    for path in sorted(raw.glob(f"*-{args.seed}-parent-{args.cap}.pdbqt")):
        rec, cid, *_ = path.stem.split("-")
        for rank, mdl in enumerate(path.read_text().split("ENDMDL")[:-1], 1):
            frozen.write_text(
                "\n".join(
                    l for l in mdl.splitlines() if l and not l.startswith("MODEL")
                )
                + "\n"
            )
            cli = [
                "--score_only",
                "--cpu",
                "1",
                "--seed",
                str(args.seed),
                "--unbound_energy",
                "0",
                "--pose_coordinates",
                str(frozen),
                "--receptor",
                str(data / (rec + ".pdbqt")),
                "--config",
                str(data / (rec + ".box")),
                "--ligand",
                str(data / (cid + ".pdbqt")),
            ]
            terms = []
            start = time.perf_counter()
            for binary in filter(None, (args.parent, args.candidate, args.control)):
                out = subprocess.check_output([str(binary), *cli], text=True)
                terms.append(
                    [
                        l
                        for l in out.splitlines()
                        if "(kcal/mol)" in l or "Ligand - " in l or "Flex   - " in l
                    ]
                )
            if len(terms[0]) != 10 or any(t != terms[0] for t in terms[1:]):
                raise ValueError((path, rank, terms))
            rows.append(
                {
                    "key": path.stem,
                    "rank": rank,
                    "pose_sha256": sha(frozen),
                    "terms": " ".join(
                        re.search(r":\s+([+.\d-]+)", l)[1] for l in terms[0]
                    ),
                    "seconds": round(time.perf_counter() - start, 6),
                }
            )
    write(args.work / "fixed_scores.csv", rows)
    frozen.unlink()


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=("run", "fixed"))
    p.add_argument("repo", type=Path)
    p.add_argument("work", type=Path)
    for key in ("upstream", "counts", "parent", "candidate", "control"):
        p.add_argument("--" + key, type=Path)
    p.add_argument("--seed", type=int, default=912743)
    p.add_argument("--cap", type=int, default=12000)
    p.add_argument("--cpu", type=int)
    a = p.parse_args()
    if a.cpu is not None:
        os.sched_setaffinity(0, {a.cpu})
    a.work.mkdir(parents=True, exist_ok=True)
    (cells if a.action == "run" else fixed)(a)

"""End-to-end runs, counterbalanced within each ligand/seed/budget block."""

import argparse
import csv
import hashlib
import itertools
import json
import os
import random
import re
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "macrocycle_labels"))
from run_labels import REFINE, TRAJECTORY


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def write(path, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, rows[0], lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def cases(work, study):
    if study == "kd":
        d = work / "prepared"
        for r in read(HERE.parent / "search_backtracking/assay/panel.csv"):
            pdb, cid = r["pdb"], r["cid"]
            yield (
                pdb,
                cid,
                d / (pdb + ".pdbqt"),
                d / (cid + ".pdbqt"),
                d / (pdb + ".box"),
            )
        return
    b = work / "tutorial/solution"
    yield (
        "5QCT",
        "E4Y",
        b / "BACE_1_receptor.pdbqt",
        b / "BACE_1_ligand.pdbqt",
        b / "BACE_1_receptor_vina_box.txt",
    )
    h = work / "hsp/1YET"
    yield (
        "1YET",
        "GDM",
        h / "receptor.pdbqt",
        h / "ligand.pdbqt",
        work / "hsp/docked/box.txt",
    )
    for pdb in ("1FKF", "1FKB"):
        d = work / "extra" / pdb
        meta = json.loads((d / "preparation.json").read_text())
        if (
            meta["receptor_status"] != "ready"
            or " CG0" not in (d / "ligand.pdbqt").read_text()
        ):
            print(
                f"{pdb} ineligible: {meta['receptor_status']}, check closure tag",
                flush=True,
            )
            continue
        box = d / "box.txt"
        box.write_text(
            "".join(
                f"center_{k} = {v}\nsize_{k} = 26\n"
                for k, v in zip("xyz", meta["center"])
            )
        )
        yield (pdb, meta["ligand"], d / "receptor.pdbqt", d / "ligand.pdbqt", box)


def run(a):
    os.sched_setaffinity(0, {a.cpu})
    a.output.mkdir(parents=True, exist_ok=True)
    raw = a.output / "raw"
    raw.mkdir(exist_ok=True)
    arms = {"upstream": a.upstream, "parent": a.parent, "candidate": a.candidate}
    seeds, caps, modes = (
        ((740317,), (12000,), 9)
        if a.study == "kd"
        else ((458921, 740317), (3000, 12000), 5)
    )
    if a.study in ("tail", "tail-cost"):
        arms.pop("upstream")
        seeds = (*seeds, 191009)
    panel = list(cases(a.work, a.study))
    if a.study == "tail-cost":
        panel = [c for c in panel if c[0] == "1FKB"]
        caps = (2000, 3000)
    cells = list(itertools.product(panel, seeds, caps))
    random.Random(740317).shuffle(cells)
    rows = []
    for idx, ((pdb, cid, rec, lig, box), seed, cap) in enumerate(cells):
        keys = list(arms)
        j = idx % len(keys)
        for arm in keys[j:] + keys[:j]:
            key = f"{pdb}_{cid}_{seed}_{cap}_{arm}"
            f = raw / (key + ".pdbqt")
            args = [
                "--receptor",
                rec,
                "--ligand",
                lig,
                "--config",
                box,
                "--seed",
                str(seed),
                "--cpu",
                "1",
                "--exhaustiveness",
                "1",
                "--max_evals",
                str(cap),
                "--num_modes",
                str(modes),
                "--energy_range",
                "100",
                "--verbosity",
                "0",
            ]
            cmd = [arms[arm], *args, "--out", f]
            if arm != "upstream":
                cmd.append("--search_counters")
            start = time.perf_counter_ns()
            log = subprocess.check_output(cmd, text=True, stderr=subprocess.STDOUT)
            elapsed = time.perf_counter_ns() - start
            (raw / (key + ".log")).write_text(log)
            if arm == "upstream":
                replay = raw / (key + ".replay")
                log = subprocess.check_output(
                    [a.upstream_counts, *args, "--out", replay],
                    text=True,
                    stderr=subprocess.STDOUT,
                )
                if f.read_bytes() != replay.read_bytes():
                    raise RuntimeError("Upstream counting replay differs: " + key)
                replay.unlink()
                (raw / (key + ".counts")).write_text(log)
            t, e = TRAJECTORY.search(log), REFINE.search(log)
            values = re.findall(
                r"^REMARK VINA RESULT:\s+([-+.\d]+)", f.read_text(), re.M
            )
            rows.append(
                dict(
                    key=key,
                    pdb=pdb,
                    compound=cid,
                    seed=seed,
                    cap=cap,
                    arm=arm,
                    elapsed_ns=elapsed,
                    first_score=values[0],
                    modes=len(values),
                    hunt_calls=t[1],
                    hunt_evals=t[2],
                    second_calls=t[3],
                    second_evals=t[4],
                    explicit_calls=e[1],
                    explicit_evals=e[2],
                    joint=int(t[2]) + int(t[4]) + int(e[2]),
                    pool_sha256=hashlib.sha256(f.read_bytes()).hexdigest(),
                )
            )
    write(a.output / "runs.csv", rows)
    provenance = {
        k: hashlib.sha256(v.read_bytes()).hexdigest()
        for k, v in {**arms, "upstream_counts": a.upstream_counts}.items()
    }
    (a.output / "binaries.json").write_text(json.dumps(provenance, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("study", choices=("macro", "kd", "tail", "tail-cost"))
    for arg in ("work", "output", "upstream", "upstream_counts", "parent", "candidate"):
        p.add_argument(arg, type=Path)
    p.add_argument("--cpu", type=int, default=min(os.sched_getaffinity(0)))
    run(p.parse_args())

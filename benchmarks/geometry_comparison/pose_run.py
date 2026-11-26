"""Paired cognate and reciprocal docking on frozen experimental structures."""

import argparse
import csv
import hashlib
import itertools
import json
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "paired_crossdocking"))
from run import alignment, coords


def load_meta(work, pdb):
    meta = json.loads((work / pdb / "metadata.json").read_text())
    source = (work / pdb / (pdb + ".pdb")).read_text().splitlines()
    meta["ca"] = {
        line[22:27].strip() + ":" + line[17:20]: coords(line)
        for line in source
        if line.startswith("ATOM  ")
        and line[21] == meta["chain"]
        and line[12:16].strip() == "CA"
        and line[16] in " A"
    }
    return meta


def align_with_numbering(query, receiver):
    # Use author numbering first. For constructs with disjoint numbering,
    # identify a single gapless offset from unique sequence fragments. Neither
    # ligand coordinates nor docking outputs determine correspondence.
    try:
        return (*alignment(query, receiver), 0)
    except ValueError:

        def residues(meta):
            return {
                int(key.split(":")[0]): key.split(":")[1]
                for key in meta["ca"]
                if key.split(":")[0].lstrip("-").isdigit()
            }

        def fragments(aa):
            result = {}
            for i in aa:
                fragment = tuple(aa.get(j) for j in range(i, i + 12))
                if None not in fragment:
                    result.setdefault(fragment, []).append(i)
            return result

        source, dest = residues(query), residues(receiver)
        a, b = fragments(source), fragments(dest)
        offsets = {
            b[k][0] - a[k][0]
            for k in a.keys() & b.keys()
            if len(a[k]) == len(b[k]) == 1
        }
        if len(offsets) != 1:
            raise
        (offset,) = offsets
        if offset == 0:
            raise
        shared = [i for i in source if i + offset in dest]
        if len(shared) < 30 or sum(
            source[i] == dest[i + offset] for i in shared
        ) < 0.95 * len(shared):
            raise
        ca = {
            f"{i + offset}:{source[i]}": query["ca"][f"{i}:{source[i]}"]
            for i in shared
            if source[i] == dest[i + offset]
        }
        return (*alignment(dict(query, ca=ca), receiver), offset)


def write_csv(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def measure(job, exhaustiveness=3):
    (
        work,
        output,
        binary,
        query,
        receiver,
        ref,
        backbone,
        n_ca,
        offset,
        seed,
        label,
        cap,
    ) = job
    mode = "armijo" if label == "armijo" else "geometry"
    q, r = query["pdb"], receiver["pdb"]
    prefix = output / f"{q}-{r}-{seed}-{label}"
    pose = prefix.with_suffix(".pdbqt")
    cmd = [
        str(binary),
        "--receptor",
        str(work / r / "receptor.pdbqt"),
        "--ligand",
        str(work / q / "ligand.pdbqt"),
        "--pose_selection",
        "late",
        "--refinement_optimizer",
        mode,
        "--exhaustiveness",
        str(exhaustiveness),
        "--max_evals",
        str(cap),
        "--num_modes",
        "9",
        "--energy_range",
        "100",
        "--search_counters",
        "--cpu",
        "1",
        "--seed",
        str(seed),
        "--out",
        str(pose),
    ]
    for ax, pos in zip("xyz", receiver["center"]):
        cmd += ["--center_" + ax, str(pos), "--size_" + ax, "24"]
    start = time.monotonic()
    rproc = subprocess.run(
        cmd, capture_output=True, text=True, timeout=300, check=False
    )
    wall = time.monotonic() - start
    log = rproc.stdout + rproc.stderr
    prefix.with_suffix(".log").write_text(log)
    rproc.check_returncode()
    counters = re.search(r"Search counters: .*explicit_evals=(\d+) ", log)
    search = [int(value) for value in re.findall(r"(?:hunt|second)_evals=(\d+)", log)]
    if counters is None or len(search) != 2 * exhaustiveness:
        raise ValueError("Missing effort counters in " + str(prefix))
    cell = {
        "target": query["target"],
        "query": q,
        "receiver": r,
        "seed": seed,
        "arm": label,
        "kind": "redock" if q == r else "cross",
        "aligned_ca": n_ca,
        "ca_rmsd": backbone,
        "numbering_offset": offset,
        "max_evals": cap,
    }
    modes = []
    for block in pose.read_text().split("ENDMDL"):
        atoms = {
            int(l[6:11]): coords(l)
            for l in block.splitlines()
            if l.startswith(("ATOM  ", "HETATM"))
        }
        if atoms:
            xyz = np.array([atoms[i] for i in query["serial"]])
            rmsd = np.sqrt(((ref - xyz) ** 2).sum(axis=2).mean(axis=1)).min()
            score = float(
                next(
                    l.split()[3]
                    for l in block.splitlines()
                    if l.startswith("REMARK VINA RESULT:")
                )
            )
            modes.append(dict(**cell, rank=len(modes) + 1, rmsd_A=rmsd, score=score))
    if not modes:
        raise ValueError("Empty output at " + str(prefix))
    cell.update(
        top1_A=modes[0]["rmsd_A"],
        oracle_A=min(m["rmsd_A"] for m in modes),
        score=modes[0]["score"],
        modes=len(modes),
        search_evals=sum(search),
        explicit_evals=int(counters[1]),
        wall_s=wall,
        pdbqt_sha256=hashlib.sha256(pose.read_bytes()).hexdigest(),
    )
    return cell, modes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "work", type=Path, help="Prepared fixtures containing manifest.json"
    )
    parser.add_argument("binary", type=Path, help="Analytic build")
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--finite-binary", type=Path, help="142c676 finite-difference control"
    )
    parser.add_argument("--max-evals", type=int, default=12000)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--seeds", type=int, nargs="+", default=[63419, 63421])
    args = parser.parse_args()
    work, out = args.work.resolve(), args.output.resolve()
    (out / "poses").mkdir(parents=True, exist_ok=True)
    rows = json.loads((work / "manifest.json").read_text())
    meta = [load_meta(work, row["pdb"]) for row in rows if row["status"] == "ready"]
    jobs, exclusions = [], []
    for query, receiver in itertools.product(meta, repeat=2):
        if query["target"] != receiver["target"]:
            continue
        same = query["pdb"] == receiver["pdb"]
        rot, tr, rms, n, offset = np.eye(3), np.zeros(3), 0.0, 0, 0
        if not same:
            try:
                rot, tr, rms, n, offset = align_with_numbering(query, receiver)
            except ValueError as error:
                exclusions.append([query["pdb"], receiver["pdb"], str(error)])
                continue
        ref = np.array(query["mapped"]) @ rot + tr
        for i, seed in enumerate(args.seeds):
            arms = [("armijo", args.binary), ("analytic", args.binary)]
            if args.finite_binary:
                arms.append(("finite", args.finite_binary))
            for label, binary in arms[i % len(arms) :] + arms[: i % len(arms)]:
                jobs.append(
                    (
                        work,
                        out / "poses",
                        binary.resolve(),
                        query,
                        receiver,
                        ref,
                        rms,
                        n,
                        offset,
                        seed,
                        label,
                        args.max_evals,
                    )
                )
    (out / "alignment_exclusions.json").write_text(
        json.dumps(exclusions, indent=2) + "\n"
    )
    with ThreadPoolExecutor(args.jobs) as pool:
        results = list(pool.map(measure, jobs))
    write_csv(out / "cells.csv", [cell for cell, _ in results])
    write_csv(out / "modes.csv", [m for _, modes in results for m in modes])


if __name__ == "__main__":
    main()

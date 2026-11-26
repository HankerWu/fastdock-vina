"""Run and score fixed-budget paired Late docking from graph-generated query geometries."""

import argparse
import csv
import itertools
import json
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from time import monotonic

import gemmi
import numpy as np

SEEDS = (51913, 51941, 51973)


def write(path, records):
    with path.open("w", newline="") as stream:
        w = csv.DictWriter(stream, fieldnames=records[0])
        w.writeheader()
        w.writerows(records)


def aligned_receivers(a, b):
    # The same canonical residues may have different PDB author numbering.
    left, right = list(a["ca"]), list(b["ca"])
    seq1 = [k.split(":")[-1] for k in left]
    seq2 = [k.split(":")[-1] for k in right]
    spec = gemmi.align_string_sequences(seq1, seq2, []).cigar_str()
    i = j = 0
    ca1, ca2 = {}, {}
    for n, op in re.findall(r"(\d+)([MID])", spec):
        n = int(n)
        if op == "M":
            for z in range(n):
                if seq1[i + z] == seq2[j + z]:
                    ca1[left[i + z]] = a["ca"][left[i + z]]
                    ca2[left[i + z]] = b["ca"][right[j + z]]
        if op in "MI":
            i += n
        if op in "MD":
            j += n
    if i != len(left) or j != len(right):
        raise ValueError("alignment truncated")
    return dict(a, ca=ca1), dict(b, ca=ca2)


def main(repo, work, binary, arm, cap, njobs, seeds, label, cross_only):
    sys.path.insert(0, str(repo.resolve() / "benchmarks/paired_crossdocking"))
    import run as base

    meta = json.loads((work / "metadata.json").read_text())
    output = work / f"runs_{cap}"
    output.mkdir(exist_ok=True)

    def cell(item):
        query, receiver, seed = item
        a, b = meta[query], meta[receiver]
        tag = f"{query}_{receiver}_{seed}_{arm}"
        output_file = output / (tag + ".pdbqt")
        rot, tr, backbone, n = base.alignment(*aligned_receivers(a, b))
        refs = np.asarray(a["mapped"]) @ rot + tr
        cmd = [
            str(binary.resolve()),
            "--receptor",
            str(work / receiver / "receptor.pdbqt"),
            "--ligand",
            str(work / query / "ligand.pdbqt"),
            "--pose_selection",
            "late",
            "--exhaustiveness",
            "4",
            "--max_evals",
            str(cap),
            "--refine_candidates",
            "40",
            "--num_modes",
            "10",
            "--energy_range",
            "100",
            "--cpu",
            "1",
            "--search_counters",
            "--seed",
            str(seed),
            "--out",
            str(output_file),
        ]
        if arm == "balanced":
            cmd += ["--balanced_torsion_search"]
        if arm == "prepared":
            cmd += ["--prepared_torsion_starts"]
        for ax, c in zip("xyz", b["center"]):
            cmd += ["--center_" + ax, str(c), "--size_" + ax, "24"]
        log = output / (tag + ".log")
        started = monotonic()
        with log.open("w") as f:
            subprocess.run(
                cmd, stdout=f, stderr=subprocess.STDOUT, check=True, timeout=600
            )
        duration = monotonic() - started
        costs = []
        for line in log.read_text().splitlines():
            if line.startswith("Search trajectory:"):
                costs.append({k: int(v) for k, v in re.findall(r"(\w+)=(\d+)", line)})
        if len(costs) != 4:
            raise ValueError("wrong search budget")
        rows = []
        for block in output_file.read_text().split("ENDMDL"):
            atoms = {
                int(line[6:11]): base.coords(line)
                for line in block.splitlines()
                if line.startswith(("ATOM  ", "HETATM"))
            }
            if not atoms:
                continue
            xyz = np.array([atoms[i] for i in a["serial"]])
            rms = float(
                np.sqrt(np.mean(np.sum((xyz - refs) ** 2, axis=2), axis=1)).min()
            )
            score = float(
                next(
                    line.split()[3]
                    for line in block.splitlines()
                    if line.startswith("REMARK VINA RESULT:")
                )
            )
            rows.append(
                dict(
                    target=a["target"],
                    ligand=query,
                    receiver=receiver,
                    kind="cognate" if query == receiver else "cross",
                    seed=seed,
                    arm=arm,
                    seconds=duration,
                    cap=cap,
                    rank=len(rows) + 1,
                    score=score,
                    rmsd=rms,
                    ca_rmsd=backbone,
                    n_ca=n,
                    hunt=sum(c["hunt_evals"] for c in costs),
                    second=sum(c["second_evals"] for c in costs),
                )
            )
        if not rows:
            raise ValueError("no output modes")
        return rows

    tasks = [
        (q, r, k)
        for q, r in itertools.product(meta, repeat=2)
        if meta[q]["target"] == meta[r]["target"] and (not cross_only or q != r)
        for k in seeds
    ]
    with ThreadPoolExecutor(njobs) as ex:
        rows = [row for group in ex.map(cell, tasks) for row in group]
    write(work / f"modes-{arm}-{cap}{label}.csv", rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for field in ("repo", "work", "binary"):
        p.add_argument(field, type=Path)
    p.add_argument("arm", choices=["native", "balanced", "prepared"])
    p.add_argument("cap", type=int)
    p.add_argument("--jobs", type=int, default=3)
    p.add_argument("--seeds", type=int, nargs="+", default=SEEDS)
    p.add_argument(
        "--label", default="", help="Suffix for this run table, e.g. -bridge"
    )
    p.add_argument("--cross-only", action="store_true")
    a = p.parse_args()
    main(a.repo, a.work, a.binary, a.arm, a.cap, a.jobs, a.seeds, a.label, a.cross_only)

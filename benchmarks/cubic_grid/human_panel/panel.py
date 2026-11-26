"""Cognate and reciprocal cross-docking on six independent human targets."""

import argparse
import csv
import itertools
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

PANELS = {
    "EPAS1": {"4GHI": "0X3", "5TBM": "79A"},
    "FKBP1A": {"1FKG": "SB3", "1FKH": "SBX"},
    "RBP4": {"1RBP": "RTL", "3FMZ": "2T1"},
    "FABP4": {"2NNQ": "T4B", "3FR2": "8CA"},
    "DYRK1A": {"3ANR": "HRM", "7AKL": "RK5"},
    "JAK2": {"6VS3": "R6V", "6WTN": "RXT"},
}
SEEDS = [14591, 30553, 75431]
CAP = 24000


def prepare_entry(entry, work, base):
    target, pdb, ligand = entry
    directory = work / pdb
    raw = work / "raw"
    directory.mkdir(exist_ok=True, parents=True)
    raw.mkdir(exist_ok=True, parents=True)
    d = {"target": target, "pdb": pdb, "ligand": ligand}
    try:
        pdbfile, ccd = raw / (pdb + ".pdb"), raw / (ligand + ".cif")
        d["pdb_sha256"] = base.fetch(
            pdbfile, f"https://files.rcsb.org/download/{pdb}.pdb"
        )
        d["ccd_sha256"] = base.fetch(
            ccd, f"https://files.rcsb.org/ligands/download/{ligand}.cif"
        )
        lines = pdbfile.read_text().splitlines()
        d["resolution"] = next(
            line[23:30].strip()
            for line in lines
            if line.startswith("REMARK   2 RESOLUTION.")
        )
        records = base.dep_records(lines, ligand)
        observed = np.array(list(records[0].values()))
        d["center"] = observed.mean(axis=0).tolist()
        # Chain A is sufficient only for a ligand isolated from excluded groups.
        omitted = []
        for line in lines:
            if not line.startswith(("HETATM", "ATOM  ")) or line[76:78].strip() in (
                "H",
                "D",
            ):
                continue
            if line[17:20] == "HOH" or (
                line[21] == "A"
                and (line.startswith("ATOM  ") or line[17:20].strip() == ligand)
            ):
                continue
            distance = float(np.linalg.norm(observed - base.coords(line), axis=1).min())
            omitted.append((distance, line[17:27].strip()))
        d["nearest_nonwater_omitted"] = min(omitted) if omitted else None
        if omitted and min(omitted)[0] < 8:
            raise ValueError("chain A omits nonwater group within 8 A")
        d.update(base.ligand(ccd, records, directory / "ligand.pdbqt"))
        # receptor() excludes unresolved loops, rebuilds absent atoms on resolved residues.
        d["ca"] = {
            l[22:27].strip() + ":" + l[17:20]: base.coords(l)
            for l in lines
            if l.startswith("ATOM  ")
            and l[21] == "A"
            and l[12:16].strip() == "CA"
            and l[16] in " A"
        }
        if pdb == "6WTN":
            # GLN1003, ASN874 and multiple other side chains are missing in the
            # deposition. Rebuilding the pocket in successive tests yielded invalid
            # PDBQT topology involving Lys883 (observed ligand distance 6.80 A).
            # Its ligand and receptor frame remain usable for cross-docking.
            d.update(
                status="query_only",
                reason="6WTN Lys883 missing sidechain within 6.81 A, rebuilt receptor fails topology checks",
            )
        else:
            d.update(base.receptor(lines, directory, pdb, observed))
            d["status"] = "ready"
    except (ValueError, RuntimeError, OSError, subprocess.SubprocessError) as exc:
        d.update(status="excluded", error=str(exc))
    (directory / "preparation.json").write_text(json.dumps(d, indent=2) + "\n")
    print(pdb, d["status"], d.get("error", ""), flush=True)
    return d


def run_entry(item, work, binary, base, arm, cap):
    query, receptor, seed = item
    key = f"{query['pdb']}_{receptor['pdb']}_{seed}_{arm}_{cap}"
    runs = work / "runs"
    runs.mkdir(exist_ok=True)
    out = runs / (key + ".pdbqt")
    rot, trans, ca_rmsd, ca_n = base.alignment(query, receptor)
    references = np.stack(query["mapped"]) @ rot + trans
    cmd = [
        str(binary),
        "--receptor",
        str(work / receptor["pdb"] / "receptor.pdbqt"),
        "--ligand",
        str(work / query["pdb"] / "ligand.pdbqt"),
        "--pose_selection",
        "late",
        "--cpu",
        "1",
        "--exhaustiveness",
        "3",
        "--max_evals",
        str(cap),
        "--num_modes",
        "9",
        "--refine_candidates",
        "15",
        "--energy_range",
        "100",
        "--search_counters",
        "--seed",
        str(seed),
        "--out",
        str(out),
    ]
    for axis, position in zip("xyz", receptor["center"]):
        cmd += ["--center_" + axis, str(position), "--size_" + axis, "24"]
    if arm == "cubic":
        cmd.append("--cubic_search_grid")
    start = time.monotonic()
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=300, check=False)
    elapsed = time.monotonic() - start
    (runs / (key + ".log")).write_text(r.stdout + r.stderr)
    r.check_returncode()
    (runs / (key + ".seconds")).write_text(str(elapsed))
    rows = []
    for i, model in enumerate(out.read_text().split("ENDMDL")[:-1], 1):
        ls = model.splitlines()
        xyz = {
            int(l[6:11]): base.coords(l)
            for l in ls
            if l.startswith(("HETATM", "ATOM  "))
        }
        ordered = np.array([xyz[s] for s in query["serial"]])
        rmsd = np.sqrt(((ordered - references) ** 2).sum(-1).mean(-1)).min()
        score = float(
            next(l.split()[3] for l in ls if l.startswith("REMARK VINA RESULT:"))
        )
        rows.append(
            {
                "target": query["target"],
                "query": query["pdb"],
                "receiver": receptor["pdb"],
                "kind": "cognate" if query["pdb"] == receptor["pdb"] else "cross",
                "seed": seed,
                "arm": arm,
                "cap": cap,
                "model": i,
                "score": score,
                "rmsd": rmsd,
                "aligned_ca": ca_n,
                "pocket_rmsd": ca_rmsd,
            }
        )
    if not rows:
        raise ValueError(key)
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=("prepare", "native", "cubic"))
    p.add_argument("--repo", type=Path, required=True)
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--binary", type=Path)
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--cap", type=int, default=CAP)
    a = p.parse_args()
    work = a.work.resolve()
    os.environ.setdefault("OPENMM_CPU_THREADS", "1")
    sys.path.insert(0, str(a.repo.resolve() / "benchmarks/paired_crossdocking"))
    import run as base

    entries = [
        (t, pdb, lig) for t, group in PANELS.items() for pdb, lig in group.items()
    ]
    if a.action == "prepare":
        with ThreadPoolExecutor(a.jobs) as pool:
            list(pool.map(lambda item: prepare_entry(item, work, base), entries))
        return
    todo = []
    for group in PANELS.values():
        ds = [
            json.loads((work / pdb / "preparation.json").read_text()) for pdb in group
        ]
        good = [d for d in ds if d["status"] == "ready"]
        queries = [d for d in ds if d["status"] in ("ready", "query_only")]
        todo.extend(itertools.product(queries, good, SEEDS))
    with ThreadPoolExecutor(a.jobs) as pool:
        jobs = pool.map(
            lambda item: run_entry(
                item, work, a.binary.resolve(), base, a.action, a.cap
            ),
            todo,
        )
        rows = [row for job in jobs for row in job]
    with (work / f"modes_{a.action}_{a.cap}.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()

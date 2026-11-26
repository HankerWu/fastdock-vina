"""Evaluate frozen pose sets. Input directories are those of the anchor research benchmarks."""

import argparse
import csv
import hashlib
import json
import re
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
from polar_endpoint import PolarEndpoint, atoms, read_models


def read(path):
    with path.open() as src:
        return list(csv.DictReader(src))


def write(path, rows):
    with path.open("w", newline="") as fp:
        w = csv.DictWriter(fp, fieldnames=rows[0], lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def jobs(target, runs, prepared, binary):
    if target.startswith("ca_holdout"):
        table = read(runs.parent / "holdout_modes.csv")
        metadata = {}
    elif target == "hcaii":
        table = read(runs / "candidates.csv")
        metadata = json.loads((prepared / "metadata.json").read_text())
    else:
        table = read(runs / "modes.csv")
        metadata = json.loads((prepared / "manifest.json").read_text())
    grouped = defaultdict(list)
    for row in table:
        r, l = row["receptor"], row["ligand"]
        slot = row.get("slot") or row["arm"]
        stem = f"{r}_{row.get('crystal', l)}_{row['seed']}_{slot}"
        exp = (
            {"AZM": 7.2, "MZM": 1.2}[l]
            if target.startswith("ca_holdout")
            else metadata[l]["Ki_nM" if target == "hcaii" else "kd_uM"]
        )
        if target.startswith("ca_holdout") and r == "1JD0":
            exp = ""
        grouped[r, l].append(
            {
                "target": target,
                "receptor": r,
                "ligand": l,
                "seed": int(row["seed"]),
                "slot": slot,
                "model": int(row["model"]),
                "rmsd": float(row["rmsd"]),
                "sign": int(row.get("same_sign_as_bound", 1)),
                "exp": exp,
                "stem": stem,
                "pool_sha256": row.get("pool_sha256", ""),
            }
        )
    return [
        (target, rec, lig, rows, runs, prepared, binary)
        for (rec, lig), rows in sorted(grouped.items())
    ]


def native(binary, prepared, rec, ligand, pose):
    cmd = [
        str(binary),
        "--cpu",
        "1",
        "--seed",
        "1",
        "--receptor",
        str(prepared / f"{rec}_receptor.pdbqt"),
        "--ligand",
        str(ligand),
        "--config",
        str(
            next(
                p
                for p in (prepared / f"{rec}_box.txt", prepared / f"{rec}.box")
                if p.exists()
            )
        ),
        "--pose_coordinates",
        str(pose),
        "--net_polar",
        "native",
        "--score_only",
        "--verbosity",
        "0",
    ]
    out = subprocess.check_output(cmd, text=True)
    vals = lambda label: [
        float(t)
        for t in re.search(r"^Net polar " + label + r": (.+)$", out, re.MULTILINE)[
            1
        ].split()
    ]
    inter, intra, _ = vals("native")
    return {
        "inter": inter,
        "intra": intra,
        "factor": vals("normalization")[0],
        "hbond": vals("terms")[0],
    }


def group(job):
    target, rec, lig, rows, runs, prepared, binary = job
    receptor = atoms(
        (prepared / f"{rec}_receptor.pdbqt").read_text().splitlines(), receptor=True
    )
    ligand_path = prepared / (
        f"{lig}.pdbqt" if target.startswith("ca_holdout") else f"{lig}_ligand.pdbqt"
    )
    calc = PolarEndpoint(receptor, atoms(ligand_path.read_text().splitlines()))
    grouped = defaultdict(list)
    for r in rows:
        grouped[r["stem"]].append(r)
    result = []
    with TemporaryDirectory() as temp:
        pose = Path(temp) / "mode.pdbqt"
        for stem, entries in sorted(grouped.items()):
            path = runs / (stem + ".pdbqt")
            data = path.read_bytes()
            checksum = hashlib.sha256(data).hexdigest()
            models = read_models(data.decode().splitlines())
            for row in entries:
                if row["pool_sha256"] and row["pool_sha256"] != checksum:
                    raise ValueError(f"Changed state pool {path}")
                lines = models[row["model"]]
                pose.write_text("\n".join(lines) + "\n")
                topology = (
                    runs / "inputs" / f"{lig}_reflected.pdbqt"
                    if row["slot"] == "second_reflected"
                    else ligand_path
                )
                terms = native(binary, prepared, rec, topology, pose)
                result.append(
                    {
                        **row,
                        "pool_sha256": checksum,
                        **terms,
                        **calc.evaluate(atoms(lines)),
                    }
                )
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "target", choices=["hcaii", "brd4", "t4l", "ca_holdout36", "ca_holdout20"]
    )
    p.add_argument("runs", type=Path)
    p.add_argument("prepared", type=Path)
    p.add_argument("binary", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--jobs", type=int, default=2)
    a = p.parse_args()
    a.output.parent.mkdir(parents=True, exist_ok=True)
    with ProcessPoolExecutor(a.jobs) as ex:
        batches = ex.map(group, jobs(a.target, a.runs, a.prepared, a.binary.resolve()))
        rows = []
        for batch in batches:
            print(
                a.target,
                batch[0]["receptor"],
                batch[0]["ligand"],
                len(batch),
                flush=True,
            )
            rows.extend(batch)
    write(a.output, rows)

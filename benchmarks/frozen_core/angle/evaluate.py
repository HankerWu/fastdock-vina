"""Rerank native poses using a frozen donor-H angular penalty; evaluate within each pool."""
import argparse
import csv
import hashlib
import io
import json
import subprocess
import sys
import tarfile
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
UPSTREAM = ROOT / "benchmarks/geometry_comparison/upstream"
sys.path.insert(0, str(UPSTREAM))
from verify import references

INPUTS = {
    "cognate": ROOT / "benchmarks/geometry_comparison/pose/inputs.tar.xz",
    "reciprocal": ROOT / "benchmarks/cubic_grid/human_panel/prepared.tar.xz",
}


def write(path, rows):
    with path.open("w", newline="") as stream:
        wr = csv.DictWriter(stream, fieldnames=rows[0])
        wr.writeheader()
        wr.writerows(rows)


def run(binary, output):
    cells, poses = [], []
    output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        for protocol, source in INPUTS.items():
            inputs = work / protocol
            inputs.mkdir()
            with tarfile.open(source) as tar:
                tar.extractall(inputs, filter="data")
            for design in ("fixed", "matched", "per_complex"):
                study = UPSTREAM / protocol
                if design != "fixed":
                    study /= design
                    caps = json.loads((study / "frozen_cap.json").read_text())
                with (study / "cells.csv").open() as stream:
                    runs = [r for r in csv.DictReader(stream) if r["arm"] == "native"]
                with tarfile.open(study / "raw.tar.xz") as pack:
                    for row in runs:
                        q, r, seed = row["query"], row["receiver"], row["seed"]
                        cap = 12000 if design == "fixed" else caps["cap"] if design == "matched" else next(
                            x["cap"] for x in caps if x["query"] == q and x["receiver"] == r)
                        base = dict(protocol=protocol, design=design, cap=cap, target=row["target"],
                                    query=q, receiver=r, seed=seed)
                        filename = f"{q}-{r}-{seed}-native.pdbqt"
                        pool = work / filename
                        pool.write_bytes(pack.extractfile(filename).read())
                        assert hashlib.sha256(pool.read_bytes()).hexdigest() == row["pool_sha256"]
                        cmd = [binary, inputs/r/"receptor.pdbqt", inputs/q/"ligand.pdbqt", pool]
                        cp = subprocess.run(cmd, capture_output=True, text=True)
                        if cp.returncode:
                            raise RuntimeError(f"{q} -> {r}. {cp.stderr}")
                        terms = list(csv.DictReader(io.StringIO(cp.stdout)))
                        meta = "metadata.json" if protocol == "cognate" else "preparation.json"
                        query, rec = [json.loads((inputs/z/meta).read_text()) for z in (q,r)]
                        ref = references(query, rec)
                        conformers = [b for b in pool.read_text().split("ENDMDL") if "MODEL" in b]
                        assert len(conformers) == len(terms) == int(row["modes"])
                        data = []
                        for rank, (part, term) in enumerate(zip(conformers, terms), 1):
                            atom_xyz = {int(l[6:11]): [float(l[k:k+8]) for k in (30,38,46)]
                                        for l in part.splitlines() if l.startswith(("ATOM  ", "HETATM"))}
                            xyz = np.array([atom_xyz[i] for i in query["serial"]])
                            rms = float(np.sqrt(((ref - xyz)**2).sum(axis=2).mean(axis=1)).min())
                            score = float(next(l for l in part.splitlines() if l.startswith("REMARK VINA RESULT:")).split()[3])
                            assert int(term["rank"]) == rank
                            penalty = float(term["score_delta"])
                            assert float(term["native_score"]) == score
                            corrected = float(term["corrected_score"])
                            assert abs(corrected - score - penalty) < 1e-9
                            entry = {**base, **term, "rmsd_A": rms, "native_score": score, "corrected_score": corrected}
                            data.append(entry)
                        assert abs(data[0]["rmsd_A"]-float(row["first"])) < 1e-9
                        assert abs(min(d["rmsd_A"] for d in data)-float(row["best"])) < 1e-9
                        ordered = sorted(data, key=lambda d: (d["corrected_score"], int(d["rank"])))
                        cells.append(dict(**base, candidates=len(data), native_top=data[0]["rmsd_A"],
                                          native_top5=min(d["rmsd_A"] for d in data[:5]),
                                          angle_top=ordered[0]["rmsd_A"], angle_rank=ordered[0]["rank"],
                                          angle_top5=min(d["rmsd_A"] for d in ordered[:5]),
                                          oracle=min(d["rmsd_A"] for d in data)))
                        poses.extend(data)
    subsets = defaultdict(list)
    for cell in cells:
        subsets[cell["protocol"], cell["design"]].append(cell)
    summary = []
    for (protocol, design), group in subsets.items():
        summary.append(dict(protocol=protocol, design=design, n=len(group),
                            **{k: sum(r[k] <= 2 for r in group) for k in
                               ("native_top","angle_top","native_top5","angle_top5","oracle")},
                            rescues=sum(r["native_top"] > 2 >= r["angle_top"] for r in group),
                            losses=sum(r["angle_top"] > 2 >= r["native_top"] for r in group),
                            switches=sum(int(r["angle_rank"]) != 1 for r in group)))
    for name, rows in (("poses",poses),("cells",cells),("summary",summary)):
        write(output / f"{name}.csv",rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("binary",type=Path)
    p.add_argument("output",type=Path)
    a = p.parse_args()
    run(a.binary.resolve(),a.output)

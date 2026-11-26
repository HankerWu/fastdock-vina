"""Independently verify the archived PDBQT, final ranks and RMSD hit counts."""
import argparse
import csv
import hashlib
import json
import math
import tarfile
from collections import defaultdict
from pathlib import Path

import numpy as np


def references(query, receiver):
    mappings = np.array(query["mapped"])
    if query["pdb"] == receiver["pdb"]:
        return mappings
    # Row-vector Kabsch reconstruction of the original quaternion-fit reference.
    keys = [k for k in sorted(set(query["ca"]) & set(receiver["ca"]))
            if all(math.dist(entry["ca"][k], entry["center"]) < 14 for entry in (query, receiver))]
    if len(keys) < 12:
        raise ValueError("Insufficient matched C-alpha residues")
    a, b = (np.array([entry["ca"][k] for k in keys]) for entry in (query, receiver))
    u, _, vt = np.linalg.svd((a - a.mean(0)).T @ (b - b.mean(0)))
    rotation = u @ np.diag([1, 1, np.linalg.det(u @ vt)]) @ vt
    return (mappings - a.mean(0)) @ rotation + b.mean(0)


def audit(results, fixtures, protocol):
    runs = list(csv.DictReader((results / "cells.csv").open()))
    expected_by_pair = {}
    same = defaultdict(dict)
    max_error, count = 0.0, 0
    summary = defaultdict(lambda: [0, 0, 0, 0, 0])
    with tarfile.open(results / "raw.tar.xz") as archive:
        for r in runs:
            q, p, seed, arm = r["query"], r["receiver"], r["seed"], r["arm"]
            if (q, p) not in expected_by_pair:
                filename = "metadata.json" if protocol == "cognate" else "preparation.json"
                query, receiver = [json.loads((fixtures / code / filename).read_text()) for code in (q, p)]
                expected_by_pair[q, p] = query["serial"], references(query, receiver)
            serial, ref = expected_by_pair[q, p]
            pool = archive.extractfile(f"{q}-{p}-{seed}-{arm}.pdbqt").read()
            if hashlib.sha256(pool).hexdigest() != r["pool_sha256"]:
                raise ValueError("Output checksum differs")
            rmsds, scores = [], []
            for model in pool.decode().split("ENDMDL"):
                pos = {int(l[6:11]): tuple(map(float, [l[30:38], l[38:46], l[46:54]]))
                       for l in model.splitlines() if l.startswith(("ATOM  ", "HETATM"))}
                if not pos:
                    continue
                scores.append(float(next(l for l in model.splitlines() if l.startswith("REMARK VINA RESULT:")).split()[3]))
                xyz = [pos[i] for i in serial]
                rmsds.append(min(math.sqrt(math.fsum(math.dist(x, y)**2 for x, y in zip(xyz, mapping)) / len(serial))
                                 for mapping in ref))
            obs = (rmsds[0], min(rmsds[:5]), min(rmsds))
            max_error = max(max_error, *(abs(v - float(r[k])) for v, k in zip(obs, ("first", "best5", "best"))))
            if len(rmsds) != int(r["modes"]) or scores[0] != float(r["native_score"]):
                raise ValueError("Mismatch between printed conformers, order or reported top score")
            same[q, p, seed][arm] = r
            for j, v in enumerate(obs):
                summary[arm][j] += int(v <= 2)
            summary[arm][3] += 1
            summary[arm][4] += len(rmsds)
            count += len(rmsds)
    if max_error > 1e-9:
        raise ValueError("Independent RMSD mismatch")
    stored = json.loads((results / "summary.json").read_text())
    for arm, row in stored.items():
        models = row.get("model_count", row.get("models"))
        if summary[arm] != [row[k] for k in ("top1", "top5", "oracle", "n")] + [models]:
            raise ValueError(f"Aggregate mismatch in {arm}")
    diffs = []
    for group in same.values():
        if set(group) != set(stored):
            raise ValueError("Incomplete matched arms")
        for field in ("hunt_evals", "second_evals"):
            if len({r[field] for name, r in group.items() if name != "official"}) != 1:
                raise ValueError("Global-search evaluation counts differ")
        if "geometry_fd" in group:
            diffs.append({f: abs(float(group["geometry_fd"][f]) - float(group["geometry_analytic"][f]))
                          for f in ("first", "best", "native_score")})
    def num_equal(a, b, fields):
        return sum(all(d[a][v] == d[b][v] for v in fields) for d in same.values())
    return dict(runs=len(runs), models=count, cells=len(same), max_independent_rmsd_error_A=max_error,
                official_native_byte_identical=num_equal("official", "native", ["pool_sha256"]),
                official_native_equal_top1_rmsd_and_score=num_equal("official", "native", ["first", "native_score"]),
                geometry_fd_analytic_byte_identical=num_equal("geometry_fd", "geometry_analytic", ["pool_sha256"]) if diffs else None,
                geometry_fd_analytic_max_abs={f:max(d[f] for d in diffs) for f in diffs[0]} if diffs else None,
                fork_hunt_and_second_evaluations_equal_in_all_cells=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("protocol", choices=("cognate", "reciprocal"))
    p.add_argument("fixtures", type=Path)
    p.add_argument("results", type=Path)
    a = p.parse_args()
    print(json.dumps(audit(a.results, a.fixtures, a.protocol), indent=2))

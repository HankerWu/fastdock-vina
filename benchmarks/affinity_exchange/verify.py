"""Independent replay from the published per-ligand inputs, without importing analyze."""

import argparse
import csv
from collections import defaultdict
import hashlib
import itertools
import json
import lzma
import math
from pathlib import Path

import numpy as np

RT = 0.00198720425864083 * 298.15
BENCH = Path(__file__).resolve().parents[1]


def read_csv(path):
    opener = lzma.open if path.suffix == ".xz" else open
    with opener(path, "rt") as handle:
        return list(csv.DictReader(handle))


def close(a, b, tolerance=2e-7):
    if not math.isclose(float(a), float(b), rel_tol=0, abs_tol=tolerance):
        raise ValueError(f"{a} differs from {b}")


def verify(out):
    for path, digest in json.loads((out / "inputs.json").read_text()).items():
        assert hashlib.sha256((BENCH / path).read_bytes()).hexdigest() == digest
    arch = BENCH / "affinity_reference_ensemble/results"
    sources = {}
    for count, name in [("1024", "primary"), ("8192", "sensitivity")]:
        sources[count] = {
            (r["target"], r["ligand"]): r
            for r in read_csv(arch / f"{name}_ligands.csv")
        }
    inputs = read_csv(out / "contrasts.csv.xz")
    groups = defaultdict(list)
    for row in inputs:
        arm, i, j = (row[k] for k in ("arm", "i", "j"))
        t, u = row["target_left"], row["target_right"]
        per_ligand = sources[row["nodes"]]
        # Compute target selectivity first, whereas the main program takes
        # within-target differences first. Add arbitrary target and ligand zeros.
        source = [per_ligand[k] for k in [(t, i), (u, i), (t, j), (u, j)]]
        v = [float(r[arm]) for r in source]
        shifted = [
            v[n] + (2.75 if n % 2 else -3.5) + (1.375 if n < 2 else -2.25)
            for n in range(4)
        ]
        predicted = (v[0] - v[1]) - (v[2] - v[3])
        close(predicted, (shifted[0] - shifted[1]) - (shifted[2] - shifted[3]), 1e-12)
        close(predicted, row["predicted"])
        kd = [float(r["kd_nM"]) for r in source]
        observed = RT * (math.log(kd[0] / kd[1]) - math.log(kd[2] / kd[3]))
        close(observed, row["experimental"])
        if arm == "size":
            close(predicted, 0, 0)
        if row["reversal_ge3"] == "True":
            # Score_ti - score_tj + delta must be strictly below (or above)
            # zero on each target; test intersection of the two open half-lines.
            lower, upper = -float("inf"), float("inf")
            for x, y in [(0, 2), (1, 3)]:
                cutoff = v[y] - v[x]
                if kd[x] < kd[y]:
                    upper = min(upper, cutoff)
                else:
                    lower = max(lower, cutoff)
            assert (lower >= upper) == (float(row["predicted"]) * observed <= 0)
        key = (row["nodes"], row["panel"], t, u, arm)
        groups[key].append(row)
    for r in read_csv(out / "summary.csv"):
        key = tuple(r[k] for k in ("nodes", "panel", "left", "right", "arm"))
        rows = groups.pop(key)
        ids = {value for row in rows for value in (row["i"], row["j"])}
        assert len(rows) == math.comb(len(ids), 2) == int(r["n"])
        errors = [float(v["predicted"]) - float(v["experimental"]) for v in rows]
        close(sum(map(abs, errors)) / len(errors), r["mae"])
        close(math.sqrt(sum(v * v for v in errors) / len(errors)), r["rmse"])
    assert not groups

    # Eliminate free ligand zeros analytically by centering observations within
    # ligand. Fit only the 2 or 3 target constants, then center again.
    for row in read_csv(out / "reference_floor.csv"):
        source = sources[row["nodes"]]
        panel = row["panel"]
        targets = ("AURKA", "CHEK1", "S6K1") if panel == "kinase" else ("Gal1", "Gal3")
        blocks = defaultdict(list)
        for (t, i), item in source.items():
            if t not in targets or (panel == "galectin_5_16" and i == "17"):
                continue
            y = RT * math.log(float(item["kd_nM"]) * 1e-9) - float(item[row["arm"]])
            blocks[i].append([y, *(float(t == name) for name in targets)])
        centered = np.vstack(
            [np.asarray(v) - np.mean(v, axis=0) for v in blocks.values()]
        )
        y, x = centered[:, 0], centered[:, 1:]
        offsets, *_ = np.linalg.lstsq(x, y, rcond=None)
        residual = y - x @ offsets
        close(np.sqrt(np.mean(residual**2)), row["best_additive_cell_rmse"])
    return len(inputs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path)
    result = verify(parser.parse_args().results)
    print(f"{result} contrast rows and all additive-reference floors verified")

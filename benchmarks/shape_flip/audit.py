"""Recompute reported endpoints from raw files, independent of the docking executables."""

import csv
import hashlib
import itertools
import json
import re
import tarfile
from decimal import Decimal
from pathlib import Path

HERE = Path(__file__).resolve().parent


def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def main():
    dev = HERE / "development"
    assay = HERE / "gal3"
    dev_rows = rows(dev / "structures.csv")
    modes = rows(dev / "models.csv")
    with tarfile.open(dev / "raw.tar.xz") as data:
        for row in dev_rows:
            raw = data.extractfile("raw/" + row["key"] + ".pdbqt").read()
            assert hashlib.sha256(raw).hexdigest() == row["sha256"]
            subset = [m for m in modes if m["key"] == row["key"]]
            assert len(subset) == raw.count(b"ENDMDL") == int(row["modes"])
            rr = [float(m["rmsd_A"]) for m in subset]
            assert (rr[0], min(rr[:5]), min(rr)) == tuple(
                float(row[k]) for k in ("top1_A", "best5_A", "oracle_A")
            )
            assert int(row["joint"]) == sum(
                int(row[k]) for k in ("hunt_evals", "second_evals", "explicit_evals")
            )
    labels = {
        r["compound"]: Decimal(r["k"])
        for r in rows(HERE.parents[0] / "accepted_retention/affinity/inputs.csv")
        if r["panel"] == "Gal3"
    }
    score = {}
    with tarfile.open(assay / "raw.tar.xz") as data:
        for row in rows(assay / "runs.csv"):
            raw = data.extractfile("raw/" + row["key"] + ".pdbqt").read()
            scores = re.findall(rb"^REMARK VINA RESULT:\s+(\S+)", raw, re.MULTILINE)
            assert len(scores) == int(row["modes"])
            first = Decimal(scores[0].decode())
            assert first == Decimal(row["first_score"])
            score[row["arm"], row["receptor"], row["compound"]] = first
    ids = sorted(labels, key=labels.get)
    reliable = [
        (i, j) for i, j in itertools.combinations(ids, 2) if labels[j] >= 3 * labels[i]
    ]
    assert len(score) == 376 and len(reliable) == 403
    credits = {}
    for arm in ("parent", "always", "adaptive", "upstream"):
        ss = {c: sum(score[arm, rec, c] for rec in ("1KJL", "3T1L")) for c in ids}
        credits[arm] = sum(
            (ss[i] < ss[j]) + 0.5 * (ss[i] == ss[j]) for i, j in reliable
        )
        record = next(
            r
            for r in rows(assay / "ranks.csv")
            if (r["arm"], r["receptor"], r["min_ratio"], r["match_charge_size"])
            == (arm, "mean", "3", "False")
        )
        assert credits[arm] == float(record["credits"])
    # Compound-level deletion keeps dependence among all its incident edges.
    sensitivity = {}
    for arm in ("always", "adaptive"):

        def credit(a, i, j):
            ss = [sum(score[a, rec, c] for rec in ("1KJL", "3T1L")) for c in (i, j)]
            return (ss[0] < ss[1]) + 0.5 * (ss[0] == ss[1])

        contrasts = [
            sum(
                credit(arm, i, j) - credit("parent", i, j)
                for i, j in reliable
                if c not in (i, j)
            )
            for c in ids
        ]
        sensitivity[arm] = [min(contrasts), max(contrasts)]
    print(
        json.dumps(
            {
                "development_conditions": len(dev_rows),
                "development_poses": len(modes),
                "fixed_development_modes": len(rows(dev / "fixed_pools.csv")),
                "gal3_runs": len(score),
                "gal3_fixed": len(rows(assay / "fixed.csv")),
                "gal3_assay_pairs": len(reliable),
                "gal3_consensus": credits,
                "leave_one_ligand_out_deltas": sensitivity,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

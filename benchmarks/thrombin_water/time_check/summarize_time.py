"""Relate measured whole-process elapsed/CPU time to crystal recovery in the confirmation seeds."""

import argparse
import csv
import hashlib
import itertools
import json
import os
import sys
from pathlib import Path
from statistics import median

import numpy as np


# At install time this directory is benchmarks/thrombin_water/time_check or soft_search/time_check.
def audit_dir():
    for p in Path(__file__).resolve().parents:
        if (p / "independent_audit.py").is_file():
            return p
    raise FileNotFoundError("Run alongside the parent thrombin benchmark")


sys.path.insert(0, str(audit_dir()))
from independent_audit import get_serials, positions, read_csv, refs_from_ccd


def save(path, rows):
    with path.open("w") as h:
        w = csv.DictWriter(h, fieldnames=rows[0])
        w.writeheader()
        w.writerows(rows)


def main(inp, work, old=None):
    out = work.resolve()
    os.chdir(inp.resolve())
    rows = read_csv(out / "validation.csv")
    series = {r["id"]: r for r in read_csv("series.csv")}
    ref_cache = {}
    evaluated = []
    identical = 0
    for row in rows:
        rec, lig, seed, method = (
            row[k] for k in ("receiver", "ligand", "seed", "phase")
        )
        exp = series[lig]
        if (rec, lig) not in ref_cache:
            top = Path("prepared", lig + ".pdbqt").read_text().splitlines()
            refs = refs_from_ccd(exp["PDB"], exp["CCD"], top, rec)
            ref_cache[rec, lig] = (refs, get_serials(top))
        refs, serials = ref_cache[rec, lig]
        path = out / (row["receipt"] + ".pdbqt")
        rmsds = []
        if old and method in ("base", "soft"):
            prev = (
                old
                / f"{rec}_{lig}_{seed}_{'hard' if method == 'base' else 'soft'}.pdbqt"
            )
            assert hashlib.sha256(prev.read_bytes()).hexdigest() == row["output_sha256"]
            identical += 1
        for m in path.read_text().split("ENDMDL"):
            if not m.strip():
                continue
            atoms = {
                int(l[6:11]): positions(l)
                for l in m.splitlines()
                if l.startswith(("ATOM", "HETATM"))
            }
            xyz = np.array([atoms[i] for i in serials])
            rmsds.append(
                min(np.sqrt(np.mean(np.sum((xyz - rr) ** 2, axis=1))) for rr in refs)
            )
        curr = dict(
            **row,
            kind="cognate" if rec == exp["PDB"] else "cross",
            top=rmsds[0],
            best=min(rmsds),
            n_models=len(rmsds),
        )
        curr["wall_seconds"] = float(row["wall_seconds"])
        curr["child_cpu_seconds"] = float(row["child_cpu_seconds"])
        evaluated.append(curr)
    save(out / "measured_cells.csv", evaluated)
    summary = []
    for kind, method in itertools.product(
        ("all", "cross"), ("base", "soft", "hard_equal")
    ):
        items = [
            c
            for c in evaluated
            if c["phase"] == method and (kind == "all" or kind == c["kind"])
        ]
        summary.append(
            {
                "kind": kind,
                "method": method,
                "n": len(items),
                "cap": items[0]["cap"],
                "top": sum(c["top"] <= 2 for c in items),
                "oracle": sum(c["best"] <= 2 for c in items),
                "total_wall": sum(c["wall_seconds"] for c in items),
                "total_cpu": sum(c["child_cpu_seconds"] for c in items),
                "median_wall": median(c["wall_seconds"] for c in items),
            }
        )
    save(out / "measured_summary.csv", summary)
    diff = []
    for comparison in ("base", "hard_equal"):
        subset = {
            (c["receiver"], c["ligand"], c["seed"]): c
            for c in evaluated
            if c["phase"] == comparison and c["kind"] == "cross"
        }
        soft = {
            (c["receiver"], c["ligand"], c["seed"]): c
            for c in evaluated
            if c["phase"] == "soft" and c["kind"] == "cross"
        }
        gain = sum(c["top"] > 2 and soft[k]["top"] <= 2 for k, c in subset.items())
        loss = sum(c["top"] <= 2 and soft[k]["top"] > 2 for k, c in subset.items())
        tt = [c["wall_seconds"] / soft[k]["wall_seconds"] for k, c in subset.items()]
        diff.append(
            {
                "hard_policy": comparison,
                "soft_gains": gain,
                "soft_losses": loss,
                "n": len(tt),
                "ratio_p25": np.quantile(tt, 0.25),
                "ratio_median": np.median(tt),
                "ratio_p75": np.quantile(tt, 0.75),
            }
        )
    save(out / "measured_switches.csv", diff)
    print(
        json.dumps(
            {"summary": summary, "switched": diff, "previous_sha256_equal": identical},
            indent=2,
            default=lambda a: a.item(),
        )
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("inputs", type=Path)
    p.add_argument("work", type=Path)
    p.add_argument("--previous-confirm", type=Path)
    a = p.parse_args()
    main(a.inputs, a.work, a.previous_confirm.resolve() if a.previous_confirm else None)

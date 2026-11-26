"""Unique within-table, within-target Ki/Kd edges after seed aggregation."""

import argparse
import csv
import itertools
import statistics
from pathlib import Path


def summarize(folder):
    results = []
    edges = []
    for path in sorted(folder.glob("*_rows.csv")):
        rows = list(csv.DictReader(path.open()))
        for target, assay, endpoint in sorted(
            {(r["target"], r["assay"], r["type"]) for r in rows}
        ):
            selected = [
                r
                for r in rows
                if (r["target"], r["assay"], r["type"]) == (target, assay, endpoint)
            ]
            groups = {}
            for r in selected:
                if r["status"] == "ok":
                    groups.setdefault((r["ligand"], r["arm"]), {})[r["seed"]] = float(
                        r["score"]
                    )
            seeds = {r["seed"] for r in selected}
            ki = {
                r["ligand"]: float(r["ki_M"])
                for r in selected
                if all(
                    set(groups.get((r["ligand"], a), {})) == seeds
                    for a in ("native", "soft", "mixed")
                )
            }
            pairs = [
                (u, v)
                for u, v in itertools.combinations(sorted(ki, key=ki.get), 2)
                if ki[v] / ki[u] >= 3
            ]
            for arm in ("native", "soft", "mixed"):
                correct = ties = 0
                for u, v in pairs:
                    d = statistics.median(groups[v, arm].values()) - statistics.median(
                        groups[u, arm].values()
                    )
                    correct += d > 0
                    ties += d == 0
                    edges.append(
                        {
                            "target": target,
                            "assay": assay,
                            "endpoint": endpoint,
                            "arm": arm,
                            "strong": u,
                            "weak": v,
                            "fold": ki[v] / ki[u],
                            "delta_score": d,
                        }
                    )
                results.append(
                    {
                        "target": target,
                        "assay": assay,
                        "endpoint": endpoint,
                        "arm": arm,
                        "n_paired_ligands": len(ki),
                        "n_seeds": len(seeds),
                        "pairs": len(pairs),
                        "correct": correct,
                        "ties": ties,
                    }
                )
    for name, data in [("summary", results), ("edges", edges)]:
        with (folder / f"{name}.csv").open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=data[0])
            writer.writeheader()
            writer.writerows(data)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("folder", type=Path)
    summarize(p.parse_args().folder)

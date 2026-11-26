"""Reconstruct the cyclic fold comparison independently from the 637 candidate rows."""
import argparse
import csv
import statistics
from pathlib import Path

SEEDS = (95701, 95713, 95717)


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def audit(work):
    poses = read(work / "energies.csv")
    folds = read(work / "folds.csv")
    stats = read(work / "summary.csv")
    assert len(poses) == 637 and len(folds) == 216
    for f in folds:
        source_arms = f["plan"].split("+")
        i = SEEDS.index(int(f["first_seed"]))
        seeds = (SEEDS[i], SEEDS[(i + 1) % 3])
        assert seeds[1] == int(f["second_seed"])
        combined = []
        for a, k in zip(source_arms, seeds):
            part = [r for r in poses if (r["target"], r["query"], r["receiver"]) ==
                    (f["target"], f["query"], f["receiver"]) and int(r["seed"]) == k and r["source_arm"] == a]
            assert 0 < len(part) <= 9
            combined += sorted(part, key=lambda r: int(r["rank"]))
        field = f["scoring"] + "_common_score"
        winner = sorted(combined, key=lambda r: float(r[field]))[0]
        oracle = min(float(r["rmsd_A"]) for r in combined)
        assert len(combined) == int(f["candidates"])
        assert (winner["source_arm"], winner["rank"], winner["seed"]) == (
            f["top_source"], f["top_source_rank"], f["top_source_seed"])
        assert winner[field] == f["score"] and float(winner["rmsd_A"]) == float(f["top1_A"])
        assert abs(float(f["oracle_A"]) - oracle) < 1e-12
    for s in stats:
        sample = [f for f in folds if all(f[k] == s[k] for k in ("target", "kind", "scoring", "plan"))]
        assert len(sample) == int(s["folds"])
        for value, key in (("oracle_hits", "oracle_A"), ("top_hits", "top1_A")):
            assert sum(float(f[key]) <= 2 for f in sample) == int(s[value])
        for k in ("source_wall_s", "explicit_evals", "search_evals", "candidates"):
            assert abs(statistics.mean(float(r[k]) for r in sample) - float(s["mean_" + k])) < 1e-9
    print(f"Verified {len(poses)} candidates, {len(folds)} selections and {len(stats)} summary strata.")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    audit(p.parse_args().work)

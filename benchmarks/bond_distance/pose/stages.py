"""Replay the passive search/refinement stage observer on the three affected structures."""

import argparse
import csv
import json
import sys
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("repo", "work", "parent", "fixed"):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    bench = a.repo.resolve() / "benchmarks"
    sys.path.insert(0, str(bench / "pair_gradient/holdout"))
    sys.path.insert(0, str(bench / "pair_gradient/pose"))
    import paired
    from run_holdout import stage_data

    w = a.work.resolve()
    dest = w / "trace"
    (dest / "runs").mkdir(parents=True, exist_ok=True)
    data = {
        p.parent.name: json.loads(p.read_text()) for p in w.glob("*/preparation.json")
    }
    rows, stages, results = [], [], []
    for filename in ("cells.csv", "confirmation.csv"):
        with (w / filename).open() as f:
            rows.extend(
                row
                for row in csv.DictReader(f)
                if row["arm"] in ("parent", "fixed")
                and data[row["query"]]["ligand"] == "3JW"
            )
    for old in rows:
        q, r, key, arm = old["query"], old["receiver"], old["key"], old["arm"]
        for code in (q, r):
            if not (dest / code).exists():
                (dest / code).symlink_to(w / code, target_is_directory=True)
        rr, _ = paired.execute(
            dest,
            data,
            q,
            r,
            int(old["seed"]),
            int(old["cap"]),
            arm,
            getattr(a, arm),
            False,
        )
        if (dest / "runs" / (key + ".pdbqt")).read_bytes() != (
            w / "runs" / (key + ".pdbqt")
        ).read_bytes():
            raise ValueError("Observer changed output " + key)
        for s in ("hunt", "second", "explicit"):
            for field in ("evals", "calls"):
                k = s + "_" + field
                if int(rr[k]) != int(old[k]):
                    raise ValueError("Observer changed work " + key)
        items = stage_data(
            (dest / "runs" / (key + ".log")).read_text(), key, data[q], data[r]
        )
        stages.extend(items)
        results.append(
            {
                "key": key,
                "target": data[q]["target"],
                "arm": arm,
                **{
                    s: min(v["rmsd_A"] for v in items if v["stage"] == s)
                    for s in ("search", "before", "after")
                },
                "returned": old["oracle_A"],
                "top1": old["top1_A"],
            }
        )
    paired.write(w / "stage_modes.csv", stages)
    paired.write(w / "stage_oracles.csv", results)


if __name__ == "__main__":
    main()

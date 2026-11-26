"""Independently recompute the two affinity series from the saved first MODEL remarks."""

import argparse
import csv
import itertools
from decimal import Decimal
from pathlib import Path
from statistics import median

HERE = Path(__file__).resolve().parent
SEEDS = (271828, 314159, 161803, 141421)


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def verify(base, work, results):
    fx = read(HERE.parent / "affinity_partition/fxa_series.csv")
    labels = {
        "fxa": {r["ligand"]: Decimal(r["ki_nm"]) for r in fx if r["included"] == "1"},
        "trypsin": {
            r["ligand"]: Decimal(r["kd_uM"])
            for r in read(base / "fdv-ionpair/trypsin/series.csv")
        },
    }
    rows = [
        r
        for r in read(results / "affinity.csv")
        if r["receptor"] == "all" and r["reference"] == "reported"
    ]
    tested = 0
    for row in rows:
        panel, arm = row["panel"], row["arm"]
        is_try = panel == "trypsin"
        data = labels["trypsin" if is_try else "fxa"]
        seeds = SEEDS[:2] if panel == "inverted" else SEEDS
        recs = ("3PTB", "1BTY") if is_try else ("1FJS", "2P16")
        scores = {}
        for rec, seed, lig in itertools.product(recs, seeds, data):
            name = f"{rec}_{lig}_{seed}"
            original = base / "fdv-ionpair"
            if arm.startswith("mixed"):
                path = work / panel / (name + "_mixed.pdbqt")
            elif arm == "native" and not is_try and seed in SEEDS[:2]:
                path = (
                    base
                    / ("fdv-fxa1" if rec == "1FJS" else "fdv-fxa")
                    / f"{lig}_{seed}.pdbqt"
                )
            else:
                folder = panel if arm != "native" or panel != "inverted" else "fxa"
                path = original / folder / f"{name}_{arm}.pdbqt"
            with path.open() as f:
                assert next(f).split() == ["MODEL", "1"]
                score = next(f).split()
            assert score[:3] == ["REMARK", "VINA", "RESULT:"]
            scores[lig, rec, seed] = Decimal(score[3])
            tested += 1
        n, wins = 0, Decimal(0)
        for a, b in itertools.permutations(data, 2):
            if data[b] < 3 * data[a]:
                continue
            n += 1
            diff = median(
                scores[b, r, s] - scores[a, r, s]
                for r, s in itertools.product(recs, seeds)
            )
            wins += 1 if diff > 0 else Decimal(".5") if diff == 0 else 0
        assert (wins, n) == (Decimal(row["correct"]), int(row["pairs"])), (row, wins)
    print(
        f"Raw first MODELs independently reproduce all {len(rows)} reported chemical pair scores ({tested} file reads)."
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("base", "work", "results"):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    verify(a.base, a.work, a.results)

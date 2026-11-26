"""Collapse repeated receptor/seed observations to one decision per measured ligand pair."""

import collections
import csv
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def main():
    grouped = collections.defaultdict(lambda: collections.defaultdict(list))
    for target, path, gap_field in (
        ("T4L", ROOT / "transfer_holdout_t4l/pairs.csv", "score_gap"),
        ("BRD4", ROOT / "transfer_holdout_brd4/pairs.csv", "score_gap"),
        ("hCAII", ROOT / "hcaii/scaffold/pairs.csv", "score_weak_minus_strong"),
    ):
        for row in csv.DictReader(path.open()):
            if row.get("selection", "native") != "native":
                continue
            arm = "native" if row["arm"] == "native" else "scaffold"
            compound_pair = row["strong"], row["weak"]
            grouped[(target, *compound_pair)][arm].append(float(row[gap_field]))
    out = []
    for (target, strong, weak), arms in sorted(grouped.items()):
        row = {"target": target, "strong": strong, "weak": weak}
        for arm in ("native", "scaffold"):
            gaps = arms[arm]
            delta = statistics.median(gaps)
            row.update(
                {
                    arm + "_n": len(gaps),
                    arm + "_median": round(delta, 6),
                    arm + "_correct": (int(delta > 0) + 0.5 * (delta == 0)),
                    arm + "_min": round(min(gaps), 6),
                    arm + "_max": round(max(gaps), 6),
                }
            )
        values = arms["native"] + arms["scaffold"]
        row["same_sign_all"] = min(values) > 0 or max(values) < 0
        out.append(row)
    with (ROOT / "compound_pairs.csv").open("w") as f:
        writer = csv.DictWriter(f, fieldnames=out[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(out)
    for target in sorted({r["target"] for r in out}):
        rows = [r for r in out if r["target"] == target]
        agree = [r for r in rows if r["same_sign_all"]]
        print(
            target,
            len(rows),
            sum(r["native_correct"] for r in rows),
            sum(r["scaffold_correct"] for r in rows),
            len(agree),
            sum(r["native_correct"] for r in agree),
        )


if __name__ == "__main__":
    main()

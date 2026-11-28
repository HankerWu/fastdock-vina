"""Append a four-version measurement snapshot, or check the archived rows."""

import argparse
import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
POSE = HERE.parent / "geometry_comparison/upstream/four_versions"
LEDGER = HERE / "ledger.csv"
FIELDS = "round_source role engine_source task cohort target assay receptor planned observed mean_top1_A median_top1_A p95_A mean_best5_A mean_oracle_A mean_regret_A rho tau_b correct_3x pairs_3x cpu_s wall_s evidence_file".split()


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def snapshot():
    builds = json.loads((POSE / "builds.json").read_text())["versions"]
    round_source = builds["candidate"]["revision"]
    rows = []

    def add(arm, path, **values):
        row = dict(
            round_source=round_source,
            role=arm,
            engine_source=builds[arm]["revision"],
            evidence_file=str(path.relative_to(HERE.parent.parent)),
            **values,
        )
        rows.append({k: str(row.get(k, "")) for k in FIELDS})

    for path in (POSE / "results/summary.csv", HERE / "external/results/summary.csv"):
        for r in read(path):
            if r["target"] == "ALL":
                add(
                    r["arm"],
                    path,
                    task="pose",
                    cohort=r["cohort"],
                    target="ALL",
                    planned=r["scheduled"],
                    observed=r["valid"],
                    mean_top1_A=r["mean_first"],
                    median_top1_A=r["median"],
                    p95_A=r["p95"],
                    mean_best5_A=r["mean_best5"],
                    mean_oracle_A=r["mean_best"],
                    mean_regret_A=r["mean_regret"],
                    cpu_s=r["cpu_s"],
                    wall_s=r["wall_s_sum"],
                )
    path = HERE / "affinity/results/assays.csv"
    for r in read(path):
        if r["stage"] == "nominal":
            add(
                r["engine"],
                path,
                task="affinity",
                cohort="nominal",
                target=r["target"],
                assay=r["assay"],
                receptor=r["receptor"],
                planned=r["expected_compounds"],
                observed=r["observed_compounds"],
                rho=r["rho"],
                tau_b=r["tau_b"],
                correct_3x=r["correct_3x"],
                pairs_3x=r["pairs_3x"],
            )
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--append",
        action="store_true",
        help="add a new source snapshot, preserving older rounds",
    )
    args = p.parse_args()
    fresh = snapshot()
    if not fresh:
        raise ValueError("empty snapshot")
    old = read(LEDGER) if LEDGER.exists() else []
    existing = [r for r in old if r["round_source"] == fresh[0]["round_source"]]
    if existing:
        if existing != fresh:
            raise ValueError(
                "different measurements for an existing round; keep the original snapshot"
            )
        print(f"Verified {len(fresh)} records for {fresh[0]['round_source']}")
    elif args.append:
        with LEDGER.open("w", newline="") as f:
            w = csv.DictWriter(f, FIELDS, lineterminator="\n")
            w.writeheader()
            w.writerows(old + fresh)
    else:
        raise ValueError(
            "round missing from ledger; use --append after confirming inputs and evidence"
        )


if __name__ == "__main__":
    main()

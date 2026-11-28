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


def snapshot(build_file, pose_tables, assay_file):
    builds = json.loads(build_file.read_text())["versions"]
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

    for path, prefix in pose_tables:
        for r in read(path):
            if r["target"] == "ALL":
                add(
                    r["arm"],
                    path,
                    task="pose",
                    cohort=prefix + r["cohort"],
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
    if assay_file is None:
        return rows
    path = assay_file
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
    bench = HERE.parent
    panel = bench / "search_trial_panel/results"
    assays = bench / "search_trial_affinity"
    series = [
        (
            POSE / "builds.json",
            [
                (POSE / "results/summary.csv", ""),
                (HERE / "external/results/summary.csv", ""),
            ],
            HERE / "affinity/results/assays.csv",
        ),
        (
            panel / "builds.json",
            [
                (panel / "summary.csv", ""),
                (assays / "controls/four_versions/summary.csv", ""),
                (bench / "search_trial_nov28/results/summary.csv", "BRD4_transfer_"),
            ],
            assays / "results/assays.csv",
        ),
        (
            bench / "stage_local_structure/fixed/builds.json",
            [
                (bench / "stage_local_structure/fixed/summary.csv", ""),
                (bench / "stage_local_structure/external/summary.csv", "ESR1_transfer_"),
            ],
            bench / "stage_local_affinity/results/assays.csv",
        ),
        (
            bench / "fork_structure/builds.json",
            [
                (bench / "fork_structure/fixed/summary.csv", ""),
                (bench / "fork_structure/external/summary.csv", "MAPK14_transfer_"),
            ],
            bench / "fork_round/affinity/assays.csv",
        ),
    ]
    old = read(LEDGER) if LEDGER.exists() else []
    additions = []
    for build_file, pose_tables, assay_file in series:
        fresh = snapshot(build_file, pose_tables, assay_file)
        if build_file == bench / "fork_structure/builds.json":
            disclosed = bench / "fork_structure/disclosed"
            fresh.extend(snapshot(build_file, [
                (disclosed / "ESR1/summary.csv", "disclosed_ESR1_"),
                (disclosed / "BRD4/summary.csv", "disclosed_BRD4_"),
            ], None))
        revision = fresh[0]["round_source"]
        existing = [r for r in old if r["round_source"] == revision]
        if existing == fresh:
            print(f"Verified {len(fresh)} records for {revision}")
        elif args.append and existing and existing == fresh[:len(existing)]:
            additions.extend(fresh[len(existing):])
        elif existing:
            raise ValueError(f"Different measurements for existing round {revision}")
        elif not args.append:
            raise ValueError(
                f"Round {revision} missing; use --append after verifying evidence"
            )
        else:
            additions.extend(fresh)
    if additions:
        with LEDGER.open("w", newline="") as f:
            w = csv.DictWriter(f, FIELDS, lineterminator="\n")
            w.writeheader()
            w.writerows(old + additions)
        print(f"Added {len(additions)} records")


if __name__ == "__main__":
    main()

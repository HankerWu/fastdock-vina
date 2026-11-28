"""Summarize physical launches (one record per worker process) and observed CPU matches."""

import argparse
import json
from pathlib import Path


def audit(work):
    settings = json.loads((work / "common_cost.json").read_text())
    total = {}
    for engine in ("official", "parent", "stable", "candidate"):
        records = [json.loads(p.read_text()) for p in sorted((work / "outputs" / engine).glob("*/*.json"))]
        counter_only = [r for r in records if r["stage"] == "counter"]
        records = [r for r in records if r["stage"] != "counter"]
        binaries = {r["binary_sha256"] for r in records}
        assert len(binaries) == 1
        confirmations = [r for r in records if r["stage"] == settings["final"][engine]]
        assert len(confirmations) == 6 and all(r["valid"] for r in confirmations)
        cost = sum(r["user_s"] + r["sys_s"] for r in confirmations)
        reference = settings["target_s"]
        total[engine] = {
            "binary_sha256": binaries.pop(),
            "unique_launches": len(records),
            "valid_launches": sum(r["valid"] for r in records),
            "cpu_s": sum(r["user_s"] + r["sys_s"] for r in records),
            "sum_job_wall_s": sum(r["wall_s"] for r in records),
            "max_rss_kib": max(r["maxrss_kib"] for r in records),
            "separate_stock_counter_replays": len(counter_only),
            "extra_stock_counter_cpu_s": sum(r["user_s"] + r["sys_s"] for r in counter_only),
            "extra_stock_counter_sum_wall_s": sum(r["wall_s"] for r in counter_only),
            "cpu_control_final": settings["final"][engine],
            "cpu_control_s": cost,
            "cpu_control_planned_ligands": 6,
            "cpu_control_valid_ligands": len(confirmations),
            "cpu_target_s": reference,
            "within_15pct": abs(cost / reference - 1) <= 0.15,
        }
    return total


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.write_text(json.dumps(audit(args.work), indent=2) + "\n")

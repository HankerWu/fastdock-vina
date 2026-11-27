"""Run rigid lookahead and both counterfactuals on the same declared structure pairs."""

import argparse
import json
import sys
from functools import partial
from collections import defaultdict
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from search import translate_search

spec = spec_from_file_location(
    "rigid_scan_panel", HERE.parent / "rigid_translation_scan/run.py"
)
scan_panel = module_from_spec(spec)
spec.loader.exec_module(scan_panel)


def run_panel(base, meta, pairs, seeds, arms, work, probe, vina, stage, budget=9):
    queries = defaultdict(list)
    for receiver, query in pairs:
        queries[receiver].append(query)
    receivers = [(r, qs, base, meta, stage) for r, qs in queries.items()]
    scan_panel.evaluate_receivers(
        receivers,
        seeds,
        arms,
        work,
        probe,
        vina,
        partial(translate_search, budget=budget),
        lambda arm: 1 if arm == "raw" else 2,
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "probe", "vina"):
        p.add_argument("--" + name, type=Path, required=True)
    p.add_argument(
        "--arms",
        nargs="+",
        choices=["raw", "gradient", "sham"],
        default=["raw", "gradient", "sham"],
    )
    data = p.add_mutually_exclusive_group(required=True)
    data.add_argument("--dev", type=Path)
    data.add_argument("--secondary", type=Path)
    p.add_argument("--budget", type=int, choices=[9, 24, 48], default=9)
    a = p.parse_args()
    if a.dev:
        base, stage = a.dev, "development"
        meta = json.loads((base / "metadata.json").read_text())
        receivers = ("5C7N", "4OD0")
    else:
        base, stage = a.secondary, "secondary"
        receivers = ("1A28", "1SQN", "2P54", "1K7L")
        meta = {
            r: json.loads((base / r / "preparation.json").read_text())
            for r in receivers
        }
    pairs = [(r, q) for r in receivers for q in scan_panel.CASES[r]]
    run_panel(
        base,
        meta,
        pairs,
        [43103, 43117],
        a.arms,
        a.work,
        a.probe,
        a.vina,
        stage,
        a.budget,
    )


if __name__ == "__main__":
    main()

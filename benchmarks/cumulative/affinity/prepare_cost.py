"""Replay historical source preparation in isolation and time its actual stages."""

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import resource
import shutil
import sys
import time

SOURCE = Path(__file__).resolve().parents[2] / "geometry_comparison/affinity"


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("frozen", type=Path, help="extracted historical inputs")
    p.add_argument("scratch", type=Path, help="new directory for regeneration")
    p.add_argument("receipt", type=Path)
    args = p.parse_args()
    args.scratch.mkdir(parents=True)
    shutil.copytree(args.frozen / "raw", args.scratch / "raw")
    start = time.perf_counter(), time.process_time()
    spec = importlib.util.spec_from_file_location(
        "historical_prepare_affinity", SOURCE / "prepare_affinity.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    stages = []

    def wrap(name):
        original = getattr(module, name)

        def timed(*values):
            wall, cpu = time.perf_counter(), time.process_time()
            result = original(*values)
            label = values[1] if len(values) > 1 else None
            if isinstance(label, dict):
                label = label["cid"]
            stages.append(
                dict(
                    step=name,
                    id=label,
                    cpu_s=time.process_time() - cpu,
                    wall_s=time.perf_counter() - wall,
                )
            )
            return result

        setattr(module, name, timed)

    for step in ("cases", "ligand", "prepare_receptor"):
        wrap(step)
    sys.argv = [str(SOURCE / "prepare_affinity.py"), str(args.scratch)]
    module.main()
    wall, cpu = time.perf_counter() - start[0], time.process_time() - start[1]

    def sha(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()

    comparisons = []
    for path in sorted((args.scratch / "prepared").glob("*")):
        frozen = args.frozen / "prepared" / path.name
        comparisons.append(
            dict(file=path.name, generated=sha(path), frozen=sha(frozen) if frozen.is_file() else None)
        )
    for name in ("panel.csv", "exclusions.csv"):
        comparisons.append(
            dict(file=name, generated=sha(args.scratch / name), frozen=sha(SOURCE / name))
        )
    receipt = dict(
        wall_s=wall,
        cpu_s=cpu,
        maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        stages=stages,
        comparisons=comparisons,
    )
    args.receipt.write_text(json.dumps(receipt, indent=2) + "\n")


if __name__ == "__main__":
    main()

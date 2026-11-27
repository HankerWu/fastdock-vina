"""Five interleaved timing replicates for every eligible retinoid cell, with byte parity."""

import argparse
import itertools
import json
import random
import subprocess
import sys
from pathlib import Path

from panel import BENCH, SEEDS, holdout


def main(work, binaries):
    timer = BENCH / "rigid_flip/process_cost.py"
    meta = {
        q: json.loads((work / q / "preparation.json").read_text())
        for q in ("1FBY", "2LBD")
    }
    dest = work / "timing"
    dest.mkdir(exist_ok=True)
    rows = []
    order = random.Random(511041)
    cells = list(itertools.product(range(5), sorted(meta), SEEDS, (3000, 12000)))
    order.shuffle(cells)
    for repeat, q, s, cap in cells:
        configs = [
            (arm, c)
            for arm in binaries
            for c in ({cap} if arm == "parent" else (cap, int(cap * 0.85)))
        ]
        order.shuffle(configs)
        for arm, c in configs:
            key = f"{q}_{q}_{s}_{c}_{arm}"
            cmd = holdout.args_for(work, meta, q, q, s, c, dest / "out.pdbqt")
            subprocess.run(
                [
                    sys.executable,
                    "-S",
                    timer,
                    dest / "out.log",
                    dest / "cost.json",
                    binaries[arm].resolve(),
                    *cmd,
                    "--search_counters",
                ],
                check=True,
            )
            if (dest / "out.pdbqt").read_bytes() != (
                work / "runs" / (key + ".pdbqt")
            ).read_bytes():
                raise ValueError(key)
            rows.append(
                dict(
                    key=key,
                    query=q,
                    seed=s,
                    cap=c,
                    reference_cap=cap,
                    arm=arm,
                    repeat=repeat,
                    **json.loads((dest / "cost.json").read_text()),
                )
            )
    holdout.paired.write(work / "timings.csv", rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "parent", "always", "adaptive"):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    main(
        a.work.resolve(),
        {arm: getattr(a, arm) for arm in ("parent", "always", "adaptive")},
    )

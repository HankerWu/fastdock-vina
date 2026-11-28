"""Count stock-engine joint energy/gradient calls, validating unchanged output."""

import argparse
import json
import os
from pathlib import Path
from run import cells, digest, run_cell
from evaluate import write


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("engine", choices=("official", "parent"))
    p.add_argument("binary", type=Path)
    args = p.parse_args()
    binary = args.binary.resolve(strict=True)
    os.chdir(args.work)
    root = Path("outputs") / args.engine
    (root / "counter").mkdir(exist_ok=True, parents=True)
    rows = []
    for job in cells("low"):
        state, row, seed = job
        stem = f"{state}_{row['compound']}_{seed}"
        before = json.loads((root / "nominal" / (stem + ".json")).read_text())
        run_cell(job, args.engine, binary, "counter", 12000, digest(binary))
        after = json.loads((root / "counter" / (stem + ".json")).read_text())
        identical = (
            before["valid"] and after["valid"] and before["output_sha256"] == after["output_sha256"]
        )
        counts = [after[k] for k in ("hunt_evals", "second_evals", "explicit_evals")]
        assert identical and all(n is not None and n > 0 for n in counts)
        rows.append(
            dict(
                engine=args.engine,
                stem=stem,
                identical=identical,
                hunt_evals=counts[0],
                second_evals=counts[1],
                explicit_evals=counts[2],
                cpu_s=after["user_s"] + after["sys_s"],
                wall_s=after["wall_s"],
                stock_sha256=before["binary_sha256"],
                counter_sha256=after["binary_sha256"],
            )
        )
    write(root / "counter.csv", rows)


if __name__ == "__main__":
    main()

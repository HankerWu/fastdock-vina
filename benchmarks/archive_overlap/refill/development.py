"""Compare vacancy-triggered saving to the once-full rule on the complete development panel."""

import argparse
import csv
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "development"))
from dock_panel import BENCH, counts, holdout


def main(a):
    work = a.inputs.resolve()
    output = a.output.resolve()
    raw = output / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    meta = {
        q: json.loads((work / q / "preparation.json").read_text())
        for q in holdout.CASES
    }
    for q in meta:
        if not (output / q).exists():
            (output / q).symlink_to(work / q, target_is_directory=True)
    cases = [r for r in csv.DictReader(a.cells.open()) if r["arm"] == "fixed"]
    rows, modes = [], []
    for ordinal, case in enumerate(cases):
        q, r = case["query"], case["receiver"]
        for arm in ["fixed", "hold"] if ordinal % 2 == 0 else ["hold", "fixed"]:
            key = case["key"].rsplit("_", 1)[0] + "_" + arm
            out = raw / (key + ".pdbqt")
            log = raw / (key + ".log")
            cost = raw / (key + ".cost")
            cli = holdout.args_for(work, meta, q, r, case["seed"], case["cap"], out)
            cli[cli.index("--exhaustiveness") + 1] = "4"
            cli.append("--search_counters")
            subprocess.run(
                [
                    sys.executable,
                    "-S",
                    str(BENCH / "rigid_flip/process_cost.py"),
                    str(log),
                    str(cost),
                    str(getattr(a, arm).resolve()),
                    *cli,
                ],
                check=True,
            )
            digest = hashlib.sha256(out.read_bytes()).hexdigest()
            c = counts(log.read_text())
            skipped = []
            if arm == "fixed":
                if digest != case["sha256"]:
                    raise ValueError("Fixed output replay differs")
            else:
                other = out.with_suffix(".replay")
                cli[cli.index("--out") + 1] = str(other)
                cnt = subprocess.run(
                    [str(a.counted.resolve()), *cli],
                    capture_output=True,
                    text=True,
                    check=True,
                )
                (raw / (key + ".counts")).write_text(cnt.stdout)
                if out.read_bytes() != other.read_bytes() or counts(cnt.stdout) != c:
                    raise ValueError("Counter build differs")
                other.unlink()
                skipped = [
                    int(n) for n in re.findall(r"Skipped refills: (\d+)", cnt.stdout)
                ]
                if len(skipped) != 4:
                    raise ValueError("Expected four refill counts")
            row = dict(
                key=key,
                arm=arm,
                query=q,
                receiver=r,
                seed=case["seed"],
                cap=case["cap"],
                skipped_by_trajectory=",".join(map(str, skipped)),
                sha256=digest,
                **c,
                **json.loads(cost.read_text()),
            )
            modes.extend(holdout.paired.inspect(raw, meta[q], meta[r], row))
            rows.append(row)
    holdout.paired.write(output / "cells.csv", rows)
    holdout.paired.write(output / "modes.csv", modes)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("inputs", "cells", "output", "fixed", "hold", "counted"):
        p.add_argument("--" + k, type=Path, required=True)
    main(p.parse_args())

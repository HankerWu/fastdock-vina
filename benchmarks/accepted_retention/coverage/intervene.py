"""Counterfactual replacement of the last existing slot only in the frozen opportunities."""

import argparse
import csv
import os
import re
import subprocess
from pathlib import Path

import numpy as np
from probe import command, counters, mode_distances, references, rmsd, table


def main(a):
    os.sched_setaffinity(0, {a.cpu})
    a.output.mkdir(parents=True, exist_ok=True)
    results = []
    for row in csv.DictReader((a.results / "opportunities.csv").open()):
        if not row["oracle_step"]:
            continue
        key = row["key"]
        src = a.work / "runs" / key
        dest = a.output / key
        capture = dest.with_suffix(".conf")
        capture.unlink(missing_ok=True)
        env = dict(
            os.environ,
            FD_ORACLE_ACCEPT_STEP=row["oracle_step"],
            FD_ORACLE_CONFORMATION=str(capture.resolve()),
        )
        call = subprocess.run(
            command(a.work, row, a.binary, dest.with_suffix(".pdbqt")),
            env=env,
            capture_output=True,
            check=True,
        )
        native = counters(src.with_suffix(".log").read_bytes())
        modified = counters(call.stdout)
        if native[:3] != modified[:3] or len(native) != 4 or len(modified) != 4:
            raise ValueError("Intervention changed MC or merge")
        counts = [dict(re.findall(rb"(\w+)=(\d+)", x[-1])) for x in (native, modified)]
        if counts[0][b"explicit_inputs"] != counts[1][b"explicit_inputs"]:
            raise ValueError("Intervention changed slots")
        serial, refs = references(a.work, row)
        scores = mode_distances(dest.with_suffix(".pdbqt").read_bytes(), serial, refs)

        stages = {}
        for l in call.stderr.splitlines():
            if l.startswith((b"FDTRACE explicit_input ", b"FDTRACE explicit_output ")):
                _, name, idx, atoms = l.split(maxsplit=3)
                arr = np.fromstring(atoms.replace(b" ", b","), sep=",").reshape(-1, 4)
                lookup = dict(zip(arr[:, 0].astype(int), arr[:, 1:]))
                stages[name, int(idx)] = rmsd(
                    np.array([lookup[i] for i in serial]), refs
                )
        k = int(counts[1][b"explicit_inputs"])
        before, after = (stages[x, k] for x in (b"explicit_input", b"explicit_output"))
        if abs(before - float(row["skipped"])) > 1e-9:
            raise ValueError("Restored pose differs from nominated accepted minimum")
        results.append(
            {
                "key": key,
                "step": row["oracle_step"],
                "native_top1": row["top1"],
                "native_oracle": row["oracle"],
                "input_rmsd": before,
                "after_rmsd": after,
                "oracle_top1": scores[0],
                "oracle_top5": min(scores[:5]),
                "oracle_all": min(scores),
                "final_rank_of_hit": next(
                    (i for i, v in enumerate(scores, 1) if v <= 2), 0
                ),
                "slots": k,
                "native_explicit_evals": int(counts[0][b"explicit_evals"]),
                "counterfactual_explicit_evals": int(counts[1][b"explicit_evals"]),
            }
        )
        dest.with_suffix(".log").write_bytes(call.stdout)
    table(a.output / "interventions.csv", results)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "results", "output", "binary"):
        p.add_argument(name, type=Path)
    p.add_argument("--cpu", type=int, required=True)
    main(p.parse_args())

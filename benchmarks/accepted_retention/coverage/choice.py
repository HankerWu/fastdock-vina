"""For the three oracle opportunities, inspect the selector without changing a single output."""

import argparse
import csv
import re
import subprocess
from pathlib import Path

import numpy as np
from probe import STAGES, command, counters, table


def main(a):
    rows = []
    for r in csv.DictReader((a.results / "opportunities.csv").open()):
        if not r["oracle_step"]:
            continue
        path = a.work / "runs" / r["key"]
        saved = np.load(path.with_suffix(".npz"))
        names = np.asarray(STAGES)[saved["stage"]]
        known = {
            int(i): (d, x)
            for i, d, x in zip(
                saved["step"][names == "hunt_accepted"],
                saved["rmsd"][names == "hunt_accepted"],
                saved["xyz_milli"][names == "hunt_accepted"],
            )
        }
        merged = saved["xyz_milli"][names == "merged"].astype(float) / 1000
        out = a.work / (r["key"] + ".choice.pdbqt")
        obs = subprocess.run(
            command(a.work, r, a.observer, out), capture_output=True, check=True
        )
        blind = a.work / (r["key"] + ".policy.pdbqt")
        real = subprocess.run(
            command(a.work, r, a.candidate, blind), capture_output=True, check=True
        )
        if out.read_bytes() != blind.read_bytes() or counters(obs.stdout) != counters(
            real.stdout
        ):
            raise ValueError("Decision observer changed trajectory or output")
        fields = dict(
            re.findall(
                rb"(\w+)=([^ ]+)",
                next(l for l in obs.stdout.splitlines() if b"Deferred hunt:" in l),
            )
        )
        pool = []
        for l in obs.stderr.splitlines():
            if not l.startswith(b"FDPICK "):
                continue
            _, idx, energy, nearest, record = l.split()
            d, x = known[int(idx)]
            final = np.sqrt(
                np.square(merged - x.astype(float) / 1000).sum(axis=2).mean(axis=1)
            ).min()
            pool.append(
                {
                    "key": r["key"],
                    "step": int(idx),
                    "energy": float(energy),
                    "nearest_saved": float(nearest),
                    "nearest_merged": final,
                    "sets_record": int(record),
                    "rmsd": d,
                    "is_oracle": int(idx) == int(r["oracle_step"]),
                    "is_final_record": False,
                    "policy_replaced": int(fields[b"replaced"]),
                }
            )
        if len(pool) != int(r["n_skipped"]):
            raise ValueError("Incomplete skipped states")
        recs = [row for row in pool if row["sets_record"]]
        if recs:
            recs[-1]["is_final_record"] = True
        rows.extend(pool)
    table(a.results / "choices.csv", rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("work", "results", "observer", "candidate"):
        p.add_argument(k, type=Path)
    main(p.parse_args())

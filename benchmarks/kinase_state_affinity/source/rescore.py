"""Re-read the serialized top poses in the native C++ scorer, independently of the original run."""

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess
import tempfile
from statistics import median
from collections import defaultdict
import numpy as np


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def write(path, records):
    with path.open("w") as f:
        writer = csv.DictWriter(f, fieldnames=records[0])
        writer.writeheader()
        writer.writerows(records)


def verify(work, binary, jobs):
    inputs = []
    for cohort, root in [
        ("primary", work),
        ("replication", work / "independent-seeds"),
    ]:
        for arm, folder, table in [
            ("native", "pools", "results"),
            ("mixed", "mixed-pools", "mixed-results"),
            ("soft", "soft-pools", "soft-results"),
        ]:
            inputs.extend(
                (cohort, arm, root / folder, r)
                for r in read(root / table / "cells.csv")
            )

    def one(case):
        cohort, arm, folder, r = case
        rec, lig, seed = (r[k] for k in ["state", "ligand", "replicate"])
        path = folder / f"{rec}_{lig}_{seed}.pdbqt"
        top = path.read_text().split("ENDMDL")[0].split("\n", 1)[1]
        with tempfile.TemporaryDirectory() as temp:
            model = Path(temp) / "pose.pdbqt"
            model.write_text(top)
            result = subprocess.run(
                [
                    binary,
                    "--score_only",
                    "--cpu",
                    "1",
                    "--unbound_energy",
                    r["prepared_unbound"],
                    "--ligand",
                    model,
                    "--receptor",
                    work / f"prepared/{rec}.pdbqt",
                    "--config",
                    work / f"prepared/{rec}.box",
                ],
                capture_output=True,
                text=True,
                check=True,
            )
        energy = float(
            next(
                l
                for l in result.stdout.splitlines()
                if l.startswith("Estimated Free Energy of Binding")
            )
            .split(":")[1]
            .split()[0]
        )
        native = energy + float(r["normalization"]) * (
            float(r["prepared_unbound"]) - float(r["pool_unbound"])
        )
        return dict(
            cohort=cohort,
            arm=arm,
            target=r["target"],
            ligand=lig,
            state=rec,
            seed=seed,
            saved_native=float(r["native"]),
            readback_native=native,
            saved_fixed=float(r["energy"]),
            readback_fixed=energy,
        )

    with ThreadPoolExecutor(jobs) as pool:
        records = list(pool.map(one, inputs))
    write(work / "results/readback-cells.csv", records)
    kd = {(r["target"], r["cid"]): float(r["kd_nm"]) for r in read(work / "panel.csv")}
    output = []
    for score in ["native", "fixed"]:
        for target in sorted({t for t, lig in kd}):
            ligands = [lig for t, lig in kd if t == target]
            for arm in ["native", "mixed", "soft"]:
                rows = [r for r in records if r["arm"] == arm and r["target"] == target]
                states = sorted({r["state"] for r in rows})
                vals = defaultdict(list)
                for row in rows:
                    for tag in ["saved", "readback"]:
                        vals[tag, row["ligand"], row["state"]].append(
                            row[tag + "_" + score]
                        )
                same = changed = 0
                for a in ligands:
                    for b in ligands:
                        if kd[target, b] < 3 * kd[target, a]:
                            continue
                        signs = []
                        for tag in ["saved", "readback"]:
                            energies = np.array(
                                [
                                    [median(vals[tag, c, s]) for s in states]
                                    for c in (a, b)
                                ]
                            )
                            ka, kb = np.exp(
                                -energies / (0.00198720425864083 * 298.15)
                            ).sum(1)
                            signs.append(np.sign(ka - kb))
                        if signs[0] == signs[1]:
                            same += 1
                        else:
                            changed += 1
                output.append(
                    dict(
                        target=target,
                        arm=arm,
                        score=score,
                        unchanged=same,
                        flipped=changed,
                        max_abs_pose_energy_delta=max(
                            abs(r["saved_" + score] - r["readback_" + score])
                            for r in rows
                        ),
                    )
                )
    write(work / "results/readback-summary.csv", output)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("binary", type=Path)
    p.add_argument("--jobs", default=3, type=int)
    args = p.parse_args()
    verify(args.work.resolve(), args.binary.resolve(), args.jobs)

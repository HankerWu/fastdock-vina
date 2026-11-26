"""Rescore serialized poses with a prepared topology, then compare unique Kd pairs."""

import argparse
import csv
import itertools
import math
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from statistics import median
import json

from run import ARMS, SEEDS, STATES


def write(path, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0])
        w.writeheader()
        w.writerows(rows)


def num(s, tag):
    return float(
        next(l for l in s.splitlines() if l.startswith(tag)).split(":", 1)[1].split()[0]
    )


def compare(delta):
    return 1 if delta > 1e-8 else 0 if delta < -1e-8 else 0.5


def agg(e, kind):
    if kind == "mean":
        return sum(e) / len(e)
    if kind == "exponential":
        t = 0.00198720425864083 * 298.15
        m = min(e)
        return m - t * math.log(sum(math.exp((m - x) / t) for x in e) / len(e))
    return e[int(kind)]


def main(work, binary, jobs):
    rows = list(csv.DictReader((work / "panel.csv").open()))
    prepared = work / "prepared"
    refs = {}
    for cid in sorted({r["cid"] for r in rows}):
        vals = []
        for u in (0, 100):
            cmd = [
                binary,
                "--ligand",
                prepared / (cid + ".pdbqt"),
                "--receptor",
                prepared / "2E9P.pdbqt",
                "--score_only",
                "--autobox",
                "--unbound_energy",
                str(u),
                "--cpu",
                "1",
            ]
            log = subprocess.check_output(cmd, text=True)
            vals.append(num(log, "Estimated Free Energy of Binding"))
        refs[cid] = (
            num(log, "(2) Final Total Internal Energy"),
            (vals[0] - vals[1]) / 100,
        )

    def cell(z):
        r, state, k, arm = z
        cid = r["cid"]
        path = work / "pools" / f"{state}_{cid}_{k}_{arm}.pdbqt"
        pool = path.read_text()
        if not pool.endswith("ENDMDL\n"):
            raise ValueError(path)
        top = pool.split("ENDMDL")[0]
        u = num(top, "REMARK UNBOUND:")
        score = num(top, "REMARK VINA RESULT:")
        total = num(top, "REMARK INTER + INTRA:")
        pre, slope = refs[cid]
        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / "top.pdbqt"
            f.write_text(top.split("\n", 1)[1])
            cmd = [
                binary,
                "--ligand",
                prepared / (cid + ".pdbqt"),
                "--receptor",
                prepared / (state + ".pdbqt"),
                "--pose_coordinates",
                f,
                "--score_only",
                "--cpu",
                "1",
                "--unbound_energy",
                str(u),
            ]
            center = json.loads((prepared / (state + ".json")).read_text())["center"]
            for ax, coord in zip("xyz", center):
                cmd.extend(["--center_" + ax, str(coord), "--size_" + ax, "26"])
            proc = subprocess.run(cmd, text=True, capture_output=True)
            if proc.returncode != 0:
                raise RuntimeError((path, proc.stderr))
            rb = num(proc.stdout, "Estimated Free Energy of Binding")
        return dict(
            target=r["target"],
            cid=cid,
            state=state,
            seed=k,
            arm=arm,
            pool_score=score,
            readback=rb,
            prepared_ref_score=slope * (total - pre),
            pool_unbound=u,
            prepared_unbound=pre,
            normalization=slope,
            modes=pool.count("ENDMDL"),
        )

    cases = [
        (r, s, k, a)
        for r in rows
        for s in STATES[r["target"]]
        for k in SEEDS
        for a in ARMS
    ]
    with ThreadPoolExecutor(jobs) as e:
        cells = list(e.map(cell, cases))
    write(work / "cells.csv", cells)
    idx = {(r["cid"], r["state"], r["seed"], r["arm"]): r for r in cells}
    pairs = []
    summary = []
    for target, states in STATES.items():
        molecules = [r for r in rows if r["target"] == target]
        chem_pairs = [
            (a, b)
            for a, b in itertools.permutations(molecules, 2)
            if float(a["kd_nm"]) * 3 <= float(b["kd_nm"])
        ]
        for cohort, seeds in [
            ("all", SEEDS),
            ("first", SEEDS[:2]),
            ("last", SEEDS[2:]),
        ]:
            for field in ["pool_score", "readback", "prepared_ref_score"]:
                for pooling in ["mean", "exponential", "0", "1"]:
                    outcomes = {arm: [] for arm in ARMS}
                    for a, b in chem_pairs:
                        d = dict(
                            target=target,
                            strong=a["compound"],
                            weak=b["compound"],
                            ratio=float(b["kd_nm"]) / float(a["kd_nm"]),
                            cohort=cohort,
                            field=field,
                            pooling=pooling,
                        )
                        for arm in ARMS:
                            energies = [
                                agg(
                                    [
                                        median(
                                            idx[x["cid"], s, k, arm][field]
                                            for k in seeds
                                        )
                                        for s in states
                                    ],
                                    pooling,
                                )
                                for x in (a, b)
                            ]
                            delta = energies[1] - energies[0]
                            d[arm] = delta
                            outcomes[arm].append(compare(delta))
                        pairs.append(d)
                    result = dict(
                        target=target,
                        cohort=cohort,
                        field=field,
                        pooling=pooling,
                        n=len(chem_pairs),
                        size_only=sum(
                            compare(int(a["heavy"]) - int(b["heavy"]))
                            for a, b in chem_pairs
                        ),
                    )
                    for arm, values in outcomes.items():
                        result[arm] = sum(values)
                        result[arm + "_gain"] = sum(
                            v > n for v, n in zip(values, outcomes["native"])
                        )
                        result[arm + "_loss"] = sum(
                            v < n for v, n in zip(values, outcomes["native"])
                        )
                    summary.append(result)
    write(work / "pairs.csv", pairs)
    write(work / "summary.csv", summary)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("binary", type=Path)
    p.add_argument("--jobs", type=int, default=2)
    a = p.parse_args()
    main(a.work.resolve(), a.binary.resolve(), a.jobs)

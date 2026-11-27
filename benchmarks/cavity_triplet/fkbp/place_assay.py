"""Frozen mapping of blind one-shot placement onto the paired FKBP12 AD4 runs."""

import argparse
import csv
import itertools
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tarfile
import time

import numpy as np
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OLD = ROOT / "benchmarks/receptor_linkage/aurora/internal_loops/rigid_fkbp"
sys.path.insert(0, str(OLD))
from run_panel import BENCH, RECEIVERS, SEED, SOURCE, cpu_seconds, prepare, sha, table


def inputs(work):
    hashes = prepare(work)
    originals = BENCH / "human_panel/prepared.tar.xz"
    archive = OLD / "pools.tar.xz"
    with tarfile.open(archive) as t:
        t.extractall(work, filter="data")
    meta = json.loads((SOURCE / "results/metadata.json").read_text())
    with tarfile.open(originals) as t:
        for rec in RECEIVERS:
            member = rec + "/receptor.pdbqt"
            t.extract(member, work, filter="data")
            assert sha(work / member) == meta["receptors"][rec]
    for path in (originals, archive, SOURCE / "results/metadata.json"):
        hashes[str(path.relative_to(ROOT))] = sha(path)
    (work / "inputs_fkbp.json").write_text(json.dumps(hashes, indent=2) + "\n")
    return meta


def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def invoke(work, exe, rec, compound, arm, env):
    parent = Path("pools") / arm
    (work / parent).mkdir(exist_ok=True)
    stem = parent / f"{rec}_{compound}"
    source = (work / "pools/base" / (stem.name + ".conf")).read_text()
    conf = stem.with_suffix(".conf")
    (work / conf).write_text(source.replace("pools/base/", str(parent) + "/"))
    settings = {k: v for k, v in os.environ.items() if not k.startswith("FD_")}
    start, cpu = time.perf_counter(), cpu_seconds()
    with (work / stem.with_suffix(".log")).open("w") as f:
        subprocess.run(
            [str(exe), "--config", str(conf), "--search_counters"],
            cwd=work,
            env=settings | {k: str(v) for k, v in env.items()},
            stdout=f,
            stderr=subprocess.STDOUT,
            check=True,
        )
    wall_s, cpu_s = time.perf_counter() - start, cpu_seconds() - cpu
    pool = work / stem.with_suffix(".pdbqt")
    log = (work / stem.with_suffix(".log")).read_text()
    counters = [
        dict(re.findall(r"(\w+)=(\d+)", line))
        for line in log.splitlines()
        if line.startswith("Search trajectory:")
    ]
    assert len(counters) == 3
    return dict(
        arm=arm,
        receiver=rec,
        compound=compound,
        seed=SEED,
        wall_s=wall_s,
        cpu_s=cpu_s,
        hunt_evals=sum(int(s["hunt_evals"]) for s in counters),
        second_evals=sum(int(s["second_evals"]) for s in counters),
        steps=sum(int(s["steps"]) for s in counters),
        modes=pool.read_text().count("ENDMDL"),
        sha256=sha(pool),
        executable_sha256=sha(exe),
    )


def sources(work, exe):
    recorded = []
    capture = work / "capture"
    capture.mkdir(exist_ok=True)
    for rec, n in itertools.product(RECEIVERS, range(4, 12)):
        p = capture / f"{rec}_{n}"
        row = invoke(
            work,
            exe,
            rec,
            n,
            "capture",
            dict(
                FD_SOURCE_STEP=1,
                FD_SOURCE_TRAJECTORY=1,
                FD_SOURCE_XYZ=p.with_suffix(".xyz"),
                FD_SOURCE_Q=p.with_suffix(".q"),
            ),
        )
        row["historical_sha256"] = sha(work / "pools/base" / (p.name + ".pdbqt"))
        assert row["sha256"] == row["historical_sha256"]
        recorded.append(row)
    table(work / "capture.csv", recorded)


def placements(work, meta):
    sys.path.insert(0, str(ROOT / "tools"))
    from pocket_triangle import Pocket

    pocket_rows, timings = [], []
    plans = work / "plans"
    plans.mkdir(exist_ok=True)
    for rec in RECEIVERS:
        start, cpu = time.perf_counter(), time.process_time()
        pocket = Pocket.from_pdbqt(
            work / rec / "receptor.pdbqt", meta["centers_A"][rec], 24
        )
        pocket_rows.append(
            dict(
                receiver=rec,
                wall_s=time.perf_counter() - start,
                cpu_s=time.process_time() - cpu,
                grid_points=pocket.n_grid,
                burial_points=pocket.n_burial,
                sites=len(pocket.sites),
                triples=len(pocket.triangles),
            )
        )
        np.save(plans / (rec + "_sites.npy"), pocket.sites)
        for n in range(4, 12):
            stem = rec + "_" + str(n)
            lines = (work / "capture" / (stem + ".xyz")).read_text().splitlines()
            serials = np.array(lines[0].lstrip("# ").split(), int)
            source = np.array(lines[1].split()[3:], float).reshape(-1, 3)
            q = np.loadtxt(work / "capture" / (stem + ".q"))[2:]
            start, cpu = time.perf_counter(), time.process_time()
            timing = {}
            res = pocket.place(source, SEED, timing)
            wall_s, cpu_s = time.perf_counter() - start, time.process_time() - cpu
            timings.append(
                dict(
                    receiver=rec,
                    compound=n,
                    wall_s=wall_s,
                    cpu_s=cpu_s,
                    field_queries_per_arm=len(source) * len(res["center"]),
                    rigid_fits=len(res["center"]),
                    **timing,
                )
            )
            arrays = dict(serials=serials, xyz=source, q=q, **res)
            for arm, key in (("match", "matched"), ("shuffle", "sham")):
                rotm = res[key + "_rotation"]
                winner = res[key + "_best"]
                rot = Rotation.from_matrix(rotm[winner])
                origin = rot.apply(q[:3] - source.mean(axis=0)) + res["center"][winner]
                quat = (rot * Rotation.from_quat(q[[4, 5, 6, 3]])).as_quat()
                target = np.r_[1, 1, origin, quat[[3, 0, 1, 2]]]
                (plans / (stem + "_" + arm + ".input")).write_text(
                    " ".join(["1", "1", *[f"{x:.17g}" for x in target[2:]]]) + "\n"
                )
            np.savez_compressed(plans / (stem + ".npz"), **arrays)
    table(work / "pocket.csv", pocket_rows)
    table(work / "placement.csv", timings)


def pairs(work, exe):
    collected = []
    for i, (rec, n) in enumerate(itertools.product(RECEIVERS, range(4, 12))):
        stem = f"{rec}_{n}"
        for arm in ("match", "shuffle") if i % 2 == 0 else ("shuffle", "match"):
            p = work / "pools" / arm / stem
            row = invoke(
                work,
                exe,
                rec,
                n,
                arm,
                dict(
                    FD_RIGID_INTERVENTION=work
                    / "plans"
                    / (stem + "_" + arm + ".input"),
                    FD_RIGID_LOG=p.with_suffix(".stages"),
                    FD_RIGID_COORDS=p.with_suffix(".xyz"),
                    FD_RIGID_HISTORY=p.with_suffix(".history"),
                ),
            )
            recorded = np.array(
                p.with_suffix(".stages").read_text().splitlines()[0].split()[1:], float
            )
            np.testing.assert_array_equal(
                recorded, np.loadtxt(work / "capture" / (stem + ".q"))[2:]
            )
            collected.append(row)
    for arm in ("match", "shuffle"):
        table(work / f"runs-{arm}.csv", [r for r in collected if r["arm"] == arm])


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=("prepare", "capture", "place", "run"))
    p.add_argument("work", type=lambda s: Path(s).resolve())
    p.add_argument("--exe", type=lambda s: Path(s).resolve())
    a = p.parse_args()
    if a.action == "prepare":
        inputs(a.work)
    if a.action == "capture":
        sources(a.work, a.exe)
    if a.action == "place":
        placements(a.work, json.loads((SOURCE / "results/metadata.json").read_text()))
    if a.action == "run":
        pairs(a.work, a.exe)

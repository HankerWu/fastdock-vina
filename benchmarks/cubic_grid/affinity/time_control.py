"""Calibrate native search budgets using a fresh seed, then time all paired cells serially."""

import csv
import hashlib
import itertools
import os
import resource
import statistics
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from time import monotonic
from run import SEEDS, STATES, box_args

work, binary, out = map(lambda p: Path(p).resolve(), sys.argv[1:4])
out.mkdir(exist_ok=True, parents=True)
core = max(os.sched_getaffinity(0))
os.sched_setaffinity(0, {core})
panel = list(csv.DictReader((work / "panel.csv").open()))
prepared = work / "prepared"
results = []
fields = [
    "target",
    "cid",
    "seed",
    "arm",
    "budget",
    "wall",
    "cpu",
    "score",
    "sha256",
    "core",
]


def write(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


def trial(r, seed, arm, cap):
    state = STATES[r["target"]]
    stem = f"{state}_{r['cid']}_{seed}_{arm}_{cap}"
    dest = out / (stem + ".pdbqt")
    log = out / (stem + ".log")
    cmd = [
        binary,
        "--ligand",
        prepared / (r["cid"] + ".pdbqt"),
        "--receptor",
        prepared / (state + ".pdbqt"),
        "--pose_selection",
        "late",
        "--cpu",
        "1",
        "--seed",
        str(seed),
        "--exhaustiveness",
        "2",
        "--num_modes",
        "9",
        "--energy_range",
        "100",
        "--max_evals",
        str(cap),
        "--search_counters",
        "--out",
        dest,
        *box_args(prepared, state),
    ]
    if arm == "cubic":
        cmd.append("--cubic_search_grid")
    start = monotonic()
    t0 = resource.getrusage(resource.RUSAGE_CHILDREN)
    with log.open("w") as f:
        subprocess.run(cmd, check=True, stdout=f, stderr=subprocess.STDOUT)
    t1 = resource.getrusage(resource.RUSAGE_CHILDREN)
    elapsed = monotonic() - start
    cpu = t1.ru_stime + t1.ru_utime - t0.ru_stime - t0.ru_utime
    data = dest.read_bytes()
    if seed in SEEDS and cap == 12000:
        assert (
            data
            == (work / "pools" / f"{state}_{r['cid']}_{seed}_{arm}.pdbqt").read_bytes()
        )
    score = float(
        next(
            line.split(":", 1)[1].split()[0]
            for line in data.decode().splitlines()
            if line.startswith("REMARK VINA RESULT:")
        )
    )
    return dict(
        target=r["target"],
        cid=r["cid"],
        seed=seed,
        arm=arm,
        budget=cap,
        wall=round(elapsed, 6),
        cpu=round(cpu, 6),
        score=score,
        sha256=hashlib.sha256(data).hexdigest(),
        core=core,
    )


if "--extend" in sys.argv[4:]:
    with (out / "runs.csv").open() as f:
        results = list(csv.DictReader(f))
    for row in results:
        for key in ("wall", "cpu", "score"):
            row[key] = float(row[key])
else:
    cal = []
    budget = {}
    for target in STATES:
        rows = sorted(
            [r for r in panel if r["target"] == target],
            key=lambda r: (int(r["heavy"]), r["cid"]),
        )
        for n, i in enumerate((0, len(rows) // 2, len(rows) - 1)):
            methods = [("native", 8000), ("native", 20000), ("cubic", 12000)]
            for method, cap in methods[n:] + methods[:n]:
                cal.append(trial(rows[i], 78317, method, cap))

        def duration(arm, cap):
            return sum(
                r["wall"]
                for r in cal
                if r["target"] == target and r["arm"] == arm and r["budget"] == cap
            )

        t8, t20, tc = (
            duration("native", 8000),
            duration("native", 20000),
            duration("cubic", 12000),
        )
        assert t20 > t8
        budget[target] = max(
            8000, round((8000 + 12000 * (tc - t8) / (t20 - t8)) / 500) * 500
        )
    write(out / "calibration.csv", cal)
    for i, (seed, r) in enumerate(itertools.product(SEEDS, panel)):
        methods = [
            ("native", 12000),
            ("cubic", 12000),
            ("native_time", budget[r["target"]]),
        ]
        k = i % 3
        for arm, cap in methods[k:] + methods[:k]:
            results.append(trial(r, seed, arm, cap))
for seed, r in itertools.product(SEEDS, panel):
    results.append(trial(r, seed, "native_double", 24000))
write(out / "runs.csv", results)
summary = []
for t in STATES:
    rows = [r for r in panel if r["target"] == t]
    kd = {r["cid"]: float(r["kd_nm"]) for r in rows}
    pair = [(u, v) for u in kd for v in kd if kd[v] >= 3 * kd[u]]
    for arm in ("native", "cubic", "native_time", "native_double"):
        batch = [r for r in results if r["target"] == t and r["arm"] == arm]
        scores = {
            cid: statistics.mean(r["score"] for r in batch if r["cid"] == cid)
            for cid in kd
        }
        total = sum(
            1
            if scores[a] + 1e-8 < scores[b]
            else 0
            if scores[b] + 1e-8 < scores[a]
            else 0.5
            for a, b in pair
        )
        summary.append(
            dict(
                target=t,
                arm=arm,
                cap=batch[0]["budget"],
                runs=len(batch),
                wall_s=round(sum(r["wall"] for r in batch), 3),
                cpu_s=round(sum(r["cpu"] for r in batch), 3),
                correct=total,
                pairs=len(pair),
            )
        )
write(out / "time_summary.csv", summary)
print((out / "time_summary.csv").read_text())

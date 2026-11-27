"""Run the frozen cross-dock searches on the prepared four-target transfer panel."""

import argparse
import hashlib
import itertools
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "mobile_water"))
from replay import save

SEEDS = (43103, 43117, 43133)
BUDGETS = (3000, 12000)
TARGETS = {"BACE1", "EPHX2", "BRD9", "CHEK1"}


def main(binary, work, jobs):
    meta = json.loads((work / "metadata.json").read_text())
    if len(meta) != 8 or {v["target"] for v in meta.values()} != TARGETS:
        raise ValueError(
            "Expected both original receptors for each of the four targets"
        )
    for cap in BUDGETS:
        (work / f"runs_{cap}").mkdir(exist_ok=True)
    tasks = [
        (q, r, cap, seed, arm)
        for q, r in itertools.permutations(sorted(meta), 2)
        if meta[q]["target"] == meta[r]["target"]
        for cap, seed, arm in itertools.product(BUDGETS, SEEDS, ("native", "soft"))
    ]

    def dock(task):
        q, r, cap, seed, arm = task
        prefix = Path(f"runs_{cap}/{q}_{r}_{seed}_{arm}")
        args = [
            str(binary),
            "--receptor",
            f"{r}/receptor.pdbqt",
            "--ligand",
            f"{q}/ligand.pdbqt",
            "--pose_selection",
            "late",
            "--exhaustiveness",
            "4",
            "--refine_candidates",
            "40",
            "--num_modes",
            "10",
            "--energy_range",
            "100",
            "--max_evals",
            str(cap),
            "--cpu",
            "1",
            "--seed",
            str(seed),
            "--search_counters",
            "--out",
            str(prefix.with_suffix(".pdbqt")),
        ]
        if arm == "soft":
            args.append("--soft_receptor_search")
        for ax, value in zip("xyz", meta[r]["center"]):
            args.extend([f"--center_{ax}", str(value), f"--size_{ax}", "24"])
        start = time.monotonic()
        with (work / prefix.with_suffix(".log")).open("w") as f:
            proc = subprocess.Popen(args, cwd=work, stdout=f, stderr=subprocess.STDOUT)
            _, status, usage = os.wait4(proc.pid, 0)
            proc.returncode = os.waitstatus_to_exitcode(status)
        return dict(
            target=meta[q]["target"],
            query=q,
            receiver=r,
            cap=cap,
            seed=seed,
            search=arm,
            exit=proc.returncode,
            wall_s=time.monotonic() - start,
            cpu_s=usage.ru_utime + usage.ru_stime,
            max_rss_kib=usage.ru_maxrss,
        )

    with ThreadPoolExecutor(jobs) as pool:
        rows = list(pool.map(dock, tasks))
    save(work / "search_cost.csv", rows)
    (work / "binary.json").write_text(
        json.dumps(
            dict(
                sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                workers=jobs,
                version=subprocess.check_output(
                    [str(binary), "--version"], text=True
                ).strip(),
            ),
            indent=2,
        )
        + "\n"
    )
    if any(r["exit"] for r in rows):
        raise RuntimeError("Failed transfer runs are listed in search_cost.csv")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("binary", type=Path)
    p.add_argument("work", type=Path)
    p.add_argument("--jobs", type=int, default=2)
    a = p.parse_args()
    main(a.binary.resolve(), a.work.resolve(), a.jobs)

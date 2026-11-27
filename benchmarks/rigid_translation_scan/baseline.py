"""Fresh current-Pro and upstream cross-docks on external PGR and PPARalpha."""

import argparse
import csv
import itertools
import json
import os
import subprocess
import time
from pathlib import Path

import numpy as np

from run import CASES, paired, rms, save


def main():
    p = argparse.ArgumentParser()
    for f in ("hold", "dev", "work", "pro", "upstream"):
        p.add_argument("--" + f, type=Path, required=True)
    a = p.parse_args()
    a.work.mkdir(exist_ok=True, parents=True)
    cells = []
    costs = []
    for f in ("geometry", "cost"):
        with (a.dev / (f + ".csv")).open() as stream:
            rows = list(csv.DictReader(stream))
        if f == "cost":
            for r in rows:
                if r["receiver"] != r["query"]:
                    costs.append(
                        {
                            k: r[k]
                            for k in (
                                "receiver",
                                "query",
                                "cap",
                                "seed",
                                "engine",
                                "wall_s",
                                "cpu_s",
                                "peak_rss_kib",
                                "search_evals",
                                "explicit_evals",
                            )
                        }
                    )
        else:
            for key, group in itertools.groupby(
                sorted(
                    [r for r in rows if r["receiver"] != r["query"]],
                    key=lambda r: (
                        r["receiver"],
                        r["query"],
                        r["cap"],
                        r["seed"],
                        r["engine"],
                        int(r["model"]),
                    ),
                ),
                key=lambda r: (
                    r["receiver"],
                    r["query"],
                    r["cap"],
                    r["seed"],
                    r["engine"],
                ),
            ):
                scores = [float(r["rmsd"]) for r in group]
                cells.append(
                    dict(zip(("receiver", "query", "cap", "seed", "engine"), key))
                    | {
                        "stage": "development",
                        "top": scores[0],
                        "top5": min(scores[:5]),
                        "oracle": min(scores),
                    }
                )
    for rec in ["1A28", "1SQN", "2P54", "1K7L"]:
        (query,) = CASES[rec]
        meta = {
            e: json.loads((a.hold / e / "preparation.json").read_text())
            for e in [rec, query]
        }
        rot, shift, *_ = paired.alignment(meta[query], meta[rec])
        refs = np.array(meta[query]["mapped"]) @ rot + shift
        for seed, cap, engine in itertools.product(
            [43103, 43117], [3000, 12000], ("pro", "upstream")
        ):
            path = a.work / f"{query}_{rec}_{cap}_{seed}_{engine}.pdbqt"
            cmd = [
                str(getattr(a, engine)),
                "--receptor",
                str(a.hold / rec / "receptor.pdbqt"),
                "--ligand",
                str(a.hold / query / "ligand.pdbqt"),
                "--cpu",
                "1",
                "--seed",
                str(seed),
                "--exhaustiveness",
                "4",
                "--num_modes",
                "9",
                "--max_evals",
                str(cap),
                "--energy_range",
                "100",
                "--out",
                str(path),
            ]
            for axis, val in zip("xyz", meta[rec]["center"]):
                cmd += ["--center_" + axis, str(val), "--size_" + axis, "24"]
            if engine == "pro":
                cmd += ["--search_counters"]
            log = path.with_suffix(".log")
            start = time.perf_counter()
            with log.open("w") as f:
                proc = subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT)
                _, s, u = os.wait4(proc.pid, 0)
                proc.returncode = os.waitstatus_to_exitcode(s)
            if proc.returncode:
                raise RuntimeError(str(path))
            import re

            output = log.read_text()
            search = re.findall(r"hunt_evals=(\d+).*?second_evals=(\d+)", output)
            finals = re.findall(r"explicit_evals=(\d+)", output)
            costs.append(
                {
                    "receiver": rec,
                    "query": query,
                    "cap": cap,
                    "seed": seed,
                    "engine": engine,
                    "wall_s": time.perf_counter() - start,
                    "cpu_s": u.ru_utime + u.ru_stime,
                    "peak_rss_kib": u.ru_maxrss,
                    "search_evals": sum(int(x) + int(y) for x, y in search)
                    if search
                    else "",
                    "explicit_evals": sum(map(int, finals)) if finals else "",
                }
            )
            rmsds = []
            for block in path.read_text().split("ENDMDL"):
                if "MODEL" not in block:
                    continue
                atoms = {
                    int(l[6:11]): [float(l[j : j + 8]) for j in (30, 38, 46)]
                    for l in block.splitlines()
                    if l.startswith(("ATOM", "HETATM"))
                }
                rmsds.append(
                    rms(np.array([atoms[i] for i in meta[query]["serial"]]), refs)
                )
            cells.append(
                {
                    "receiver": rec,
                    "query": query,
                    "cap": cap,
                    "seed": seed,
                    "engine": engine,
                    "stage": "transfer",
                    "top": rmsds[0],
                    "top5": min(rmsds[:5]),
                    "oracle": min(rmsds),
                }
            )
        save(a.work / "baseline_cells.csv", cells)
        save(a.work / "baseline_cost.csv", costs)
        print(rec, flush=True)


if __name__ == "__main__":
    main()

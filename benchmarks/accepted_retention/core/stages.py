"""Trace only actual swaps, including the displaced slot, with full-precision coordinates."""

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

from replay import PRIOR, paired

sys.path.insert(0, str(PRIOR / "holdout"))
from run_holdout import stage_data


def run(args):
    work = args.work.resolve()
    runs = work / "runs"
    out = work / "stage_runs"
    out.mkdir(exist_ok=True)
    meta = {
        p.name: json.loads((p / "preparation.json").read_text())
        for p in work.iterdir()
        if (p / "preparation.json").exists()
    }
    paired.renin_correspondence(work, meta)
    all_records = []
    with (work / "cells.csv").open() as f:
        for cell in csv.DictReader(f):
            if (cell["arm"], cell["replaced"]) != ("candidate", "1"):
                continue
            s, r = cell["query"], cell["receiver"]
            slot = int(cell["slot"])
            paired_stages = {}
            for arm in ("parent", "candidate"):
                key = "_".join((s, r, cell["seed"], cell["cap"], arm))
                cmd = [
                    getattr(args, arm),
                    "--cpu",
                    "1",
                    "--exhaustiveness",
                    "1",
                    "--seed",
                    cell["seed"],
                    "--max_evals",
                    cell["cap"],
                    "--num_modes",
                    "9",
                    "--energy_range",
                    "100",
                    "--ligand",
                    work / s / "ligand.pdbqt",
                    "--receptor",
                    work / r / "receptor.pdbqt",
                    "--out",
                    out / (key + ".pdbqt"),
                    "--search_counters",
                ]
                for axis, v in zip("xyz", meta[r]["center"]):
                    cmd += ["--center_" + axis, str(v), "--size_" + axis, "24"]
                proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
                log = proc.stdout + proc.stderr
                (out / (key + ".log")).write_text(log)
                if (out / (key + ".pdbqt")).read_bytes() != (
                    runs / (key + ".pdbqt")
                ).read_bytes():
                    raise ValueError("Instrumented output differs " + key)
                frozen = (runs / (key + ".log")).read_text()
                if paired.counts(log) != paired.counts(frozen):
                    raise ValueError("Instrumented count mismatch " + key)
                records = stage_data(log, key, meta[s], meta[r])
                paired_stages[arm] = {
                    (x["stage"], int(x["input_index"])): x for x in records
                }
                for x in records:
                    x["arm"] = arm
                    x["swapped_slot"] = slot
                all_records += records
            pn, cn = paired_stages["parent"], paired_stages["candidate"]
            if pn.keys() != cn.keys():
                raise ValueError("Mismatched explicit input count")
            for stage, i in pn:
                if (
                    i != slot
                    and pn[stage, i]["coords_sha256"] != cn[stage, i]["coords_sha256"]
                ):
                    raise ValueError("Other slot changed")
    paired.write(work / "stages.csv", all_records)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    for arm in ("parent", "candidate"):
        p.add_argument("--" + arm, type=Path, required=True)
    run(p.parse_args())

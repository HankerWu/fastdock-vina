"""Replay the first redirected proposal at every full-budget development condition."""

import argparse
import csv
import json
from pathlib import Path
import re
import subprocess
import sys

import numpy as np

BENCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCH / "rigid_flip"))
from holdout import args_for, paired


def run(args):
    work = args.work.resolve()
    out = args.out.resolve()
    out.mkdir(exist_ok=True, parents=True)
    with (work / "cells.csv").open() as f:
        cells = list(csv.DictReader(f))
    stages = []
    for cell in cells:
        arm = cell["arm"]
        if arm not in {"always", "adaptive"} or int(cell["cap"]) != 12000:
            continue
        q, r = cell["query"], cell["receiver"]
        meta = {
            s: json.loads((work / s / "preparation.json").read_text()) for s in (q, r)
        }
        path = out / (cell["key"] + ".pdbqt")
        cli = args_for(work, meta, q, r, cell["seed"], cell["cap"], path)
        result = subprocess.run(
            [getattr(args, arm), *cli, "--search_counters"],
            check=True,
            capture_output=True,
            text=True,
        )
        if path.read_bytes() != (work / "runs" / path.name).read_bytes():
            raise ValueError("Observer altered the final output: " + path.name)
        path.unlink()
        (out / (cell["key"] + ".log")).write_text(result.stdout + result.stderr)
        rot, shift, _, _ = paired.base.alignment(meta[q], meta[r])
        ref = np.array(meta[q]["mapped"]) @ rot + shift
        branch = {
            int(i): (float(a), float(b), int(ok), int(n))
            for i, a, b, ok, n in re.findall(
                r"Silent branch arm=(\d) initial=(\S+) final=(\S+) reference=\S+ accepted=(\d) evals=(\d+)",
                result.stdout,
            )
        }
        for m in re.finditer(
            r"Silent replay step=(\d+) stage=(\w+) source_rms=(\S+) max_displacement=(\S+)\n(MODEL[^$]*?ENDMDL)",
            result.stdout,
            re.DOTALL,
        ):
            step, stage, displacement, largest, block = m.groups()
            xyz = paired.atoms(block.splitlines())
            pose = np.array([xyz[i][0] for i in meta[q]["serial"]])
            rms = np.sqrt(np.square(ref - pose).sum(axis=2).mean(axis=1)).min()
            info = (
                branch[stage.startswith("redirected")]
                if stage.endswith("hunt")
                else [""] * 4
            )
            stages.append(
                dict(
                    key=cell["key"],
                    step=step,
                    stage=stage,
                    displacement_A=float(displacement),
                    max_displacement_A=float(largest),
                    crystal_rmsd_A=float(rms),
                    initial_energy=info[0],
                    final_energy=info[1],
                    accepted=info[2],
                    hunt_evals=info[3],
                )
            )
    paired.write(out / "stages.csv", stages)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("out", type=Path)
    for name in ("always", "adaptive"):
        parser.add_argument("--" + name, required=True)
    run(parser.parse_args())

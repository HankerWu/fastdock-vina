"""Prepare and redock the two predeclared ER-beta cases with an external Vina binary."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time


BENCH = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(BENCH / "pair_gradient/pose"))
import prepare

CASES = {"1QKM": ("ESR2", "GEN", "A"), "5TOA": ("ESR2", "EST", "A")}


def main(work, binary):
    os.environ.setdefault("OPENMM_CPU_THREADS", "1")
    prepare.CASES.update(CASES)
    for code in CASES:
        if not binary:
            start = time.perf_counter()
            prepare.prepare(work, code)
            (work / code / "prep_seconds.txt").write_text(f"{time.perf_counter()-start:.6f}\n")
            continue
        site = work / code
        meta = json.loads((site / "preparation.json").read_text())
        if meta["receptor_status"] != "ready":
            continue
        for cap in (12000, 60000):
            for seed in (52711, 52727):
                name = f"{code}-{cap}-{seed}"
                output = site / (name + ".pdbqt")
                cmd = [str(binary), "--ligand", str(site / "ligand.pdbqt"), "--receptor", str(site / "receptor.pdbqt"),
                       "--cpu", "1", "--exhaustiveness", "4", "--num_modes", "9", "--energy_range", "100",
                       "--max_evals", str(cap), "--seed", str(seed), "--out", str(output)]
                for axis, center in zip("xyz", meta["center"]):
                    cmd.extend(["--center_" + axis, str(center), "--size_" + axis, "24"])
                start = time.perf_counter()
                with (site / (name + ".log")).open("w") as log:
                    subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
                (site / (name + ".json")).write_text(json.dumps(dict(
                    cmd=[Path(arg).name if "/" in arg else arg for arg in cmd],
                    wall_s=time.perf_counter()-start,
                    binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                    pool_sha256=hashlib.sha256(output.read_bytes()).hexdigest()), indent=2)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("--binary", type=Path)
    args = parser.parse_args()
    main(args.work.resolve(), args.binary.resolve() if args.binary else None)

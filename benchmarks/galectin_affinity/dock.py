"""Execute a completely crossed, single-core paired search panel."""

import argparse
import hashlib
import itertools
import json
import random
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

STATES = ("1GZW", "4Q26", "1KJL", "3T1L")
SEEDS = (271828, 314159, 161803, 141421)
ARMS = ("native", "balanced", "anchored", "balanced_anchored", "prepared")


def run(work, binary, jobs, arms):
    dest = work / "pools"
    dest.mkdir(exist_ok=True)
    binary_hash = hashlib.sha256(binary.read_bytes()).hexdigest()
    cases = list(itertools.product(STATES, range(5, 18), SEEDS, arms))
    random.Random(2026).shuffle(cases)
    from sugar import crystal_core, ligand_core

    for state, compound in itertools.product(STATES, range(5, 18)):
        xyz = crystal_core(work, state)
        serials = ligand_core(work, compound)
        anchor = work / "prepared" / f"{state}_{compound}.anchor"
        anchor.write_text(
            "".join(
                f"{i} {x:.4f} {y:.4f} {z:.4f}\n" for i, (x, y, z) in zip(serials, xyz)
            )
        )

    def one(case):
        state, compound, seed, arm = case
        name = f"{state}_{compound}_{seed}_{arm}"
        out, log, meta = [dest / (name + s) for s in (".pdbqt", ".log", ".json")]
        args = [
            binary,
            "--receptor",
            work / "prepared" / f"{state}.pdbqt",
            "--ligand",
            work / "prepared" / f"{compound}.pdbqt",
            "--config",
            work / "prepared" / f"{state}.box",
            "--seed",
            str(seed),
            "--cpu",
            "1",
            "--exhaustiveness",
            "2",
            "--max_evals",
            "25000",
            "--pose_selection",
            "late",
            "--num_modes",
            "9",
            "--energy_range",
            "100",
            "--verbosity",
            "2",
            "--out",
            out,
        ]
        if arm in ("balanced", "balanced_anchored"):
            args.append("--balanced_torsion_search")
        if arm in ("anchored", "balanced_anchored"):
            args.extend(
                ["--search_anchor", work / "prepared" / f"{state}_{compound}.anchor"]
            )
        if arm == "prepared":
            args.append("--prepared_torsion_starts")
        start = time.monotonic()
        with log.open("w") as outf:
            subprocess.run(args, stdout=outf, stderr=subprocess.STDOUT, check=True)
        seconds = time.monotonic() - start
        data = out.read_bytes()
        if not data.endswith(b"ENDMDL\n"):
            raise ValueError(f"Incomplete docking: {name}")
        meta.write_text(
            json.dumps(
                {
                    "binary_sha256": binary_hash,
                    "out_sha256": hashlib.sha256(data).hexdigest(),
                    "wall_seconds": seconds,
                }
            )
            + "\n"
        )

    with ThreadPoolExecutor(max_workers=jobs) as pool:
        for _ in pool.map(one, cases):
            pass


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("binary", type=Path)
    p.add_argument("--jobs", type=int, default=3)
    p.add_argument("--arms", choices=ARMS, nargs="+", default=ARMS)
    args = p.parse_args()
    run(args.work, args.binary.resolve(), args.jobs, args.arms)

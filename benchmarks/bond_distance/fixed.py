"""Rescore every parent candidate using original topology and printed coordinates."""

import argparse
import hashlib
import re
import subprocess
import tempfile
import time
from pathlib import Path

from run import cases, read, write

TERMS = (
    "total",
    "inter",
    "ligand_receptor",
    "ligand_flex",
    "intra",
    "ligand",
    "flex_receptor",
    "flex_flex",
    "torsion",
    "unbound",
)


def main(a):
    cells = {
        (rec, cid): (protein, ligand, box)
        for rec, cid, protein, ligand, box in cases(a.work, a.study)
    }
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        pose = Path(tmp) / "pose.pdbqt"
        for r in read(a.output / "runs.csv"):
            if r["arm"] != "parent":
                continue
            protein, ligand, box = cells[r["pdb"], r["compound"]]
            pool = a.output / "raw" / (r["key"] + ".pdbqt")
            models = re.findall(r"^MODEL .*?^ENDMDL", pool.read_text(), re.M | re.S)
            assert len(models) == int(r["modes"])
            for rank, block in enumerate(models, 1):
                pose.write_text(block + "\n")
                for mode in ("common_zero", "self") if rank == 1 else ("common_zero",):
                    args = [
                        "--receptor",
                        protein,
                        "--ligand",
                        ligand,
                        "--config",
                        box,
                        "--pose_coordinates",
                        pose,
                        "--score_only",
                        "--cpu",
                        "1",
                        "--seed",
                        r["seed"],
                    ]
                    if mode == "common_zero":
                        args += ["--unbound_energy", "0"]
                    for arm, binary in (
                        ("parent", a.parent),
                        ("candidate", a.candidate),
                    ):
                        start = time.perf_counter_ns()
                        p = subprocess.run(
                            [binary, *args], capture_output=True, text=True
                        )
                        elapsed = time.perf_counter_ns() - start
                        values = re.findall(
                            r":\s*([-+.\deE]+)\s+\(kcal/mol\)", p.stdout
                        )
                        if p.returncode or len(values) != len(TERMS):
                            raise RuntimeError(
                                f"{r['key']}/{rank}/{arm}: {p.stdout} {p.stderr}"
                            )
                        rows.append(
                            dict(
                                key=r["key"],
                                rank=rank,
                                mode=mode,
                                arm=arm,
                                elapsed_ns=elapsed,
                                pose_sha256=hashlib.sha256(
                                    pose.read_bytes()
                                ).hexdigest(),
                                **dict(zip(TERMS, values)),
                            )
                        )
    write(a.output / "fixed.csv", rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("study", choices=("kd", "macro"))
    for arg in ("work", "output", "parent", "candidate"):
        p.add_argument(arg, type=Path)
    main(p.parse_args())

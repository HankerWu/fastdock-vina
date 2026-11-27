"""Move every newly appended atom outside the box; replay all eight blind runs."""

import argparse
import csv
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import replay


def run(work, binary, out):
    expected = "6a56ec3db8c23b2e8bea4f65eb94fe3a4ccc2b272b4ae27c0da8113c078b4304"
    if replay.digest(binary) != expected:
        raise ValueError(
            "This exact-replay control requires the recorded production executable."
        )
    out.mkdir(parents=True, exist_ok=True)
    raw = out / "raw_far"
    raw.mkdir(exist_ok=True)
    dry = work / replay.PROTEINS["dry"]
    released = work / replay.PROTEINS["1OL7_released"]
    src = released.read_bytes()
    assert src.startswith(dry.read_bytes())
    suffix = src[len(dry.read_bytes()) :].decode().splitlines()
    assert suffix and all(l.startswith(("ATOM  ", "HETATM")) for l in suffix)
    lines, moved = [], []
    for line in suffix:
        xyz = np.array(replay.geom.xyz(line)) + np.array([100.0, 0.0, 0.0])
        lines.append(line[:30] + "".join(f"{v:8.3f}" for v in xyz) + line[54:])
        moved.append(xyz)
    rec = raw / "far.pdbqt"
    rec.write_bytes(dry.read_bytes() + ("\n".join(lines) + "\n").encode())
    center = json.loads((work / "prepared/2X6E.json").read_text())["center"]
    # This lower bound to the entire 26-A cube includes all grid-interior points.
    distance = np.linalg.norm(
        np.maximum(np.abs(np.array(moved) - center) - 13, 0), axis=1
    ).min()
    assert distance > 8.5
    box = [
        s
        for c, v in zip("xyz", center)
        for s in ("--center_" + c, v, "--size_" + c, 26)
    ]
    cells = []
    with (work / "base/cells.csv").open() as f:
        original = {
            (r["drug"], int(r["seed"])): r
            for r in csv.DictReader(f)
            if r["arm"] == "dry"
        }
    for (drug, seed), row in original.items():
        dst = raw / f"far_{drug}_{seed}.pdbqt"
        dry_out = work / "base/blind" / f"dry_{drug}_{seed}.pdbqt"
        assert replay.digest(dry_out) == row["pool_sha256"]
        log = dst.with_suffix(".log")
        cost = replay.invoke(
            [
                binary,
                "--receptor",
                rec,
                "--ligand",
                work / "prepared" / f"{drug}.pdbqt",
                *box,
                "--seed",
                seed,
                "--exhaustiveness",
                2,
                "--cpu",
                1,
                "--num_modes",
                9,
                "--energy_range",
                100,
                "--max_evals",
                25000,
                "--out",
                dst,
                "--pose_selection",
                "late",
                "--search_counters",
            ],
            log,
        )
        counters = sum(
            map(
                int,
                re.findall(r"\b(?:hunt|second|explicit)_evals=(\d+)", log.read_text()),
            )
        )
        cells.append(
            dict(
                drug=drug,
                seed=seed,
                dry_sha256=row["pool_sha256"],
                far_sha256=replay.digest(dst),
                identical=dst.read_bytes() == dry_out.read_bytes(),
                dry_evals=row["evals"],
                far_evals=counters,
                **cost,
            )
        )
    replay.geom.write_csv(out, "far.csv", cells)
    meta = {
        "binary_sha256": expected,
        "dry_sha256": replay.digest(dry),
        "released_sha256": replay.digest(released),
        "far_receptor_sha256": replay.digest(rec),
        "n_added": len(suffix),
        "n_added_heavy": sum(l.split()[-1] not in ("H", "HD") for l in suffix),
        "min_to_cube_A": float(distance),
        "displacement_A": [100, 0, 0],
    }
    (out / "far.json").write_text(json.dumps(meta, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for a in ("work", "binary", "out"):
        p.add_argument(a, type=Path)
    a = p.parse_args()
    run(a.work.resolve(), a.binary.resolve(), a.out.resolve())

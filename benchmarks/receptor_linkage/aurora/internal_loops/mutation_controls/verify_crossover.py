"""Verify every common-start graph intervention against the crystal and diagonal trace."""

import argparse
import json
import struct
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import check_predictions as audit


def check(work, crossed, traces):
    mol, query_order = audit.frame.read_query(work, "at9283")
    ref, _ = audit.refmol(work, "at9283")
    mappings = ref.GetSubstructMatches(mol, uniquify=False, maxMatches=10000)
    coords = ref.GetConformer().GetPositions()[np.array(mappings)]
    report = []
    total = 0
    for source in ("inferred", "complete"):
        arms = {
            target: np.load(crossed / f"{source}_on_{target}.npz")
            for target in ("inferred", "complete")
        }
        diagonal = arms[source]
        with (traces / f"at9283_{source}_61231.bin").open("rb") as f:
            (n,) = struct.unpack("<I", f.read(4))
            serials = np.fromfile(f, "<u4", n)
        assert len(serials) == len(query_order) == 28
        dtype = np.dtype(
            [
                ("tr", "<u4"),
                ("step", "<u4"),
                ("phase", "<u4"),
                ("event", "<u4"),
                ("local", "<u4"),
                ("e", "<f8"),
                ("xyz", "<f8", (n, 3)),
            ]
        )
        trace = np.memmap(
            traces / f"at9283_{source}_61231.bin",
            dtype=dtype,
            mode="r",
            offset=4 * (n + 1),
        )
        flag = (trace["event"] == 2) & (trace["phase"] == 1) & (trace["local"] == 1)
        starts, ends = trace[flag], trace[trace["event"] == 4]
        orig_order = [serials.tolist().index(s) for s in diagonal["serials"]]
        for a, b, xyz in (
            (diagonal["start_xyz"], starts, "proposal"),
            (diagonal["return_xyz"], ends, "local_return"),
        ):
            np.testing.assert_array_equal(a, b["xyz"][:, orig_order])
            np.testing.assert_array_equal(b["tr"], diagonal["trajectory"])
            np.testing.assert_array_equal(b["step"], diagonal["step"])
        for key in ("serials", "generalized", "trajectory", "step", "start_xyz"):
            np.testing.assert_array_equal(arms["inferred"][key], arms["complete"][key])
        for target, data in arms.items():
            ref_xyz = coords[:, [query_order.index(s) for s in data["serials"]]]

            def rms(x):
                return np.min(
                    [
                        np.sqrt(np.mean(np.sum((x - y) ** 2, axis=-1), axis=-1))
                        for y in ref_xyz
                    ],
                    axis=0,
                )

            initial = rms(data["start_xyz"])
            out = rms(data["return_xyz"])
            error = max(
                np.max(np.abs(initial - data["start_rmsd"])),
                np.max(np.abs(out - data["end_rmsd"])),
            )
            assert error < 1e-8
            near = np.flatnonzero(out <= 2)
            best = int(out.argmin())
            report.append(
                dict(
                    source=source,
                    target=target,
                    starts=len(initial),
                    best_return_A=float(out[best]),
                    best_trajectory=int(data["trajectory"][best]),
                    best_step=int(data["step"][best]),
                    exact_diagonal_coordinate_replay=(source == target),
                    max_RMSD_difference_A=error,
                    start_le2=int(np.count_nonzero(initial <= 2)),
                    return_le2=[
                        dict(
                            trajectory=int(data["trajectory"][i]),
                            step=int(data["step"][i]),
                            proposal_A=float(initial[i]),
                            local_return_A=float(out[i]),
                        )
                        for i in near
                    ],
                )
            )
            total += len(out)
    print(json.dumps(dict(crossed_local_minimizations=total, arms=report), indent=2))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("work", "crossed", "traces"):
        p.add_argument(k, type=Path)
    check(*(v.resolve() for v in vars(p.parse_args()).values()))

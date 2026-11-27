"""Independently check MC parent lineage, unlimited-pivot minima and the 2 Å classifications."""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent), str(HERE.parent / "basin_trace")]
from replay_basin import read_trace
import check_predictions as audit


def constrained_error(a, b):
    # For an energy E and move length squared L this is
    # min_R [E(R) + lambda*(L(R)-4)].  align_vectors uses column rotations.
    def dual(mult, details=False):
        rot, _ = Rotation.align_vectors(b + mult * a, a)
        trial = rot.apply(a)
        err = ((trial - b) ** 2).sum(1).mean()
        length = ((trial - a) ** 2).sum(1).mean()
        return (err, length) if details else err + mult * (length - 4.0)

    err, length = dual(0.0, True)
    if err > 4 + 1e-9:
        return False
    if length <= 4:
        return True
    hi = 1.0
    while dual(hi, True)[1] > 4:
        hi *= 2
    optimum = minimize_scalar(
        lambda v: -dual(v), bounds=(0, hi), method="bounded", options={"xatol": 1e-12}
    )
    if -optimum.fun > 4 + 1e-9:
        return False
    for extra in (1.0e-10, 1.0e-8, 1.0e-6, 1.0e-4):
        err, length = dual(optimum.x + extra, True)
        if err <= 4 and length <= 4:
            return True
    raise RuntimeError("Threshold not certified")


def verify(work, streams, results):
    num = near = checks = 0
    worst = 0.0
    with np.load(results / "parents.npz") as archive:
        for case in archive:
            drug = case.split("_")[0]
            serials, trace = read_trace(streams / (case + ".bin"))
            rows = archive[case]
            ref, _ = audit.refmol(work, drug)
            query, idx = audit.frame.read_query(work, drug)
            refs = ref.GetConformer().GetPositions()[
                np.array(ref.GetSubstructMatches(query, uniquify=False))[
                    :, [idx.index(s) for s in serials]
                ]
            ]
            history = []
            for i in np.flatnonzero(np.isin(trace["event"], (1, 4, 5, 6))):
                if trace["event"][i] == 4:
                    history.append(live)
                else:
                    live = i
            np.testing.assert_array_equal(rows["parent_event"], history)
            num += len(rows)
            # Check the strongest unlimited example for each trajectory, directly in Cartesian coordinates.
            for tr in (1, 2):
                ids = np.flatnonzero(trace["traj"][history] == tr)
                r = rows[ids[np.argmin(rows["any_pivot_A"][ids])]]
                x = trace["xyz"][r["parent_event"]]
                best = np.inf
                for y in refs:
                    for p in x:
                        rot, _ = Rotation.align_vectors(y - p, x - p)
                        fit = rot.apply(x - p) + p
                        best = min(best, np.sqrt(((fit - y) ** 2).sum(1).mean()))
                worst = max(worst, abs(r["any_pivot_A"] - best))
                checks += 1
            # All rows too far by the RMS triangle inequality cannot yield a direct hit.
            assert not np.any(rows["bounded_pivot_le2"][rows["raw_A"] > 4])
            for i in np.unique(rows["parent_event"][rows["raw_A"] <= 4]):
                x = trace["xyz"][i]
                expected = rows["bounded_pivot_le2"][rows["parent_event"] == i][0]
                actual = any(constrained_error(x - p, y - p) for y in refs for p in x)
                assert expected == actual
                near += 1
    assert worst < 1e-8
    return dict(
        parent_visits=num,
        threshold_tests_for_parent_events_within_4=near,
        unlimited_witnesses=checks,
        max_independent_witness_error_A=worst,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("work", "streams", "results"):
        parser.add_argument(key, type=Path)
    print(json.dumps(verify(**vars(parser.parse_args())), indent=2))

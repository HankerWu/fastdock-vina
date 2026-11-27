"""Crystal-informed geometric reach of one rigid proposal from each actual MC parent."""

import argparse
import csv
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent), str(HERE.parent / "basin_trace")]
import check_predictions as audit
from replay_basin import read_trace


def squared_norm(x):
    return np.einsum("...i,...i->...", x, x)


def optimum(cross):
    """Maximize <R,cross> over proper row-vector rotations."""
    u, singular, vh = np.linalg.svd(cross)
    sign = np.linalg.det(u @ vh)
    u[..., :, -1] *= sign[..., None]
    value = singular[..., 0] + singular[..., 1] + sign * singular[..., 2]
    return value, u @ vh


def constrained(x, y, radius):
    """Primal witness and Lagrangian lower bound for one pivot and atom map."""
    xx, xy = x.T @ x / len(x), x.T @ y / len(x)
    tr = np.trace(xx)
    total = tr + squared_norm(y).mean()
    orig = squared_norm(x - y).mean()
    lower, upper = 0.0, orig
    lo, hi = 0.0, 1.0

    def evaluate(mult):
        peak, r = optimum(xy + mult * xx)
        dual = total + mult * (2 * tr - radius**2) - 2 * peak
        move2 = squared_norm(x @ r - x).mean()
        err2 = squared_norm(x @ r - y).mean()
        return dual, move2, err2

    while evaluate(hi)[1] > radius**2:
        hi *= 2
    for _ in range(56):
        mid = (lo + hi) / 2
        dual, move2, err2 = evaluate(mid)
        lower = max(lower, dual)
        if move2 <= radius**2:
            hi, upper = mid, min(upper, err2)
        else:
            lo = mid
    return max(0.0, lower), upper


def parent_indices(trace):
    """MC tmp changes on accept and again on second minimization in the default chain."""
    indices = []
    event = trace["event"]
    for tr in (1, 2):
        in_tr = trace["traj"] == tr
        starts = np.flatnonzero(
            in_tr & (trace["phase"] == 1) & (event == 2) & (trace["local"] == 1)
        )
        updates = np.flatnonzero(in_tr & np.isin(event, (1, 5, 6)))
        assert event[updates[0]] == 1
        assert np.array_equal(trace["step"][starts], np.arange(1, len(starts) + 1))
        indices.extend(updates[np.searchsorted(updates, starts) - 1])
    return np.array(indices)


def measure(x, references):
    """Visit every heavy-atom pivot and symmetry, then resolve only ambiguous 2 A decisions."""
    n, atoms, _ = x.shape
    cx = x.mean(1)
    px = cx[:, None] - x
    dx = x - cx[:, None]
    xx0 = squared_norm(dx).mean(1)
    xx = xx0[:, None] + squared_norm(px)
    current = np.full(n, np.inf)
    free = current.copy()
    pivot = current.copy()
    for y in references:
        cy = y.mean(0)
        dy = y - cy
        current = np.minimum(current, squared_norm(x - y).mean(1))
        py = cy - x
        yy0 = squared_norm(dy).mean()
        yy = yy0 + squared_norm(py)
        cov = dx.transpose(0, 2, 1) @ dy / atoms
        val, _ = optimum(cov)
        free = np.minimum(free, xx0 + yy0 - 2 * val)
        joint = cov[:, None] + px[:, :, :, None] * py[:, :, None, :]
        peaks, _ = optimum(joint)
        pivot = np.minimum(pivot, np.min(xx + yy - 2 * peaks, axis=1))
    current, free, pivot = (np.sqrt(np.maximum(0, t)) for t in (current, free, pivot))
    assert np.all(free <= pivot + 1e-8) and np.all(pivot <= current + 1e-8)
    # Quotient RMSD obeys the triangle inequality. If D-2>2 even an arbitrary
    # internal change with this RMS displacement cannot hit the geometric target.
    hit = current <= 2
    possible = np.flatnonzero((~hit) & (current <= 4 + 1e-9) & (pivot <= 2 + 1e-9))
    for i in possible:
        best_lo, best_hi = np.inf, np.inf
        for y in references:
            for anchor in x[i]:
                a, b = x[i] - anchor, y - anchor
                _, rot = optimum(a.T @ b / atoms)
                err2 = squared_norm(a @ rot - b).mean()
                if err2 > 4 + 1e-9:
                    lower, upper = err2, np.inf
                elif squared_norm(a @ rot - a).mean() <= 4:
                    lower = upper = err2
                else:
                    lower, upper = constrained(a, b, 2.0)
                best_lo, best_hi = min(best_lo, lower), min(best_hi, upper)
                if best_hi <= 4:
                    break
            if best_hi <= 4:
                break
        assert best_hi <= 4 or best_lo > 4 + 1e-9, (best_lo, best_hi)
        hit[i] = best_hi <= 4
    return current, free, pivot, hit


def main(work, streams, dest):
    dest.mkdir(parents=True, exist_ok=True)
    start, cpu = time.perf_counter(), time.process_time()
    meta = audit.rows(HERE.parent / "basin_trace/results/runs.csv")
    counts, witnesses, arrays = [], [], {}
    for row in meta:
        name, graph, seed = (row[k] for k in ("drug", "graph", "seed"))
        case = f"{name}_{graph}_{seed}"
        stream = streams / f"{case}.bin"
        with stream.open("rb") as f:
            assert hashlib.file_digest(f, "sha256").hexdigest() == row["raw_sha256"]
        serials, tr = read_trace(stream)
        mol, order = audit.frame.read_query(work, name)
        ref, _ = audit.refmol(work, name)
        mappings = np.array(
            ref.GetSubstructMatches(mol, uniquify=False, maxMatches=10000)
        )
        assert len(mappings) == {"at9283": 4, "danusertib": 8}[name]
        refs = ref.GetConformer().GetPositions()[
            mappings[:, [order.index(s) for s in serials]]
        ]
        idx = parent_indices(tr)
        idx_unique, inv = np.unique(idx, return_inverse=True)
        raw, free, pivot, hit = measure(tr["xyz"][idx_unique], refs)
        with np.load(HERE.parent / "basin_trace/results" / (case + ".npz")) as previous:
            np.testing.assert_allclose(
                raw, previous["rmsd_A"][idx_unique], rtol=0, atol=2e-10
            )
        arrays[case] = np.rec.fromarrays(
            [idx, raw[inv], free[inv], pivot[inv], hit[inv]],
            names="parent_event,raw_A,any_rigid_A,any_pivot_A,bounded_pivot_le2",
        )
        for k in range(2):
            rows = np.flatnonzero(tr["traj"][idx] == k + 1)
            for weighting, take in (
                ("step_visits", rows),
                (
                    "distinct_parent_events",
                    rows[np.unique(idx[rows], return_index=True)[1]],
                ),
            ):
                jj = inv[take]
                counts.append(
                    dict(
                        case=case,
                        trajectory=k + 1,
                        weighting=weighting,
                        n=len(jj),
                        min_parent_A=min(raw[jj]),
                        min_free_A=min(free[jj]),
                        min_pivot_A=min(pivot[jj]),
                        raw_le2=int(sum(raw[jj] <= 2)),
                        free_le2=int(sum(free[jj] <= 2)),
                        unlimited_pivot_le2=int(sum(pivot[jj] <= 2)),
                        bounded_pivot_le2=int(sum(hit[jj])),
                        within_4=int(sum(raw[jj] <= 4)),
                    )
                )
            j = rows[np.argmin(pivot[inv[rows]])]
            i = inv[j]
            witnesses.append(
                dict(
                    case=case,
                    trajectory=k + 1,
                    next_step=int(j - rows[0] + 1),
                    parent_event=int(idx[j]),
                    parent_source_step=int(tr["step"][idx[j]]),
                    parent_source_event=int(tr["event"][idx[j]]),
                    raw_A=raw[i],
                    free_A=free[i],
                    unlimited_pivot_A=pivot[i],
                    bounded_pivot_le2=bool(hit[i]),
                )
            )
    for name, rows in (("parent_counts", counts), ("witnesses", witnesses)):
        with (dest / (name + ".csv")).open("w") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
    np.savez_compressed(dest / "parents.npz", **arrays)
    print(
        json.dumps(
            dict(
                wall_s=time.perf_counter() - start,
                cpu_s=time.process_time() - cpu,
                visited=sum(len(v) for v in arrays.values()),
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("work", "streams", "dest"):
        p.add_argument(arg, type=Path)
    main(**vars(p.parse_args()))

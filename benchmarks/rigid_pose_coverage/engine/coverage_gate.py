"""Pre-registered whole-ligand gates for rigid receptor cavity correspondences.

The prediction path never reads a ligand reference or assay. Positions use angstroms
and float64. Original pocket.triangles row numbers are retained in every bank.
"""
from pathlib import Path
import sys
from time import perf_counter, process_time

import numpy as np
from scipy.spatial.distance import cdist, pdist

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools"))
from pocket_triangle import _ligand_triangle

ARMS = ("edge", "feasible_edge", "coverage", "uniform", "uniform_work")
BUDGET = 1024


def _geometry(pocket, xyz, serials):
    x, s = np.asarray(xyz, dtype=np.float64), np.asarray(serials, dtype=int)
    if x.shape != (len(s), 3) or len(s) != len(np.unique(s)) or not np.isfinite(x).all():
        raise ValueError("expected ordered, finite source coordinates and unique heavy atom serials")
    order = np.argsort(s)
    anchor = order[_ligand_triangle(x[order])]
    residual = ((pocket.edges - pdist(x[anchor])) ** 2).sum(axis=1)
    return x, anchor, residual


def candidates(pocket, xyz, serials):
    """Fit all legal ordered triples. Return arrays indexed by original triangle ID."""
    x, anchor, residual = _geometry(pocket, xyz, serials)
    center, rot = pocket.fit_triangles(x, x[anchor], pocket.triangles)
    centered = x - x.mean(axis=0)
    lo, hi = pocket.center - pocket.half_width, pocket.center + pocket.half_width
    eligible = np.empty(len(center), dtype=bool)
    # Bound temporary storage for long ligands. All atoms, including distal ones,
    # must be inside the requested closed box. No energy, crystal or pH criterion.
    for i in range(0, len(center), 512):
        take = slice(i, i + 512)
        coords = centered @ rot[take].transpose(0, 2, 1) + center[take, None, :]
        eligible[take] = ((lo <= coords) & (coords <= hi)).all(axis=(1, 2))
    return dict(triangle_id=np.arange(len(center)), center=center, matched_rotation=rot,
                site_triangles=pocket.triangles, ligand_triangle=anchor,
                edge_residual_A2=residual, eligible=eligible)


def embedding(xyz, centers, rotations):
    """Squared Euclidean row distance equals mean squared *labeled* atom displacement."""
    x = np.asarray(xyz, dtype=np.float64)
    x = x - x.mean(axis=0)
    moment = x.T @ x / len(x)
    values, axes = np.linalg.eigh(moment)
    # Moment is a Gram matrix. Planar or collinear sources can have zero eigenvalues.
    factor = axes * np.sqrt(np.maximum(values, 0))
    rotational = np.asarray(rotations) @ factor
    return np.ascontiguousarray(np.column_stack((centers, rotational.reshape(len(centers), 9))))


def _farthest(features, first, budget):
    """Return local indices and achieved squared covering radius over the finite domain."""
    chosen = np.empty(budget, dtype=int)
    nearest = np.full(len(features), np.inf)
    current = first
    for j in range(budget):
        chosen[j] = current
        squared = cdist(features, features[current:current + 1], metric="sqeuclidean")[:, 0]
        np.minimum(nearest, squared, out=nearest)
        nearest[current] = -np.inf  # distinct correspondences even at zero distance
        if j + 1 < budget:
            current = int(np.argmax(np.round(nearest, 10)))
    return chosen, max(0., float(np.max(nearest)))


def select(pool, xyz, seed, arm, budget=BUDGET):
    """Return original triangle IDs in screening order and squared coverage radius.

    The radius is reported for arms executing the farthest-first loop. For uniform_work
    it describes the *discarded* coverage bank, never the uniform gate's performance.
    """
    if arm not in ARMS or budget < 1:
        raise ValueError("unknown gate or invalid budget")
    ids = pool["triangle_id"] if arm == "edge" else np.flatnonzero(pool["eligible"])
    if len(ids) < budget:
        raise ValueError(f"{arm}: {len(ids)} available correspondences, require {budget}")
    if arm in ("edge", "feasible_edge"):
        priority = np.argsort(np.round(pool["edge_residual_A2"][ids], 10), kind="stable")
        return ids[priority[:budget]], np.nan
    radius = np.nan
    if arm in ("coverage", "uniform_work"):
        first = int(np.argmin(np.round(pool["edge_residual_A2"][ids], 10)))
        features = embedding(xyz, pool["center"][ids], pool["matched_rotation"][ids])
        chosen, radius = _farthest(features, first, budget)
        if arm == "coverage":
            return np.sort(ids[chosen]), radius
    rng = np.random.Generator(np.random.PCG64(seed))
    return np.sort(rng.choice(ids, size=budget, replace=False)), radius


def place(pocket, xyz, serials, ligand, bonds, seed, arm):
    """Adapter for the native FIFO search controller and the AD4 replay controller."""
    wall, cpu = perf_counter(), process_time()
    x, serials = np.asarray(xyz, dtype=np.float64), np.asarray(serials)
    if arm == "edge":
        # Preserve legacy matrix batch and summation arithmetic, and pay only for
        # its original 1,024 fits. Feasibility would require additional fits.
        x, anchor, residual = _geometry(pocket, x, serials)
        ids, radius = select(dict(triangle_id=np.arange(len(residual)), edge_residual_A2=residual),
                             x, seed, arm)
        centers, rot = pocket.fit_triangles(x, x[anchor], pocket.triangles[ids])
        feasible_count, fitted_count = -1, len(ids)
    else:
        pool = candidates(pocket, x, serials)
        ids, radius = select(pool, x, seed, arm)
        centers, rot = pool["center"][ids], pool["matched_rotation"][ids]
        anchor, residual = pool["ligand_triangle"], pool["edge_residual_A2"]
        feasible_count, fitted_count = int(pool["eligible"].sum()), len(residual)
    gate_wall_s, gate_cpu_s = perf_counter() - wall, process_time() - cpu
    wall, cpu = perf_counter(), process_time()
    poses = (x - x.mean(axis=0)) @ rot.transpose(0, 2, 1) + centers[:, None, :]
    shape = pocket._score(poses)
    return dict(center=centers, ligand_triangle=anchor, site_triangles=pocket.triangles[ids],
                triangle_id=ids, edge_residual_A2=residual[ids], matched_rotation=rot,
                matched_score=shape, matched_best=int(np.argmin(np.round(shape, 10))),
                # Legacy schema, only used to compare banks with their parent.
                shift=int(np.random.default_rng(seed).integers(1, len(ids))),
                anchor_serials=serials[anchor], component_serials=np.sort(serials), anchor_fallback=False,
                fitted_count=fitted_count, feasible_count=feasible_count,
                distance_evals=feasible_count * len(ids) if arm in ("coverage", "uniform_work") else 0,
                radius_A2=radius, gate_wall_s=gate_wall_s, gate_cpu_s=gate_cpu_s,
                shape_queries=len(ids)*len(x), shape_wall_s=perf_counter()-wall, shape_cpu_s=process_time()-cpu)

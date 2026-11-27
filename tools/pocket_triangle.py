"""Blind rigid placement by receptor cavity triplets. Coordinates are in angstroms."""

from itertools import permutations
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from scipy.spatial import cKDTree, distance


def _ligand_triangle(x):
    """Most separated atom pair, then the atom farthest from their line."""
    a, b = np.triu_indices(len(x), k=1)
    if len(x) < 3:
        raise ValueError("three noncollinear heavy atoms are required")
    d = np.linalg.norm(x[b] - x[a], axis=1)
    i = np.argmax(np.round(d, 10))
    if d[i] < 1e-8:
        raise ValueError("all ligand heavy atoms coincide")
    direction = (x[b[i]] - x[a[i]]) / d[i]
    heights = np.linalg.norm(np.cross(x - x[a[i]], direction), axis=1)
    c = np.argmax(np.round(heights, 10))
    if heights[c] * d[i] < 1e-8:
        raise ValueError("ligand triangle is collinear")
    return np.array([a[i], b[i], c])


class Pocket:
    """Cavity samples, their ordered triangles and a heavy-atom distance field."""

    @classmethod
    def from_pdbqt(cls, receptor, center, box_size):
        xyz = [
            [float(l[k : k + 8]) for k in (30, 38, 46)]
            for l in Path(receptor).read_text().splitlines()
            if l.startswith(("ATOM  ", "HETATM"))
            and l.split()[-1] not in ("H", "HD", "HS")
        ]
        return cls(xyz, center, box_size)

    def __init__(self, receptor_xyz, center, box_size):
        self.tree = cKDTree(receptor_xyz)
        self.center = np.array(center, dtype=float)
        self.half_width = np.broadcast_to(np.asarray(box_size) / 2, (3,))
        half_cells = np.floor(self.half_width / 1.5).astype(int)
        axes = [np.arange(-h, h + 1) * 1.5 for h in half_cells]
        grid = self.center + np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)
        nearest = self.tree.query(grid)[0]
        eligible = (nearest >= 3) & (nearest <= 5.5)
        points = grid[eligible]
        depth = self.tree.query_ball_point(points, 8.0, return_length=True)
        sort_key = np.round(abs(nearest[eligible] - 3.8), 10)
        order = np.lexsort((np.arange(len(points)), sort_key, -depth))
        chosen = []
        for i in order:
            if all(np.linalg.norm(points[i] - points[j]) >= 2.5 for j in chosen):
                chosen.append(i)
                if len(chosen) == 32:
                    break
        self.n_grid = len(grid)
        self.n_burial = len(points)
        self.sites = points[chosen]
        triangles = np.asarray(list(permutations(range(len(self.sites)), 3)), dtype=int).reshape(-1, 3)
        if len(triangles) == 0:
            raise ValueError("fewer than three cavity sites")
        a, b, c = self.sites[triangles].transpose(1, 0, 2)
        self.triangles = triangles[np.linalg.norm(np.cross(b - a, c - a), axis=1) >= 1e-8]
        if not len(self.triangles):
            raise ValueError("cavity sites are collinear")
        d = distance.cdist(self.sites, self.sites)
        a, b, c = self.triangles.T
        self.edges = np.stack((d[a, b], d[a, c], d[b, c]), axis=1)

    def _score(self, xyz):
        d = self.tree.query(xyz.reshape(-1, 3))[0].reshape(xyz.shape[:-1])
        u = d - 3.8
        violation = np.maximum(abs(xyz - self.center) - self.half_width, 0)
        return (np.maximum(-u, 0) ** 2 - np.exp(-u * u)).sum(axis=1) + (
            violation * violation
        ).sum(axis=(1, 2))

    def place(self, xyz, seed, timings=None, arms=("matched", "sham")):
        """Return all blind candidates. Apply each matrix to a centered column vector."""
        start = perf_counter(), process_time()
        x = np.asarray(xyz, dtype=float)
        atoms = _ligand_triangle(x)
        src = x[atoms]
        residual = ((self.edges - distance.pdist(src)) ** 2).sum(axis=1)
        priority = np.argsort(np.round(residual, 10), kind="stable")[:1024]
        triples = self.triangles[priority]
        dest = self.sites[triples]
        if len(dest) < 2:
            raise ValueError("cannot form a permutation sham with fewer than two candidates")
        dest_center = dest.mean(axis=1)
        src_center = src.mean(axis=0)
        cov = (src - src_center).T @ (dest - dest_center[:, None, :])
        u, _, vt = np.linalg.svd(cov)
        u[:, :, -1] *= np.linalg.det(u @ vt)[:, None]
        matched = (u @ vt).transpose(0, 2, 1)
        centers = (matched @ (x.mean(axis=0) - src_center)) + dest_center
        shift = int(np.random.default_rng(seed).integers(1, len(dest)))
        sham = np.roll(matched, shift, axis=0)
        centered = x - x.mean(axis=0)
        out = dict(
            center=centers, ligand_triangle=atoms, site_triangles=triples,
            edge_residual_A2=residual[priority], shift=shift,
        )
        if timings is not None:
            timings["match_wall_s"] = perf_counter() - start[0]
            timings["match_cpu_s"] = process_time() - start[1]
        rotations = {"matched": matched, "sham": sham}
        for arm in arms:
            rot = rotations[arm]
            start = perf_counter(), process_time()
            poses = centered @ rot.transpose(0, 2, 1) + centers[:, None, :]
            score = self._score(poses)
            out[arm + "_rotation"] = rot
            out[arm + "_score"] = score
            out[arm + "_best"] = int(np.argmin(np.round(score, 10)))
            if timings is not None:
                timings[arm + "_wall_s"] = perf_counter() - start[0]
                timings[arm + "_cpu_s"] = process_time() - start[1]
        return out

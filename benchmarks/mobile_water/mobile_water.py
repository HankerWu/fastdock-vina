"""Marginalize the translational position of an experimentally resolved site oxygen.

Distances are in Angstrom, isotropic B in Angstrom squared, and scores in kcal/mol.
This fixed-coordinate research rescorer uses a geometric site measure.
"""

from dataclasses import dataclass
from functools import lru_cache

import numpy as np

RT = 0.00198720425864083 * 298.15
PRESENT = 0.75


def accessible(distance):
    fraction = np.clip((distance - 2.0) / 0.6, 0.0, 1.0)
    return fraction**2 * (3.0 - 2.0 * fraction)


@lru_cache(None)
def gaussian_grid(order):
    """Cartesian product cubature for three independent standard normals."""
    if order < 1:
        raise ValueError("Quadrature order must be positive")
    x, w = np.polynomial.hermite.hermgauss(order)
    nodes = np.stack(
        np.meshgrid(x * np.sqrt(2), x * np.sqrt(2), x * np.sqrt(2), indexing="ij"),
        axis=-1,
    ).reshape(-1, 3)
    weights = np.einsum("i,j,k->ijk", w, w, w).ravel() / np.pi**1.5
    return nodes, weights


def penalty(availability, prior=PRESENT):
    if not 0 < prior < 1:
        raise ValueError("Occupied-state prior must be between zero and one")
    return -RT * np.log1p(-prior * (1.0 - np.clip(availability, 0.0, 1.0)))


@dataclass
class Site:
    center: np.ndarray
    samples: np.ndarray
    weights: np.ndarray
    radius: float

    @classmethod
    def prepare(cls, center, b_iso, receptor, order=11):
        """Receptor is a cKDTree of heavy atoms (or None in free space)."""
        if not np.isfinite(b_iso) or b_iso < 0:
            raise ValueError("Expected finite nonnegative B factor")
        center = np.asarray(center, dtype=float)
        nodes, weights = gaussian_grid(order)
        displacements = nodes * np.sqrt(b_iso / (8.0 * np.pi**2))
        samples = center + displacements
        if receptor is not None:
            weights = weights * accessible(receptor.query(samples)[0])
        if weights.sum() <= 0:
            raise ValueError("No allowed reference mass at hydration site")
        return cls(
            center,
            samples,
            weights / weights.sum(),
            float(np.linalg.norm(displacements, axis=1).max()),
        )

    def evaluate(self, ligand):
        """Return point and mobile penalties on an identical ligand KD tree."""
        distance = ligand.query(self.center)[0]
        fixed = penalty(accessible(distance))
        if distance >= self.radius + 2.6:
            return float(fixed), 0.0
        distances = ligand.query(self.samples)[0]
        mobile = penalty(np.dot(self.weights, accessible(distances)))
        return float(fixed), float(mobile)

"""One bounded rigid-body step on trilinear receptor fields before pruning poses.

All coordinates and map spacing are in Angstrom. Atomwise fields are in kcal/mol.
The heavy-atom centroid defines rotation and the heavy-atom RMS tangent metric
defines the trust radius. Flexible, explicit scoring is left to the caller.
"""

import numpy as np
from scipy.spatial.transform import Rotation


def field_values(field, xyz, origin, spacing, gradient=True):
    """Interpolate a single XS map at N positions, with analytic spatial derivatives.

    Outside the closed map domain use +inf (ineligible); the returned gradient
    is zero there. At internal cell faces use the derivative from the upper cell.
    """
    fractional = (xyz - origin) / spacing
    limit = np.asarray(field.shape) - 1
    inside = ((xyz >= origin) & (xyz <= origin + spacing * limit)).all(axis=1)
    bounded = np.clip(fractional, 0, limit)
    cell = np.minimum(bounded.astype(np.intp), limit - 1)
    x, y, z = cell.T
    a, b, c = (bounded - cell).T
    v000, v001 = field[x, y, z], field[x, y, z + 1]
    v010, v011 = field[x, y + 1, z], field[x, y + 1, z + 1]
    v100, v101 = field[x + 1, y, z], field[x + 1, y, z + 1]
    v110, v111 = field[x + 1, y + 1, z], field[x + 1, y + 1, z + 1]
    z00, z01 = (1 - c) * v000 + c * v001, (1 - c) * v010 + c * v011
    z10, z11 = (1 - c) * v100 + c * v101, (1 - c) * v110 + c * v111
    y0, y1 = (1 - b) * z00 + b * z01, (1 - b) * z10 + b * z11
    value = (1 - a) * y0 + a * y1
    value[~inside] = np.inf
    deriv = None
    if gradient:
        dz0 = (1 - b) * (v001 - v000) + b * (v011 - v010)
        dz1 = (1 - b) * (v101 - v100) + b * (v111 - v110)
        deriv = (
            np.column_stack(
                (
                    y1 - y0,
                    (1 - a) * (z01 - z00) + a * (z11 - z10),
                    (1 - a) * dz0 + a * dz1,
                )
            )
            / spacing
        )
        deriv[~inside] = 0
    return value, deriv


def rigid_trial(fields, origin, spacing, offsets, types, positions, arm, rng):
    """Evaluate parents and optionally a single metric-limited candidate per parent.

    offsets contain heavy atoms in one fixed orientation about their centroid.
    The returned translation and rotation vector describe the ACCEPTED move.
    """
    if arm not in ("raw", "gradient", "sham"):
        raise ValueError(f"unknown arm {arm}")
    n = len(offsets)
    if not n or len(types) != n or spacing <= 0:
        raise ValueError("expected nonempty typed heavy atoms and positive spacing")
    if not np.allclose(offsets.mean(axis=0), 0, atol=1e-8):
        raise ValueError("rotation offsets must be centered on the heavy atoms")
    force = np.zeros_like(positions, dtype=float)
    torque = np.zeros_like(force)
    parent = np.zeros(len(positions))
    for off, typ in zip(offsets, types):
        value, deriv = field_values(
            fields[typ], positions + off, origin, spacing, arm != "raw"
        )
        parent += value
        if deriv is not None:
            force -= deriv
            torque -= np.cross(off, deriv)
    if not np.all(np.isfinite(parent)):
        raise ValueError("a parent pose lies outside the map or has nonfinite energy")
    if arm == "raw":
        return dict(
            parent_energy=parent,
            energy=parent,
            trial_energy=parent,
            accepted=np.zeros(len(positions), dtype=bool),
            translation=force,
            rotvec=torque,
        )

    second = offsets.T @ offsets
    inertia = np.trace(second) * np.eye(3) - second
    eig, axes = np.linalg.eigh(inertia)
    observable = eig > 3 * np.finfo(float).eps * eig[-1]
    eig, axes = eig[observable], axes[:, observable]
    # For a unit-mass ligand the mean-square speed is |v|² + wᵀ I w / N.
    # In these whitened coordinates it is just the Euclidean squared norm.
    if arm == "gradient":
        t_dir = force / n
        w_dir = (torque @ axes) / np.sqrt(n * eig)
    else:
        direction = rng.standard_normal((len(positions), 3 + len(eig)))
        t_dir, w_dir = direction[:, :3], direction[:, 3:]
    norm = np.sqrt((t_dir * t_dir).sum(axis=1) + (w_dir * w_dir).sum(axis=1))
    scale = np.divide(spacing, norm, out=np.zeros_like(norm), where=norm > 0)
    translation = t_dir * scale[:, None]
    rotvec = ((w_dir * scale[:, None]) * np.sqrt(n / eig)) @ axes.T
    rotation = Rotation.from_rotvec(rotvec)
    trial = np.zeros_like(parent)
    center = positions + translation
    for off, typ in zip(offsets, types):
        value, _ = field_values(
            fields[typ], center + rotation.apply(off), origin, spacing, False
        )
        trial += value
    accepted = trial < parent
    translation[~accepted] = 0
    rotvec[~accepted] = 0
    return dict(
        parent_energy=parent,
        energy=np.minimum(trial, parent),
        trial_energy=trial,
        accepted=accepted,
        translation=translation,
        rotvec=rotvec,
    )

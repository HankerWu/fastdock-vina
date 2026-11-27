"""Check trilinear Cartesian forces and the rigid descent direction independently."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

from replay import atom_field, jacobian, tangent_step
from rigid_lookahead import field_values, rigid_trial


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("output", type=Path)
    args = p.parse_args()
    rng = np.random.default_rng(32131)
    field = rng.normal(size=(12, 13, 14))
    origin = np.array([2.1, -6.3, 7.1])
    h = 0.75
    # Interior points away from knots, edges and derivative ambiguities.
    indices = rng.integers(1, 10, size=(111, 3)) + rng.uniform(
        0.15, 0.85, size=(111, 3)
    )
    xyz = origin + h * indices
    en, grad = field_values(field, xyz, origin, h)
    eps = 1e-5
    fd = np.column_stack(
        [
            (
                atom_field(field, xyz + eps * d, origin, h)
                - atom_field(field, xyz - eps * d, origin, h)
            )
            / (2 * eps)
            for d in np.eye(3)
        ]
    )
    np.testing.assert_allclose(grad, fd, rtol=3e-8, atol=1e-8)
    metrics = {"max_grid_gradient_error": float(np.max(np.abs(grad - fd)))}
    np.testing.assert_allclose(en, atom_field(field, xyz, origin, h), atol=1e-13)
    # Both ends of the closed box and the adjacent representably-outside interval.
    indices = np.asarray(
        [[0, 0, 0], [11, 12, 13], [6.0, 4.2, 1.8], [-1e-4, 1, 1], [11 + 1e-4, 2, 2]]
    )
    e, g = field_values(field, origin + h * indices, origin, h)
    assert np.isfinite(e[:3]).all() and np.isposinf(e[3:]).all() and np.all(g[3:] == 0)
    # A Cartesian affine potential is exactly represented. Distinct slopes on atoms
    # produce a nontrivial torque independent of net translation.
    axes = np.moveaxis(np.indices((25, 25, 25), dtype=float), 0, -1) * h
    directions = np.array(
        [[1.0, -0.1, 0.4], [-0.5, 0.1, 1.2], [-0.4, 1.3, -0.8], [0.4, 0.3, 0.6]]
    )
    fields = {k: axes @ d for k, d in enumerate(directions)}
    center = np.array([[8.1, 7.8, 9.3]])
    zero = np.zeros(3)
    metrics["geometries"] = {}
    for name, off in {
        "point": [[0, 0, 0]],
        "line": [[-1, 0, 0], [1, 0, 0]],
        "nonlinear": [[-1, -0.5, 0], [1, -0.5, 0.4], [0, 1, -0.3], [1, 1, 1]],
    }.items():
        off = np.array(off, dtype=float)
        off -= off.mean(axis=0)
        types = list(range(len(off)))
        slopes = directions[: len(off)]
        dx = tangent_step(off, slopes[None, :, :], h)[0]
        jac = jacobian(off)
        tangent = (jac @ dx).reshape(-1, 3)
        realized = np.sqrt(np.square(tangent).sum(axis=1).mean())
        np.testing.assert_allclose(realized, h, atol=1e-13)
        test = rigid_trial(fields, zero, h, off, types, center, "gradient", rng)
        assert test["accepted"][0]
        actual = np.r_[test["translation"][0], test["rotvec"][0]]
        np.testing.assert_allclose(actual, dx, atol=1e-13)

        def energy(t):
            xyz = Rotation.from_rotvec(t * dx[3:]).apply(off) + center[0] + t * dx[:3]
            return sum(
                float(atom_field(fields[k], x[None, :], zero, h)[0])
                for k, x in zip(types, xyz)
            )

        predicted_slope = np.sum(slopes * tangent)
        finite_slope = (energy(eps) - energy(-eps)) / (2 * eps)
        np.testing.assert_allclose(finite_slope, predicted_slope, rtol=1e-7, atol=1e-8)
        assert predicted_slope < 0
        metrics["geometries"][name] = dict(
            rank=int(np.linalg.matrix_rank(jac)),
            tangent_rms=float(realized),
            slope_error=float(abs(finite_slope - predicted_slope)),
            kernel_step_error=float(np.max(np.abs(actual - dx))),
        )
    # Descent points directly out of a face. An unscorable trial must keep
    # the original finite pose, not a clamped or truncated conformation.
    field = -axes[..., 0]
    parents = np.array([[17.9, 4.0, 4.0]])
    out = rigid_trial(
        {"a": field}, zero, h, np.zeros((1, 3)), ["a"], parents, "gradient", rng
    )
    assert not out["accepted"][0] and np.isposinf(out["trial_energy"][0])
    np.testing.assert_array_equal(out["energy"], out["parent_energy"])
    np.testing.assert_array_equal(out["translation"], np.zeros((1, 3)))
    # A flat field has no descending move for either mechanism.
    for arm in ("gradient", "sham"):
        out = rigid_trial(
            {"a": np.ones((25, 25, 25))},
            zero,
            h,
            np.zeros((1, 3)),
            ["a"],
            parents,
            arm,
            rng,
        )
        assert not out["accepted"][0]
    args.output.write_text(json.dumps(metrics, indent=2) + "\n")


if __name__ == "__main__":
    main()

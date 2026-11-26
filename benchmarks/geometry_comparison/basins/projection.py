"""Fit the deposited coordinates on exactly the prepared PDBQT torsion tree."""

import numpy as np
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation


def project(template, heavy_serials, reference):
    lines = template.splitlines()
    serial_index, positions, branch_atoms, ends, stack = {}, [], [], [], []
    for line in lines:
        fields = line.split()
        if line.startswith("BRANCH "):
            ends.append(tuple(map(int, fields[1:])))
            branch_atoms.append([])
            stack.append(len(ends) - 1)
        elif line.startswith("ENDBRANCH "):
            if ends[stack.pop()] != tuple(map(int, fields[1:])):
                raise ValueError("Unbalanced branch")
        elif line.startswith(("ATOM  ", "HETATM")):
            serial = int(line[6:11])
            if serial in serial_index:
                raise ValueError("Duplicate atom serial")
            index = len(positions)
            serial_index[serial] = index
            positions.append([float(line[i : i + 8]) for i in (30, 38, 46)])
            for b in stack:
                branch_atoms[b].append(index)
    if stack:
        raise ValueError("Unclosed torsion branch")
    original = np.asarray(positions)
    axes = [(serial_index[i], serial_index[j]) for i, j in ends]
    atoms = [np.asarray(x) for x in branch_atoms]
    heavy = [serial_index[i] for i in heavy_serials]
    center = reference.mean(axis=0)
    delta = reference - center

    def place(angles):
        xyz = original.copy()
        for (a, b), ids, angle in zip(axes, atoms, angles):
            axis = xyz[b] - xyz[a]
            rotation = Rotation.from_rotvec(angle * axis / np.linalg.norm(axis))
            xyz[ids] = rotation.apply(xyz[ids] - xyz[a]) + xyz[a]
        moving = xyz[heavy]
        origin = moving.mean(axis=0)
        u, _, vt = np.linalg.svd((moving - origin).T @ delta)
        u[:, -1] *= np.sign(np.linalg.det(u @ vt))
        return (xyz - origin) @ (u @ vt) + center

    def residual(angles):
        return (place(angles)[heavy] - reference).ravel()

    rng = np.random.default_rng(43103)
    best_xyz, best_ss = place(np.zeros(len(axes))), np.inf
    if axes:
        for i in range(6):
            start = (
                np.zeros(len(axes)) if i == 0 else rng.uniform(-np.pi, np.pi, len(axes))
            )
            fit = least_squares(residual, start, max_nfev=150)
            ss = fit.fun @ fit.fun
            if ss < best_ss:
                best_xyz, best_ss = place(fit.x), ss
    result = []
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            xyz = best_xyz[serial_index[int(line[6:11])]]
            line = line[:30] + "".join(f"{v:8.3f}" for v in xyz) + line[54:]
        result.append(line)
    return "\n".join(result) + "\n"

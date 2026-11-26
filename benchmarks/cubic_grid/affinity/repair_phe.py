"""Resolve the incomplete SYK Phe382 ring against observed protein atoms."""

import numpy as np
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation


def place_ring(lines, records):
    index = {
        l[12:16].strip(): i
        for i, l in enumerate(lines)
        if l.startswith("ATOM  ") and l[22:26].strip() == "382"
    }

    def xyz(line):
        return np.array([float(line[k : k + 8]) for k in (30, 38, 46)])

    ring = [index[n] for n in ("CG", "CD1", "CD2", "CE1", "CE2", "CZ")]
    pos = np.array([xyz(lines[i]) for i in ring])
    env = np.array(
        [
            xyz(l)
            for l in records
            if l.startswith("ATOM  ") and l[22:26].strip() != "382"
        ]
    )
    ca, cb = [xyz(lines[index[n]]) for n in ("CA", "CB")]

    def rotate(coords, origin, axis, angle):
        a = axis / np.linalg.norm(axis) * np.deg2rad(angle)
        return Rotation.from_rotvec(a).apply(coords - origin) + origin

    best = None
    for chi1 in range(-180, 180, 15):
        for chi2 in range(-90, 90, 15):
            coord = rotate(pos, pos[0], pos[0] - cb, chi2)
            coord = rotate(coord, cb, cb - ca, chi1)
            d = cdist(coord, env)
            # Equal weight on each pair closer than a carbon heavy-atom contact.
            cost = float(np.maximum(0, 3.5 - d).__pow__(2).sum())
            item = (cost, abs(chi1) + abs(chi2), chi1, chi2)
            if best is None or item < best[0]:
                best = (item, coord, float(d.min()))
    for i, p in zip(ring, best[1]):
        lines[i] = lines[i][:30] + "".join(f"{v:8.3f}" for v in p) + lines[i][54:]
    return {
        "chi1_delta_deg": best[0][2],
        "chi2_delta_deg": best[0][3],
        "closest_observed_heavy_A": best[2],
        "overlap_sum_A2": best[0][0],
    }


def place_distal(lines, resid, root_name, end_name, mobile_names):
    """Rotate reconstructed distal atoms around one observed bond to clear clashes."""

    def xyz(line):
        return np.array([float(line[k : k + 8]) for k in (30, 38, 46)])

    index = {
        l[12:16].strip(): i
        for i, l in enumerate(lines)
        if l.startswith("ATOM") and int(l[22:26]) == resid
    }
    mobile = [index[n] for n in mobile_names]
    points = np.array([xyz(lines[i]) for i in mobile])
    origin = xyz(lines[index[end_name]])
    axis = origin - xyz(lines[index[root_name]])
    env = [xyz(l) for l in lines if l.startswith("ATOM") and int(l[22:26]) != resid]
    best = None
    for angle in range(-180, 180, 15):
        a = np.deg2rad(angle) * axis / np.linalg.norm(axis)
        pos = Rotation.from_rotvec(a).apply(points - origin) + origin
        dist = cdist(pos, env)
        key = (float((np.maximum(0, 3.1 - dist) ** 2).sum()), abs(angle), angle)
        if best is None or key < best[0]:
            best = key, pos, float(dist.min())
    for i, p in zip(mobile, best[1]):
        lines[i] = lines[i][:30] + "".join(f"{v:8.3f}" for v in p) + lines[i][54:]
    return {"rotation_deg": best[0][2], "min_external_A": best[2]}

"""Check water partition limits and independently recompute actual fixed-pose scores."""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from evaluate import iterate, read
from water_occupancy import RT, score_components

METHODS = ("exclusion", "contact", "joint", "replace_any", "replace_polar")


def taper(x, a, b):
    t = max(0.0, min(1.0, (x - a) / (b - a)))
    return 3 * t * t - 2 * t * t * t


def scalar(xyz, roles, sites):
    total = np.zeros(5)
    for site in sites:
        p = site["p"]
        loc = site["xyz"]
        near = [
            (math.dist(point, loc), point, role)
            for point, role in zip(xyz, roles, strict=True)
            if math.dist(point, loc) < 3.5
        ]
        if not near:
            continue
        availability = taper(min(p[0] for p in near), 2.0, 2.6)
        strength = 0.0
        any_score = 0.0
        polar_score = 0.0
        for r, point, role in near:
            geom_scores = []
            compatible = []
            for q, qrole in zip(site["partners"], site["partner_roles"], strict=True):
                a, b = np.array(point) - loc, np.array(q) - loc
                cos = np.dot(a, b) / max(1e-14, np.linalg.norm(a) * np.linalg.norm(b))
                theta = math.degrees(math.acos(max(-1.0, min(1.0, cos))))
                angular = taper(theta, 65.0, 85.0) * (1 - taper(theta, 140.0, 160.0))
                if role:
                    strength = max(
                        strength,
                        taper(r, 2.4, 2.7) * (1 - taper(r, 3.0, 3.5)) * angular,
                    )
                v, w = np.array(point) - q, np.array(loc) - q
                sep = np.linalg.norm(v)
                cos = np.dot(v, w) / max(1e-14, sep * np.linalg.norm(w))
                angle = math.degrees(math.acos(max(-1.0, min(1.0, cos))))
                h = (
                    taper(sep, 2.4, 2.7)
                    * (1 - taper(sep, 3.0, 3.5))
                    * (1 - taper(angle, 20.0, 45.0))
                )
                geom_scores.append(h)
                compatible.append(
                    h if ((role & 1 and qrole & 2) or (role & 2 and qrole & 1)) else 0.0
                )
            if r < 2.6:
                any_score = max(any_score, sum(geom_scores) / len(geom_scores))
                polar_score = max(polar_score, sum(compatible) / len(compatible))
        gain = math.exp(0.587439 * strength / RT)
        joint = availability * gain
        for k, z in enumerate(
            (
                availability,
                gain,
                joint,
                joint + (1 - availability) * any_score,
                joint + (1 - availability) * polar_score,
            )
        ):
            total[k] += -RT * math.log(1 - p + p * z)
    return total


def checks():
    sites = [
        {
            "xyz": [0, 0, 0],
            "p": 0.8,
            "partners": [[0.0, 2.85, 0.0]],
            "partner_roles": [2],
        }
    ]
    for x, role in (
        ([0, 0, 0], 0),
        ([0, 0, 0], 1),
        ([0, 0, 0], 2),
        ([1, 2, 0], 1),
        ([2.8, -0.5, 0], 1),
        ([1.15, 2.68, 0], 1),
        ([9, 9, 9], 0),
    ):
        xyz = np.array([x])
        roles = np.array([role])
        scores = score_components(sites, xyz, roles)
        np.testing.assert_allclose(
            [scores[k] for k in METHODS], scalar(xyz, roles, sites), atol=1.0e-12
        )
    at = np.zeros((1, 3))
    excluded = score_components(sites, at, np.array([0]))
    assert math.isclose(excluded["joint"], -RT * math.log(0.2))
    assert excluded["replace_any"] < 1e-12 and math.isclose(
        excluded["replace_polar"], excluded["joint"]
    )
    assert abs(score_components(sites, at, np.array([1]))["replace_polar"]) < 1e-12
    assert all(
        v == 0
        for v in score_components(
            sites, np.array([[10.0, 10.0, 10.0]]), np.array([1])
        ).values()
    )
    print("Partition, saturated clash, chemical substitution and null limits passed.")


def atom_roles(lines):
    atoms = [
        (l[77:].strip(), [float(l[i : i + 8]) for i in (30, 38, 46)])
        for l in lines
        if l.startswith(("ATOM  ", "HETATM"))
    ]
    hyd = [xyz for t, xyz in atoms if t == "HD"]
    heavy = [
        (t, xyz)
        for t, xyz in atoms
        if t not in ("HD", "H", "HS", "W", "G0", "G1", "G2", "G3")
    ]
    ranks = [
        int(t in ("OA", "NA"))
        + 2
        * int(
            t in ("O", "N", "NA", "OA") and any(math.dist(q, xyz) < 1.25 for q in hyd)
        )
        for t, xyz in heavy
    ]
    return np.array([q for _, q in heavy]), np.array(ranks)


def blocks(path):
    models = {}
    key = None
    for l in path.read_text().splitlines():
        if l.startswith("MODEL "):
            key = int(l.split()[1])
            models[key] = []
        elif l.startswith(("ATOM  ", "HETATM")):
            models[key].append(l)
    return models


def real(args):
    from functools import lru_cache

    load = lru_cache(maxsize=32)(blocks)
    saved = {
        (r["panel"], r["rec"], r["lig"], r["seed"], r["arm"], r["candidate"]): r
        for r in read(args.out / "poses.csv")
    }
    worst = 0.0
    n = 0
    for name, rec, lig, seed, arm, _, models in iterate(args):
        sites = json.loads((args.out / f"{rec}_sites.json").read_text())["sites"]
        for mid, path, index, native, _, _ in models:
            values = scalar(*atom_roles(load(path)[index]), sites) + native
            expected = saved[name, rec, lig, seed, arm, mid]
            worst = max(
                worst,
                max(abs(values[i] - float(expected[m])) for i, m in enumerate(METHODS)),
            )
            n += 1
    for row in read(args.holdout / "candidates.csv"):
        target, rec, seed, arm = row["target"], row["receptor"], row["seed"], row["arm"]
        width = int(row["box"])
        work = args.large if width == 36 else args.small
        pdb = (
            "1AZM"
            if row["ligand"] == "AZM" and target == "CAI"
            else "1BZM"
            if target == "CAI"
            else "1JD0"
        )
        sites = json.loads(
            (args.holdout / f"{target}_{rec}_{width}_sites.json").read_text()
        )["sites"]
        xyz, types = atom_roles(
            load(work / "runs" / f"{rec}_{pdb}_{seed}_{arm}.pdbqt")[int(row["model"])]
        )
        values = scalar(xyz, types, sites) + float(row["native"])
        worst = max(
            worst, max(abs(values[i] - float(row[m])) for i, m in enumerate(METHODS))
        )
        n += 1
    assert worst < 1e-10
    print(
        f"Independent scalar partition and atom-role parsing reproduced all five scores on {n} coordinates, max absolute difference {worst:.3g}."
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for k in (
        "out",
        "holdout",
        "brd_pools",
        "ca_pools",
        "state_pools",
        "large",
        "small",
    ):
        parser.add_argument("--" + k.replace("_", "-"), type=Path)
    a = parser.parse_args()
    checks()
    if a.out:
        real(a)

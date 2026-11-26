"""Rank an existing native Late pose pool with its independent search landmarks."""

import argparse
import json
import math
from pathlib import Path

from select_zinc_pose import parse_atoms, read_models, read_score_header


def read_sites(text):
    sites = {}
    for number, line in enumerate(text.splitlines(), 1):
        fields = line.split()
        if not fields:
            continue
        if len(fields) != 4 and not (len(fields) == 6 and fields[4] == "shell"):
            raise ValueError(f"Invalid landmark at line {number}")
        if not (
            1 <= len(fields[0]) <= 5 and fields[0].isascii() and fields[0].isdigit()
        ):
            raise ValueError(f"Invalid atom serial at line {number}")
        serial = int(fields[0])
        point = tuple(float(x) for x in fields[1:4])
        radius = float(fields[5]) if len(fields) == 6 else 0.0
        if (
            serial <= 0
            or serial in sites
            or not all(math.isfinite(x) for x in (*point, radius))
        ):
            raise ValueError(f"Invalid or repeated landmark at line {number}")
        if len(fields) == 6 and radius <= 0:
            raise ValueError(f"Expected a positive shell radius at line {number}")
        sites[serial] = point, radius
    if not sites:
        raise ValueError("Empty landmark file")
    return sites


def anchor_penalty(sites, coords):
    energy = 0.0
    for serial, (point, radius) in sites.items():
        d = math.dist(coords[serial], point)
        if radius:
            d = math.hypot(d, 1e-4) - math.hypot(radius, 1e-4)
        energy += 4 * d * d / (math.hypot(1, d) + 1)
    return energy / len(sites)


def ranked_pool(text, sites):
    rows = []
    for model, lines in read_models(text.splitlines()).items():
        # Only native Late headers specify the expected component ordering.
        if not lines or not lines[0].startswith("REMARK VINA RESULT:"):
            raise ValueError("Use a native scored PDBQT pose pool")
        saved = read_score_header(model, lines)["native_score"]
        prefix = "REMARK INTER + INTRA:"
        if len(lines) < 2 or not lines[1].startswith(prefix):
            raise ValueError("Missing native INTER + INTRA in opening header")
        raw = float(lines[1][len(prefix) :])
        if not math.isfinite(raw):
            raise ValueError("Nonfinite raw native energy")
        atoms = parse_atoms(lines)
        if any(
            i not in atoms
            or atoms[i].atom_type in {"H", "HD", "HS", "W", "G0", "G1", "G2", "G3"}
            for i in sites
        ):
            raise ValueError("Landmark must identify a real heavy atom in every model")
        bias = anchor_penalty(sites, {i: atoms[i].xyz for i in sites})
        rows.append(
            {
                "model": model,
                "native_affinity": saved,
                "native_raw": raw,
                "landmark_penalty": bias,
                "guided_pose_score": raw + bias,
            }
        )
    if not rows:
        raise ValueError("Empty pose pool")

    def order(key):
        return [r["model"] for r in sorted(rows, key=lambda r: (r[key], r["model"]))]

    return {
        "poses": rows,
        "native_order": [r["model"] for r in rows],
        "guided_order": order("guided_pose_score"),
        "geometry_order": order("landmark_penalty"),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pool", type=Path)
    parser.add_argument("sites", type=Path)
    args = parser.parse_args()
    print(
        json.dumps(
            ranked_pool(args.pool.read_text(), read_sites(args.sites.read_text())),
            indent=2,
            allow_nan=False,
        )
    )

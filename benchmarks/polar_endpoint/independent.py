"""Independent float64 evaluation of the HCT integrals, OBC2 Born radii and GB energy."""

import argparse
from pathlib import Path

import numpy as np
from polar_endpoint import atoms, read_models
from run import read, write
from scipy.spatial.distance import cdist


def gb_energy(molecules):
    xyz = np.concatenate([m.xyz for m in molecules])
    radii = np.concatenate([m.radii for m in molecules])
    q = np.concatenate([m.charges for m in molecules])
    offset = radii - 0.09
    size = offset * np.concatenate([m.scales for m in molecules])
    born = np.empty(len(q))
    for lo in range(0, len(q), 96):
        sl = slice(lo, lo + 96)
        d = cdist(xyz[sl], xyz)
        diagonal = d == 0
        d[diagonal] = 1
        upper = d + size
        lower = np.maximum(offset[sl, None], np.abs(d - size))
        integrals = 0.5 * (
            1 / lower
            - 1 / upper
            + 0.25 * (d - size**2 / d) * (1 / upper**2 - 1 / lower**2)
            + 0.5 * np.log(lower / upper) / d
        )
        integrals[diagonal | (upper <= offset[sl, None])] = 0.0
        psi = offset[sl] * np.sum(integrals, axis=1)
        born[sl] = 1 / (
            1 / offset[sl] - np.tanh(psi - 0.8 * psi**2 + 4.85 * psi**3) / radii[sl]
        )
    total = 0.0
    for lo in range(0, len(q), 96):
        sl = slice(lo, lo + 96)
        d2 = cdist(xyz[sl], xyz, metric="sqeuclidean")
        bij = born[sl, None] * born
        total += np.sum(q[sl, None] * q / np.sqrt(d2 + bij * np.exp(-d2 / (4 * bij))))
    # Independent CODATA conversion, kcal mol^-1 Å e^-2.
    return -0.5 * 332.06371 * (1 / 4 - 1 / 78.5) * total


def verify(target_results, prepared, runs):
    rows = read(target_results)
    sources = [rows[i] for i in (0, len(rows) // 2, -1)]
    out = []
    for row in sources:
        rec_path = prepared / (row["receptor"] + "_receptor.pdbqt")
        rec = atoms(rec_path.read_text().splitlines(), receptor=True)
        pose = read_models((runs / (row["stem"] + ".pdbqt")).read_text().splitlines())[
            int(row["model"])
        ]
        lig = atoms(pose)
        delta = gb_energy([rec, lig]) - gb_energy([rec]) - gb_energy([lig])
        d = cdist(lig.xyz, rec.xyz)
        coul = 332.06371 / 4 * np.sum(lig.charges[:, None] * rec.charges / d)
        delta_err = abs(delta - float(row["delta_gb"]))
        direct_err = abs(coul - float(row["coulomb"]))
        out.append(
            {
                "target": row["target"],
                "stem": row["stem"],
                "model": row["model"],
                "delta_gb": delta,
                "abs_delta_gb_error": delta_err,
                "coulomb": coul,
                "abs_coulomb_error": direct_err,
            }
        )
        if delta_err > 0.015 or direct_err > 0.0002:
            raise AssertionError(out[-1])
    return out


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("output", type=Path)
    p.add_argument(
        "inputs",
        nargs="+",
        help="Triples of score CSV, prepared directory, pose directory",
    )
    a = p.parse_args()
    if len(a.inputs) % 3:
        p.error("Expected triples")
    results = []
    for t in range(0, len(a.inputs), 3):
        results.extend(verify(*map(Path, a.inputs[t : t + 3])))
    write(a.output, results)

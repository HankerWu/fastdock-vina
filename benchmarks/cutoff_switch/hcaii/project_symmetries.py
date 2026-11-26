"""Fit the prepared, rigid-fragment ligand to a known crystal under every graph symmetry."""

import argparse
import json
import subprocess
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from ammonium_states import mirror_state
from coordination_probe import ligand_mapping, symmetry_mappings, write_csv
from crystal_start import atom_coordinates, project
from score_physical_pool import mapped_ligand_bonds


def fit_one(work, ligand, reflected):
    prepared = work / "prepared"
    base = (prepared / f"{ligand}_ligand.pdbqt").read_text()
    if reflected:
        base, _ = mirror_state(prepared, ligand)
    metadata = json.loads((prepared / "metadata.json").read_text())
    reference = np.array(metadata[ligand]["references"][ligand]["coords"])
    edges = mapped_ligand_bonds(base)
    _, serials, _ = ligand_mapping(prepared, ligand)
    best = None
    diagnostics = []
    for i, permutation in enumerate(symmetry_mappings(prepared, ligand)):
        target = reference[permutation]
        pose = project(base, list(serials), target, edges, analytic=True)
        coords = atom_coordinates(pose)
        rms = float(
            np.sqrt(
                np.mean(
                    np.sum(
                        (np.array([coords[s] for s in serials]) - target) ** 2, axis=1
                    )
                )
            )
        )
        diagnostics.append(
            {
                "ligand": ligand,
                "index": i,
                "fitted_mapping_rmsd": rms,
                "permutation": " ".join(str(x) for x in permutation),
            }
        )
        if best is None or rms < best[0]:
            best = rms, pose, i
    return ligand, best[1], best[2], diagnostics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--reflected",
        action="store_true",
        help="fit the other admissible ammonium state",
    )
    parser.add_argument("--jobs", type=int, default=3)
    args = parser.parse_args()
    work, output, binary = (
        args.work.resolve(),
        args.output.resolve(),
        args.binary.resolve(),
    )
    output.mkdir(parents=True, exist_ok=False)
    metadata = json.loads((work / "prepared" / "metadata.json").read_text())
    diagnostics = []
    with ProcessPoolExecutor(args.jobs) as pool:
        futures = [
            pool.submit(fit_one, work, ligand, args.reflected) for ligand in metadata
        ]
        for future in futures:
            ligand, pose, best_index, rows = future.result()
            diagnostics.extend(
                dict(**r, selected=int(r["index"] == best_index)) for r in rows
            )
            projected = output / f"{ligand}_start.pdbqt"
            projected.write_text(pose)
            subprocess.run(
                [
                    str(binary),
                    "--receptor",
                    str(work / "prepared" / f"{ligand}_receptor.pdbqt"),
                    "--config",
                    str(work / "prepared" / f"{ligand}_box.txt"),
                    "--ligand",
                    str(projected),
                    "--ligand_bonds",
                    str(work / "crystal_start" / f"{ligand}.bonds"),
                    "--local_only",
                    "--cpu",
                    "1",
                    "--seed",
                    "42",
                    "--verbosity",
                    "0",
                    "--out",
                    str(output / f"{ligand}_native_refined.pdbqt"),
                ],
                check=True,
                stdout=subprocess.PIPE,
                text=True,
            )
    write_csv(output / "fitted_mappings.csv", diagnostics)


if __name__ == "__main__":
    main()

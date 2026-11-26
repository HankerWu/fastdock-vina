"""Challenge the graph-only donor tie break and audit crystal-to-crystal transfer geometry on BRD4."""

import argparse
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from rdkit.Chem import rdFMCS

from prepare import RECEPTORS, all_correspondences, format_assignments, ligand_graph
from dock_matrix import run_one
from score_matrix import save
from make_scaffold_anchors import mcs_parameters
from select_zinc_pose import parse_atoms, read_models


def positions(entry, receptor, serials):
    coords = dict(
        zip(
            entry["heavy_serials"],
            entry["alignments"][receptor]["reference_by_match"][0],
            strict=True,
        )
    )
    return {i: coords[s] for i, s in serials.items()}


def check(prepared, baseline, output, binary, jobs):
    output.mkdir(exist_ok=True, parents=True)
    metadata = json.loads((prepared / "manifest.json").read_text())
    assigned = json.loads((baseline / "graphs.json").read_text())
    graphs = {k: ligand_graph(prepared, k) for k in metadata}
    references = []
    jobs_list = []
    for query, (q, qserial) in graphs.items():
        for donor, (d, dserial) in graphs.items():
            if query == donor:
                continue
            mcs = rdFMCS.FindMCS([q, d], mcs_parameters())
            if mcs.canceled:
                raise ValueError("Incomplete donor comparison")
            if mcs.numAtoms != assigned[query]["atoms"]:
                continue
            mappings = all_correspondences(q, d)
            for rec in RECEPTORS["brd4"]:
                donor_xyz = positions(metadata[donor], rec, dserial)
                query_xyz = positions(metadata[query], rec, qserial)
                distances = [
                    np.sqrt(
                        np.mean(
                            [
                                np.sum(
                                    np.square(np.subtract(query_xyz[a], donor_xyz[b]))
                                )
                                for a, b in m
                            ]
                        )
                    )
                    for m in mappings
                ]
                habitual = donor == assigned[query]["donor"]
                references.append(
                    dict(
                        receptor=rec,
                        query=query,
                        donor=donor,
                        usual=int(habitual),
                        atoms=mcs.numAtoms,
                        mappings=len(mappings),
                        closest_crystal_scaffold_rmsd=min(distances),
                    )
                )
                if habitual:
                    continue
                (output / f"{rec}_{query}.{donor}").write_text(
                    format_assignments(mappings, qserial, donor_xyz)
                )
                jobs_list.extend(
                    (rec, query, seed, donor)
                    for seed in (314159, 271828, 161803, 141421)
                )
    with ThreadPoolExecutor(jobs) as pool:
        list(
            pool.map(
                lambda cell: run_one(binary.resolve(), prepared, output, cell),
                jobs_list,
            )
        )
    save(output / "donor_geometry.csv", references)
    cells = []
    for rec, query, seed, donor in jobs_list:
        meta = metadata[query]
        for arm, path in (
            (donor, output / f"{rec}_{query}_{seed}_{donor}.pdbqt"),
            (assigned[query]["donor"], baseline / f"{rec}_{query}_{seed}_multi.pdbqt"),
        ):
            rmsds = []
            for lines in read_models(path.read_text().splitlines()).values():
                atoms = parse_atoms(lines)
                coords = [atoms[s].xyz for s in meta["heavy_serials"]]
                refs = meta["alignments"][rec]["reference_by_match"]
                rmsds.append(
                    np.sqrt(
                        np.square(np.subtract(coords, refs)).sum(axis=2).mean(axis=1)
                    ).min()
                )
            cells.append(
                dict(
                    receptor=rec,
                    query=query,
                    seed=seed,
                    donor=arm,
                    usual=int(arm == assigned[query]["donor"]),
                    top=rmsds[0],
                    best=min(rmsds),
                )
            )
    save(output / "donor_cells.csv", cells)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("prepared", "baseline", "output", "binary"):
        p.add_argument(key, type=Path)
    p.add_argument("--jobs", default=4, type=int)
    a = p.parse_args()
    check(a.prepared, a.baseline, a.output, a.binary, a.jobs)

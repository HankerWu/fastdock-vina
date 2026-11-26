"""Cross the two best non-query crystal scaffold donors at fixed search budget on hCA II.

Use the same two ligand states, trajectory slots, seeds, native search potential and
common-score evaluation as guidance.py. The only search input changed is the donor.
"""

import argparse
import itertools
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from evaluate import ligand_mapping, write
from guidance import ROOT, SEEDS, endpoints, pools, score_case, sources, summarize
from run import dock, mapped_ligand_bonds

sys.path.insert(0, str(ROOT / "benchmarks/mapping_ensemble"))
from prepare import (
    all_correspondences,
    choose_reference,
    declared_ca_graph,
    format_assignments,
)


def donor_anchors(cases, work):
    prepared = cases[0][3]["orig"].parent
    data = json.loads((prepared / "metadata.json").read_text())
    mapped = {lig: ligand_mapping(prepared, lig)[:2] for lig in data}
    # Strip the incidental stereochemistry of the 3D input and validate against
    # its declared PDBQT chemical graph. MCS selection consumes no coordinates.
    graphs = {
        lig: declared_ca_graph(prepared, lig, mol) for lig, (mol, _) in mapped.items()
    }
    chosen = {}
    for lig in sorted(graphs):
        remaining = graphs.copy()
        for rank in (1, 2):
            donor, _, _ = choose_reference(lig, remaining)
            del remaining[donor]
            maps = all_correspondences(graphs[lig], graphs[donor])
            chosen[lig, rank] = donor, maps

    records, paths = [], {}
    for _, receptor, ligand, _, _, _, shipped in cases:
        for rank in (1, 2):
            donor, maps = chosen[ligand, rank]
            # Only the OTHER ligand's experimental coordinates enter the prior.
            xyz = data[donor]["references"][receptor]["coords"]
            anchor = format_assignments(maps, mapped[ligand][1], xyz)
            if rank == 1 and anchor.encode() != shipped.read_bytes():
                raise ValueError(f"Known donor anchor mismatch for {receptor}/{ligand}")
            path = work / "inputs" / f"{receptor}_{ligand}_from_{donor}.anchor"
            path.write_text(anchor)
            paths[receptor, ligand, rank] = path
            records.append(
                {
                    "receptor": receptor,
                    "ligand": ligand,
                    "rank": rank,
                    "donor": donor,
                    "matched_atoms": len(maps[0]),
                    "assignments": len(maps),
                    "anchor": path.name,
                }
            )
    write(work / "donors.csv", records)
    return paths


def run(binary, base, work, workers, seeds):
    cases = [case for case in sources(base) if case[0] == "ca"]
    for name in ("inputs", "baseline", "guided"):
        (work / name).mkdir(parents=True, exist_ok=True)
    anchors = donor_anchors(cases, work)
    tasks = []
    for case in cases:
        panel, rec, lig, states, *_ = case
        for state, source in states.items():
            bonds = mapped_ligand_bonds(source.read_text())
            (work / "inputs" / f"{panel}_{lig}_{state}.bonds").write_text(
                "".join(f"{i} {j}\n" for i, j in bonds)
            )
        for (_, _, _, arm), path, (receptor, ligand, box, _, bonds, seed), _ in pools(
            case, work, seeds=seeds
        ):
            # Use guidance.py's file layout to share its scoring path. "guided"
            # denotes donor rank 2 here. Neither branch uses amide search guidance.
            anchor = anchors[rec, lig, 1 if arm == "baseline" else 2]
            tasks.append((path, (receptor, ligand, box, anchor, bonds, seed)))
    with ThreadPoolExecutor(workers) as executor:
        list(executor.map(lambda task: dock(binary, task), tasks))
        rows = list(
            itertools.chain.from_iterable(
                executor.map(
                    lambda case: score_case(binary, work, case, {}, seeds=seeds), cases
                )
            )
        )
    cells = summarize(rows)
    labels = {
        "baseline": "donor1",
        "guided": "donor2",
        "guided_slot0": "donor2_slot0",
        "guided_slot1": "donor2_slot1",
    }
    for row in rows:
        row["search"] = labels[row["search"]]
    for cell in cells:
        cell["arm"] = labels[cell["arm"]]
    write(work / "modes.csv", rows)
    write(work / "cells.csv", cells)
    write(work / "summary.csv", endpoints(cells))
    print(f"Scored {len(rows)} models from {len(tasks)} docking searches")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("binary", "base", "work"):
        p.add_argument(name, type=Path)
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument(
        "--seeds",
        nargs="+",
        type=int,
        default=SEEDS,
        help="paired nonzero base seeds (second slot adds 1009)",
    )
    args = p.parse_args()
    actual = [s + offset for s in args.seeds for offset in (0, 1009)]
    if min(actual) < 1 or max(actual) > 2**31 - 1 or len(set(actual)) != len(actual):
        p.error(
            "trajectory seeds, including the slot offsets, must be distinct positive 32-bit integers"
        )
    run(
        args.binary.resolve(),
        args.base.resolve(),
        args.work.resolve(),
        args.jobs,
        args.seeds,
    )

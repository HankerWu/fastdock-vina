"""Leave-one-ligand-out scaffold-guided search on the existing hCA II receptors."""

import argparse
import itertools
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from coordination_probe import ligand_mapping
from rdkit import Chem
from run_hcaii import LIGANDS, RECEPTORS, SEEDS, dock, evaluate
from transfer_mapping import choose_reference


def input_sites(prepared):
    graphs = {
        lig: Chem.SDMolSupplier(str(prepared / f"{lig}_graph.sdf"))[0]
        for lig in LIGANDS
    }
    data = json.loads((prepared / "metadata.json").read_text())
    sites, mapping = {}, {}
    for lig in LIGANDS:
        donor, query_idx, donor_idx = choose_reference(lig, graphs)
        _, serials, _ = ligand_mapping(prepared, lig)
        mapping[lig] = {
            "donor": donor,
            "query_atoms": query_idx,
            "donor_atoms": donor_idx,
            "heavy_atoms": graphs[lig].GetNumAtoms(),
            "matched": len(query_idx),
        }
        for rec in RECEPTORS:
            # Only the OTHER ligand's aligned experimental coordinates supply the prior.
            coords = data[donor]["references"][rec]["coords"]
            lines = [
                f"{serials[q]} " + " ".join(f"{c:.8f}" for c in coords[d])
                for q, d in zip(query_idx, donor_idx, strict=True)
            ]
            sites[rec, lig] = "\n".join(lines) + "\n"
    return sites, mapping


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("prepared", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--jobs", type=int, default=2)
    args = parser.parse_args()
    sites, mapping = input_sites(args.prepared)
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "mapping.json").write_text(json.dumps(mapping, indent=2) + "\n")
    cells = list(itertools.product(RECEPTORS, LIGANDS, SEEDS, ("native", "scaffold")))
    with ThreadPoolExecutor(args.jobs) as pool:
        list(
            pool.map(
                lambda c: dock(
                    args.binary.resolve(), args.prepared, args.output, c, sites[c[:2]]
                ),
                cells,
            )
        )
    evaluate(args.prepared, args.output, cells)


if __name__ == "__main__":
    main()

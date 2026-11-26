"""Create all distinct graph-supported correspondences for an independent crystal donor."""

import argparse
import json
import sys
from pathlib import Path

from rdkit import Chem

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "benchmarks/anchor_search"))
from transfer_mapping import choose_reference
from make_scaffold_anchors import all_correspondences, format_assignments, query_graph
from transfer_holdout import RECEPTORS, ligand_graph
from transfer_hcaii import ligand_mapping

RECEPTORS = dict(RECEPTORS, hcaii=("6WQ4", "6WQ8"))


def declared_ca_graph(prepared, ligand, geometry):
    # Keep the SDF atom ordering needed for deposited reference coordinates.
    # A 3D SDF reader assigns the incidental hand of the starting protonated
    # amine. The input SMILES for this specific panel leaves its only center open.
    declared, _ = query_graph((prepared / f"{ligand}_ligand.pdbqt").read_text())
    potential = Chem.FindPotentialStereo(declared)
    if (
        len(potential) != 1
        or potential[0].specified != Chem.StereoSpecified.Unspecified
        or declared.GetAtomWithIdx(potential[0].centeredOn).GetSymbol() != "N"
    ):
        raise ValueError(f"Expected one unspecified nitrogen in {ligand}")
    graph = Chem.Mol(geometry)
    Chem.RemoveStereochemistry(graph)
    if Chem.MolToSmiles(graph) != Chem.MolToSmiles(declared):
        raise ValueError(
            f"Source graph and declared input chemistry differ for {ligand}"
        )
    return graph


def prepare(panel, prepared, out, declared_nitrogen=False):
    if declared_nitrogen and panel != "hcaii":
        raise ValueError("The alternate nitrogen treatment applies to hCA II")
    out.mkdir(parents=True, exist_ok=True)
    file = "metadata.json" if panel == "hcaii" else "manifest.json"
    data = json.loads((prepared / file).read_text())
    if panel == "hcaii":
        mapped = {l: ligand_mapping(prepared, l)[:2] for l in data}
    else:
        mapped = {l: ligand_graph(prepared, l) for l in data}
    graphs = {l: value[0] for l, value in mapped.items()}
    declared = (
        {lig: declared_ca_graph(prepared, lig, mol) for lig, mol in graphs.items()}
        if declared_nitrogen
        else graphs
    )
    inventory = {}
    for lig in sorted(graphs):
        donor, qa, da = choose_reference(lig, graphs)
        maps = all_correspondences(declared[lig], declared[donor])
        historical = tuple(sorted(zip(qa, da, strict=True)))
        if not declared_nitrogen and historical not in maps:
            raise ValueError(f"Historical mapping missing from ensemble: {lig}/{donor}")
        inventory[lig] = dict(
            donor=donor,
            atoms=len(maps[0]),
            historical_subset=historical in maps,
            nitrogen="declared"
            if declared_nitrogen
            else "input_geometry"
            if panel == "hcaii"
            else "input_smiles",
            mappings=len(maps),
            single_indices=list(zip(qa, da, strict=True)),
            maps=maps,
        )
        for rec in RECEPTORS[panel]:
            if panel == "hcaii":
                coords = data[donor]["references"][rec]["coords"]
            else:
                ref = data[donor]["alignments"][rec]["reference_by_match"][0]
                positions = dict(zip(data[donor]["heavy_serials"], ref, strict=True))
                coords = {a: positions[s] for a, s in mapped[donor][1].items()}
            # Retain historical input order for exact single-assignment controls.
            old = [tuple(zip(qa, da, strict=True))]
            (out / f"{rec}_{lig}.single").write_text(
                format_assignments(old, mapped[lig][1], coords)
            )
            (out / f"{rec}_{lig}.multi").write_text(
                format_assignments(maps, mapped[lig][1], coords)
            )
    (out / "graphs.json").write_text(json.dumps(inventory, indent=1) + "\n")
    print(
        panel,
        [(k, v["mappings"], v["atoms"]) for k, v in inventory.items()],
        flush=True,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("panel", choices=RECEPTORS)
    p.add_argument("prepared", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument(
        "--declared-nitrogen",
        action="store_true",
        help="Use the unspecified nitrogen of the hCA II input SMILES for the multi arm; retain the historical donor and single arm.",
    )
    a = p.parse_args()
    prepare(a.panel, a.prepared, a.out, a.declared_nitrogen)

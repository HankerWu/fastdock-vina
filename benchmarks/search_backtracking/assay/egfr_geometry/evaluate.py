"""Fixed-frame redocking accuracy of recorded EGFR/erlotinib poses against deposited AQ4."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import gemmi
import numpy as np
from meeko import PDBQTMolecule, RDKitMolCreate
from rdkit import Chem
from scipy.spatial.distance import cdist

HERE = Path(__file__).resolve().parent


def pos(line):
    return [float(line[i : i + 8]) for i in (30, 38, 46)]


def graph_and_reference(work):
    ccd = gemmi.cif.read(str(HERE / "AQ4.cif")).sole_block()
    atoms = ccd.find(
        "_chem_comp_atom.", ["atom_id", "type_symbol", "charge", "pdbx_stereo_config"]
    )
    idx = {}
    mol = Chem.RWMol()
    for name, symbol, charge, stereo in atoms:
        if symbol != "H":
            assert stereo == "N"
            atom = Chem.Atom(symbol)
            atom.SetFormalCharge(int(charge))
            idx[name] = mol.AddAtom(atom)
    order = {
        "SING": Chem.BondType.SINGLE,
        "DOUB": Chem.BondType.DOUBLE,
        "TRIP": Chem.BondType.TRIPLE,
    }
    for a, b, o, s in ccd.find(
        "_chem_comp_bond.",
        ["atom_id_1", "atom_id_2", "value_order", "pdbx_stereo_config"],
    ):
        if a in idx and b in idx:
            assert s == "N"
            mol.AddBond(idx[a], idx[b], order[o.upper()])
    Chem.SanitizeMol(mol)
    deposited = (work / "raw/1M17.pdb").read_text().splitlines()
    atomrows = [
        s
        for s in deposited
        if s.startswith("HETATM") and s[17:20] == "AQ4" and s[21:26] == "A 999"
    ]
    assert all(s[16] == " " and float(s[54:60]) == 1 for s in atomrows)
    coords = {s[12:16].strip(): pos(s) for s in atomrows}
    assert (
        len(atomrows) == len(coords) == len(idx) == 29 and coords.keys() == idx.keys()
    )
    receptor = [
        s
        for s in (work / "prepared/1M17.pdbqt").read_text().splitlines()
        if s.startswith("ATOM") and s.split()[-1] not in {"H", "HD", "HS"}
    ]
    # Protein shares the original frame; no Kabsch fit is used on protein or ligand.
    ca = {
        s[21:27]: pos(s)
        for s in deposited
        if s.startswith("ATOM") and s[12:16].strip() == "CA" and s[16] in " A"
    }
    shared = [s for s in receptor if s[12:16].strip() == "CA" and s[21:27] in ca]
    largest = max(np.linalg.norm(np.array(pos(s)) - ca[s[21:27]]) for s in shared)
    assert largest < 0.005
    return (
        mol,
        np.array([coords[n] for n in idx]),
        list(idx),
        receptor,
        len(shared),
        float(largest),
    )


def assess(work, out):
    ccd, ref, names, receptor, nca, cadiff = graph_and_reference(work)
    rec = np.array([pos(s) for s in receptor])
    smiles = (work / "prepared/176870.pdbqt").read_text().splitlines()[0].split()[2]
    graph = Chem.MolFromSmiles(smiles)
    assert Chem.MolToSmiles(graph) == Chem.MolToSmiles(ccd)
    maps = ccd.GetSubstructMatches(
        graph, useChirality=True, uniquify=False, maxMatches=10000
    )
    assert 0 < len(maps) < 10000 and all(len(m) == 29 for m in maps)
    # Radius convention shared with benchmarks/pair_gradient/pose/paired.py.
    radii = {
        "C": 1.70,
        "A": 1.70,
        "N": 1.55,
        "NA": 1.55,
        "O": 1.52,
        "OA": 1.52,
        "S": 1.80,
        "SA": 1.80,
    }
    rec_rad = np.array([radii[s.split()[-1]] for s in receptor])

    def overlap_threshold(mol):
        return 0.75 * (
            np.array([radii[a.GetSymbol()] for a in mol.GetAtoms()])[:, None] + rec_rad
        )

    reference_gaps = cdist(ref, rec)
    i, j = np.unravel_index(reference_gaps.argmin(), reference_gaps.shape)
    summary = {
        "pdb": "1M17",
        "chain": "A",
        "residue": "AQ4 999",
        "cid": 176870,
        "ccd_sha256": hashlib.sha256((HERE / "AQ4.cif").read_bytes()).hexdigest(),
        "smiles": smiles,
        "heavy_atoms": len(names),
        "isomorphisms": len(maps),
        "reference_atom_names_in_query_order": [[names[i] for i in m] for m in maps],
        "matched_receptor_CA": nca,
        "max_CA_shift_A": cadiff,
        "reference_severe_pairs": int((reference_gaps < overlap_threshold(ccd)).sum()),
        "reference_closest_A": float(reference_gaps[i, j]),
        "reference_closest_ligand": names[i],
        "reference_closest_receptor": receptor[j][12:27].strip(),
        "reference_ligand_heavy_within_2A": int((reference_gaps.min(axis=1) < 2).sum()),
    }
    cells, modes = [], []
    for folder in ("experiments", "cost_check"):
        with (work / folder / "runs.csv").open() as f:
            runs = [
                r
                for r in csv.DictReader(f)
                if r["panel"] == "EGFR" and r["compound"] == "176870"
            ]
        for r in runs:
            key = r["key"]
            sections = (
                (work / folder / "raw" / (key + ".pdbqt"))
                .read_text()
                .split("ENDMDL")[:-1]
            )
            assert len(sections) == int(r["modes"])
            here = []
            for rank, section in enumerate(sections, 1):
                lines = section.strip().splitlines()
                indices = [
                    int(x)
                    for s in lines
                    if s.startswith("REMARK SMILES IDX ")
                    for x in s.split()[3:]
                ]
                mapping = dict(zip(indices[::2], indices[1::2], strict=True))
                all_xyz = {
                    int(s[6:11]): pos(s)
                    for s in lines
                    if s.startswith(("ATOM  ", "HETATM"))
                }
                points = np.array([all_xyz[mapping[k]] for k in range(1, 30)])
                # Validate the serial-to-graph interpretation with a separate Meeko decoder.
                (decoded,) = RDKitMolCreate.from_pdbqt_mol(
                    PDBQTMolecule(section.strip() + "\nENDMDL\n", skip_typing=True)
                )
                heavy = Chem.RemoveHs(decoded)
                assert Chem.MolToSmiles(heavy) == Chem.MolToSmiles(graph)
                assert (
                    np.max(np.abs(heavy.GetConformer().GetPositions() - points)) < 1e-4
                )
                err = min(
                    float(np.sqrt(np.square(points - ref[list(m)]).sum(axis=1).mean()))
                    for m in maps
                )
                distances = cdist(points, rec)
                min_dist = distances.min(axis=1)
                (score,) = [
                    float(s.split()[3])
                    for s in lines
                    if s.startswith("REMARK VINA RESULT:")
                ]
                here.append(
                    {
                        "key": key,
                        "rank": rank,
                        "score": score,
                        "rmsd_A": round(err, 6),
                        "closest_heavy_A": round(float(min_dist.min()), 6),
                        "severe_pairs": int(
                            (distances < overlap_threshold(graph)).sum()
                        ),
                        "ligand_heavy_within_2A": int((min_dist < 2).sum()),
                    }
                )
            modes.extend(here)
            dist = [m["rmsd_A"] for m in here]
            cells.append(
                {
                    "key": key,
                    "seed": r["seed"],
                    "cap": r["cap"],
                    "arm": r["arm"],
                    "modes": len(here),
                    "top1_A": dist[0],
                    "best5_A": min(dist[:5]),
                    "oracle_A": min(dist),
                    "oracle_rank": dist.index(min(dist)) + 1,
                    "top1_hit": int(dist[0] <= 2),
                    "top5_hit": int(min(dist[:5]) <= 2),
                    "oracle_hit": int(min(dist) <= 2),
                }
            )
    for name, rows in (("modes", modes), ("cells", cells)):
        with (out / (name + ".csv")).open("w") as stream:
            writer = csv.DictWriter(
                stream, fieldnames=list(rows[0]), lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
    summary.update(pools=len(cells), models=len(modes))
    (out / "reference.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("out", type=Path)
    args = p.parse_args()
    args.out.mkdir(exist_ok=True, parents=True)
    assess(args.work, args.out)

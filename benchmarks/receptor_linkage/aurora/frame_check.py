"""Recompute RMSDs with Gemmi and RDKit and test one alternative receptor frame."""

import argparse
import csv
import json
import tarfile
from pathlib import Path

import gemmi
import numpy as np
from geometry import METHODS, RECEPTORS, REFERENCES, SEEDS, summary_rows, write_csv
from rdkit import Chem
from rdkit.Chem import rdMolAlign


def ca(chain):
    return {
        (r.seqid.num, r.seqid.icode, r.name): a.pos
        for r in chain
        if r.het_flag == "A"
        for a in r
        if a.name == "CA" and a.altloc in ("\x00", "A")
    }


def untyped_graph(mol):
    mol = Chem.RWMol(mol)
    for b in mol.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE)
        b.SetIsAromatic(False)
    for a in mol.GetAtoms():
        a.SetIsAromatic(False)
        a.SetFormalCharge(0)
        a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        a.SetNumExplicitHs(0)
        a.SetNoImplicit(True)
    return mol.GetMol()


def read_query(run, drug):
    lines = (run / "prepared" / (drug + ".pdbqt")).read_text().splitlines()
    mol = Chem.MolFromSmiles(
        next(l.split()[2] for l in lines if l.startswith("REMARK SMILES "))
    )
    serial = [None] * mol.GetNumAtoms()
    for l in lines:
        if l.startswith("REMARK SMILES IDX "):
            v = list(map(int, l.split()[3:]))
            for a, b in zip(v[::2], v[1::2]):
                serial[a - 1] = b
    assert None not in serial
    mol = untyped_graph(mol)
    mol.AddConformer(Chem.Conformer(len(serial)))
    return mol, serial


def reference_graph(path, ligand, chain):
    lines = path.read_text().splitlines()
    block = [
        l
        for l in lines
        if l.startswith("HETATM")
        and l[21] == chain
        and l[17:20] == ligand
        and l[16] in " A"
        and l[76:78].strip() not in ("H", "D")
    ]
    ids = {int(l[6:11]) for l in block}
    for l in lines:
        if l.startswith("CONECT"):
            numbers = [
                int(l[k : k + 5])
                for k in range(6, len(l.rstrip()), 5)
                if l[k : k + 5].strip()
            ]
            if numbers[0] in ids:
                block.append("CONECT" + "".join(f"{i:5d}" for i in numbers if i in ids))
    mol = Chem.MolFromPDBBlock(
        "\n".join(block + ["END"]),
        sanitize=False,
        removeHs=False,
        proximityBonding=False,
    )
    return untyped_graph(mol)


def positions(path):
    result, model = [], {}
    for l in path.read_text().splitlines():
        if l.startswith(("ATOM  ", "HETATM")):
            model[int(l[6:11])] = tuple(float(l[k : k + 8]) for k in (30, 38, 46))
        elif l.startswith("ENDMDL"):
            result.append(model)
            model = {}
    return result


def check(run, refs, archive, results):
    pocket_errors = []
    alternates, alignments = [], []
    with (results / "geometry_modes.csv").open() as f:
        archived = {
            (
                r["reference"],
                r["chain"],
                r["receptor"],
                r["method"],
                int(r["seed"]),
                int(r["rank"]),
            ): float(r["rmsd"])
            for r in csv.DictReader(f)
        }
    with tarfile.open(archive) as f:
        receiver = {}
        for rec in RECEPTORS:
            s = gemmi.read_pdb_string(f.extractfile(f"raw/{rec}.pdb").read())
            s.merge_chain_parts()
            receiver[rec] = s[0]["A"].clone()
    for drug, entries in REFERENCES.items():
        mol, serials = read_query(run, drug)
        for pdb, ligand, chains in entries:
            path = refs / (pdb + ".pdb")
            src = gemmi.read_structure(str(path))
            for ch in chains:
                ref_ca = ca(src[0][ch])
                ref = reference_graph(path, ligand, ch)
                xyz = ref.GetConformer().GetPositions()
                assert (
                    mol.GetNumAtoms() == ref.GetNumAtoms()
                    and mol.GetNumBonds() == ref.GetNumBonds()
                )
                for rec, protein in receiver.items():
                    rca = ca(protein)
                    box = [
                        a.pos
                        for r in protein
                        if r.name == RECEPTORS[rec]
                        for a in r
                        if not a.element.is_hydrogen and a.altloc in ("\x00", "A")
                    ]
                    keys = sorted(rca.keys() & ref_ca.keys())
                    selected = [
                        k for k in keys if min(rca[k].dist(pos) for pos in box) <= 12
                    ]
                    for frame, anchors in [
                        ("pocket", selected),
                        ("all_observed_common_CA", keys),
                    ]:
                        fit = gemmi.superpose_positions(
                            [rca[k] for k in anchors], [ref_ca[k] for k in anchors]
                        )
                        coords = ref.GetConformer()
                        for i, v in enumerate(xyz):
                            coords.SetAtomPosition(
                                i, tuple(fit.transform.apply(gemmi.Position(*v)))
                            )
                        for k in anchors:
                            pos = fit.transform.apply(ref_ca[k])
                            gap = np.linalg.norm(
                                np.array(list(pos)) - np.array(list(rca[k]))
                            )
                            alignments.append(
                                {
                                    "reference": pdb,
                                    "chain": ch,
                                    "receptor": rec,
                                    "frame": frame,
                                    "anchors": len(anchors),
                                    "residue": k[0],
                                    "insertion": k[1].strip(),
                                    "resname": k[2],
                                    "activation_segment": int(274 <= k[0] <= 299),
                                    "CA_error_A": gap,
                                }
                            )
                        for seed in SEEDS:
                            for method in METHODS:
                                attempts = positions(
                                    run
                                    / "pools"
                                    / f"{method}_{rec}_{drug}_{seed}.pdbqt"
                                )
                                measured = []
                                for rank, mode in enumerate(attempts, 1):
                                    for i, n in enumerate(serials):
                                        mol.GetConformer().SetAtomPosition(i, mode[n])
                                    val = rdMolAlign.CalcRMS(
                                        mol,
                                        ref,
                                        symmetrizeConjugatedTerminalGroups=False,
                                    )
                                    measured.append(val)
                                    if frame == "pocket":
                                        pocket_errors.append(
                                            abs(
                                                archived[
                                                    pdb, ch, rec, method, seed, rank
                                                ]
                                                - val
                                            )
                                        )
                                if frame != "pocket":
                                    alternates.append(
                                        {
                                            "compound": drug,
                                            "reference": pdb,
                                            "chain": ch,
                                            "primary": pdb != "6CPG" and ch == "A",
                                            "receptor": rec,
                                            "seed": seed,
                                            "method": method,
                                            "count": len(measured),
                                            "top1": measured[0],
                                            "best5": min(measured[:5]),
                                            "oracle": min(measured),
                                            "regret": measured[0] - min(measured),
                                        }
                                    )
    assert len(archived) == len(pocket_errors) and max(pocket_errors) < 1e-8
    write_csv(results, "alternative_frame_cells.csv", alternates)
    write_csv(
        results, "alternative_frame_distributions.csv", list(summary_rows(alternates))
    )
    write_csv(results, "alignment_residuals.csv", alignments)
    (results / "geometry_independent.json").write_text(
        json.dumps(
            {
                "compared_mode_reference_rows": len(archived),
                "max_abs_difference_A": max(pocket_errors),
                "protein_fit": "gemmi.superpose_positions",
                "ligand_mapping": "PDB CONECT versus Meeko SMILES IDX",
                "metric": "RDKit CalcRMS; no ligand fit",
                "sensitivity": "all common observed ATOM C-alphas",
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("run", type=Path)
    p.add_argument("references", type=Path)
    p.add_argument("results", type=Path)
    p.add_argument("--kinase-archive", type=Path, required=True)
    a = p.parse_args()
    check(a.run, a.references, a.kinase_archive, a.results)

"""Independent crystal-instance audit of hCA II poses, graph mapping and central-N handedness."""

import argparse
import csv
import hashlib
import itertools
import urllib.request
from pathlib import Path

import numpy as np
from holdout_verify import backbone
from meeko import PDBQTMolecule, RDKitMolCreate
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from scipy.spatial.transform import Rotation


def audited_reference(work, code):
    path = work / "raw" / (code + "-instance.sdf")
    if not path.exists():
        url = f"https://models.rcsb.org/v1/{code.lower()}/ligand?auth_asym_id=A&auth_seq_id=303&encoding=sdf"
        with urllib.request.urlopen(url, timeout=60) as source:
            path.write_bytes(source.read())
    ref = Chem.SDMolSupplier(str(path))[0]
    Chem.RemoveStereochemistry(ref)
    # Retain deposited heavy coordinates. Reconstruct the declared tested microstate for graph mapping.
    for smarts, charge, hydrogens in (
        ("[N;H2;+0]S(=O)=O", -1, 1),
        ("[N;X3;+0;!$(NS);!$(NC=O);!$(N#*)]", 1, None),
        ("[O;H1;+0]C=O", -1, 0),
    ):
        for match in ref.GetSubstructMatches(Chem.MolFromSmarts(smarts)):
            atom = ref.GetAtomWithIdx(match[0])
            h = atom.GetTotalNumHs() + 1 if hydrogens is None else hydrogens
            atom.SetFormalCharge(charge)
            atom.SetNumExplicitHs(h)
            atom.SetNoImplicit(True)
            atom.UpdatePropertyCache(False)
    Chem.SanitizeMol(ref)
    return ref


def handedness(mol, xyz):
    (ammonium,) = mol.GetSubstructMatches(Chem.MolFromSmarts("[N+;H1;D3]"))
    idx = ammonium[0]
    ranks = Chem.CanonicalRankAtoms(mol, breakTies=False, includeChirality=False)
    neighbors = sorted(
        (a.GetIdx() for a in mol.GetAtomWithIdx(idx).GetNeighbors()),
        key=lambda i: ranks[i],
    )
    if len({ranks[i] for i in neighbors}) != 3:
        raise ValueError("Identical branches at the chosen nitrogen")
    return np.linalg.det(xyz[neighbors] - xyz[idx])


def verify(work, output):
    rows = list(csv.DictReader((output / "modes.csv").open()))
    pdbs = sorted({r["lig"] for r in rows})
    refs = {p: audited_reference(work, p) for p in pdbs}
    backbone_by_pdb = {p: backbone(work, p) for p in pdbs}
    key = lambda r: (r["rec"], r["lig"], r["seed"], r["arm"])
    errors, cells = [], []
    for (rec, lig, seed, arm), group in itertools.groupby(sorted(rows, key=key), key):
        entries = sorted(group, key=lambda r: int(r["model"]))
        crystal = Chem.Mol(refs[lig])
        ref_positions = crystal.GetConformer().GetPositions()
        src, dst = backbone_by_pdb[lig], backbone_by_pdb[rec]
        keys = sorted(src.keys() & dst.keys())
        a, b = (np.array([c[k] for k in keys]) for c in (src, dst))
        rot, _ = Rotation.align_vectors(b - b.mean(0), a - a.mean(0))
        crystal.GetConformer().SetPositions(
            rot.apply(ref_positions - a.mean(0)) + b.mean(0)
        )
        sign = handedness(crystal, ref_positions)
        name = output / f"{rec}_{lig}_{seed}_{arm}.pdbqt"
        probe = Chem.RemoveHs(
            RDKitMolCreate.from_pdbqt_mol(
                PDBQTMolecule.from_file(str(name), skip_typing=True)
            )[0]
        )
        Chem.RemoveStereochemistry(probe)
        if Chem.MolToSmiles(probe) != Chem.MolToSmiles(crystal):
            raise ValueError(
                f"Chemical graph differs {name}: {Chem.MolToSmiles(probe)} vs {Chem.MolToSmiles(crystal)}"
            )
        hits = []
        for conf, row in zip(probe.GetConformers(), entries, strict=True):
            rms = rdMolAlign.CalcRMS(
                probe, crystal, prbId=conf.GetId(), maxMatches=1000000
            )
            errors.append(abs(rms - float(row["rmsd"])))
            same = handedness(probe, conf.GetPositions()) * sign > 0
            if "sign_ok" in row and bool(int(row["sign_ok"])) != same:
                raise ValueError(f"Nitrogen handedness differs in {name} model {row['model']}")
            hits.append((rms <= 2, same))
        cells.append(
            {
                "rec": rec,
                "lig": lig,
                "seed": seed,
                "arm": arm,
                "sign_match": hits[0][1],
                "rmsd_hit": hits[0][0],
                "both_top": all(hits[0]),
                "both_pool": any(all(h) for h in hits),
            }
        )
    if max(errors) > 1e-5:
        raise ValueError(f"RMSD discrepancy {max(errors)}")
    with (output / "independent_sign.csv").open("w") as fp:
        writer = csv.DictWriter(fp, fieldnames=cells[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(cells)
    print(
        f"{len(errors)} raw poses checked; maximum absolute RMSD difference {max(errors):.8f} angstrom"
    )
    for p in pdbs:
        print(
            p,
            hashlib.sha256(
                (work / "raw" / (p + "-instance.sdf")).read_bytes()
            ).hexdigest(),
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    verify(args.work, args.output)

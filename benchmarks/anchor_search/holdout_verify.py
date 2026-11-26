"""Independently rederive all holdout RMSDs from RCSB instance SDFs and Meeko."""

import argparse
import csv
import hashlib
import itertools
import urllib.request
from pathlib import Path

import gemmi
import numpy as np
from meeko import PDBQTMolecule, RDKitMolCreate
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from scipy.spatial.transform import Rotation

# Chain A active-site residue identifiers in each deposited structure.
RESIDUES = {"1AZM": 262, "1BZM": 262, "1JD0": 1400}


def instance(work, pdb):
    path = work / "raw" / f"{pdb}-instance.sdf"
    source = (
        f"https://models.rcsb.org/v1/{pdb.lower()}/ligand?auth_asym_id=A"
        f"&auth_seq_id={RESIDUES[pdb]}&encoding=sdf"
    )
    if not path.exists():
        with urllib.request.urlopen(source, timeout=60) as response:
            path.write_bytes(response.read())
    mol = Chem.SDMolSupplier(str(path), removeHs=True)[0]
    Chem.AssignStereochemistryFrom3D(mol)
    (group,) = mol.GetSubstructMatches(Chem.MolFromSmarts("[N;H2;+0]S(=O)=O"))
    atom = mol.GetAtomWithIdx(group[0])
    atom.SetFormalCharge(-1)
    atom.SetNumExplicitHs(1)
    Chem.SanitizeMol(mol)
    return mol, hashlib.sha256(path.read_bytes()).hexdigest()


def backbone(work, pdb):
    chain = gemmi.read_structure(str(work / "raw" / f"{pdb}.pdb"))[0]["A"]
    return {
        str(r.seqid): list(r.find_atom("CA", "*").pos)
        for r in chain
        if r.het_flag == "A" and r.find_atom("CA", "*")
    }


def verify(work):
    structures = {p: instance(work, p) for p in RESIDUES}
    protein = {p: backbone(work, p) for p in RESIDUES}
    worst = 0.0
    count = 0
    with (work / "holdout_modes.csv").open() as fp:
        key = lambda row: (row["crystal"], row["receptor"], row["seed"], row["arm"])
        for (src, rec, seed, arm), iterator in itertools.groupby(
            csv.DictReader(fp), key
        ):
            rows = list(iterator)
            keys = sorted(protein[src].keys() & protein[rec].keys())
            x, y = (np.array([protein[p][k] for k in keys]) for p in (src, rec))
            transform, _ = Rotation.align_vectors(y - y.mean(0), x - x.mean(0))
            ref = Chem.Mol(structures[src][0])
            conf = ref.GetConformer()
            conf.SetPositions(
                transform.apply(conf.GetPositions() - x.mean(0)) + y.mean(0)
            )
            path = work / "runs" / f"{rec}_{src}_{seed}_{arm}.pdbqt"
            qt = PDBQTMolecule.from_file(str(path), skip_typing=True)
            probe = Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(qt)[0])
            if Chem.MolToSmiles(ref) != Chem.MolToSmiles(probe):
                raise ValueError(
                    f"Crystal vs probe stereochemistry or atom graph differs: {src}: "
                    f"{Chem.MolToSmiles(ref)} vs {Chem.MolToSmiles(probe)}"
                )
            result = [
                rdMolAlign.CalcRMS(probe, ref, prbId=c.GetId())
                for c in probe.GetConformers()
            ]
            if len(result) != len(rows):
                raise ValueError(f"Different MODEL count {path}")
            expected = [float(r["rmsd"]) for r in rows]
            worst = max(worst, max(abs(np.asarray(expected) - result)))
            count += len(result)
    if worst > 1e-5:
        raise ValueError(
            f"Independent RMSD difference exceeds tolerated rounding: {worst}"
        )
    print(
        f"Independently verified {count} poses. Maximum absolute RMSD difference {worst:.8f} angstrom."
    )
    for pdb, (_, digest) in structures.items():
        print(f"{digest}  {pdb}-instance.sdf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    verify(parser.parse_args().work)

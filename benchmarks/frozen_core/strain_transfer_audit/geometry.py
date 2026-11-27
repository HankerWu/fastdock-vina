"""Retrospective geometry audit prompted by an extreme prepared-reference energy."""

import argparse
import tarfile
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import replay
from decompose import read, read_first
from evaluate import physical
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem
from summarize import write


def inspect(mol):
    xyz = mol.GetConformer().GetPositions()
    props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
    deviations = []
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        if min(a.GetAtomicNum(), b.GetAtomicNum()) == 1:
            continue
        i, j = a.GetIdx(), b.GetIdx()
        _, _, optimal = props.GetMMFFBondStretchParams(mol, i, j)
        length = np.linalg.norm(xyz[i] - xyz[j])
        deviations.append((abs(length - optimal), i, j, length, optimal))
    worst, i, j, length, opt = max(deviations)
    return {
        "max_bond_deviation_A": worst,
        "atom1_zero_based": i,
        "atom2_zero_based": j,
        "distance_A": length,
        "mmff_r0_A": opt,
    }


def prepare(smiles, random_coords=False):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    opt = AllChem.ETKDGv3()
    opt.randomSeed = 2026
    opt.useRandomCoords = random_coords
    if AllChem.EmbedMolecule(mol, opt):
        raise ValueError(f"ETKDG failed for {smiles}")
    embedded = inspect(mol)
    if AllChem.MMFFOptimizeMolecule(mol, maxIters=1000):
        raise ValueError(f"MMFF optimization failed for {smiles}")
    (setup,) = MoleculePreparation().prepare(mol)
    text, ok, reason = PDBQTWriterLegacy.write_string(setup, add_index_map=True)
    if not ok:
        raise ValueError(reason)
    return embedded, inspect(mol), text


def main(out):
    rows = []
    inputs = read(replay.previous.BENCH / "affinity_partition/fxa_series.csv")
    with TemporaryDirectory() as tmp:
        work = Path(tmp)
        with tarfile.open(
            replay.previous.FROZEN / "affinity/fxa_frozen.tar.xz"
        ) as archive:
            archive.extractall(work, filter="data")
        for record in inputs:
            if record["included"] != "1":
                continue
            c = record["ligand"]
            text = (work / "2P16" / (c + ".pdbqt")).read_text()
            embedded, optimized, result = prepare(record["docked_smiles"])
            rows.append(dict(compound=c, stage="ETKDG", **embedded))
            rows.append(dict(compound=c, stage="MMFF94_minimized", **optimized))
            if result != text:
                raise ValueError(f"Input preparation did not reproduce: {c}")
            prep = physical.Pose(0, 0, 0, text, None, None)
            prec, _ = physical.decode_ligand(prep)
            rows.append(dict(compound=c, stage="saved_input", **inspect(prec)))
            for receptor in ("1FJS", "2P16"):
                for seed in (271828, 314159):
                    pose = read_first(work / receptor / f"{c}_{seed}.pdbqt")
                    docked, _ = physical.decode_ligand(pose)
                    rows.append(
                        dict(compound=c, stage=f"{receptor}_{seed}", **inspect(docked))
                    )
    write(out / "geometry.csv", rows)
    pooled = [r for r in read(out / "pairs.csv.xz") if r["receptor"] == "both"]
    deletes = []
    for c in sorted({r["compound"] for r in rows}):
        for status in ("incident", "remaining"):
            subset = [
                r
                for r in pooled
                if (c in (r["stronger"], r["weaker"])) == (status == "incident")
            ]
            deletes.append(
                dict(
                    compound=c,
                    subset=status,
                    pairs=len(subset),
                    **{a: sum(float(r[a]) for r in subset) for a in replay.ARMS},
                )
            )
    write(out / "influence.csv", deletes)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    main(parser.parse_args().out)

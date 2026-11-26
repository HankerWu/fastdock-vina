"""Recompute receptor transforms and all deposited-pose RMSDs by a second path."""

import argparse
import csv
import hashlib
import json
import urllib.request
from collections import defaultdict
from pathlib import Path

import gemmi
import numpy as np
from meeko import PDBQTMolecule, RDKitMolCreate
from prepare import PANEL
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from scipy.spatial.transform import Rotation


def extract(path, chain, marker):
    residues = gemmi.read_structure(str(path))[0][chain]
    pocket = []
    ca = {}
    for res in residues:
        if res.name == marker:
            for atom in res:
                if atom.element.name not in ("H", "D") and atom.altloc in ("\x00", "A"):
                    pocket.append(list(atom.pos))
        elif res.het_flag == "A":
            atom = res.find_atom("CA", "*")
            if atom:
                ca[(res.seqid.num, res.seqid.icode, res.name)] = list(atom.pos)
    return ca, np.array(pocket)


def independent_reference(raw, q, c, cc, rec, rc, marker):
    sdf = raw / (q + "-instance.sdf")
    if not sdf.exists():
        res = next(
            res
            for res in gemmi.read_structure(str(raw / (q + ".pdb")))[0][c]
            if res.name == cc
        )
        url = f"https://models.rcsb.org/v1/{q.lower()}/ligand?auth_asym_id={c}&auth_seq_id={res.seqid.num}&encoding=sdf"
        sdf.write_bytes(urllib.request.urlopen(url, timeout=60).read())
    m = Chem.SDMolSupplier(str(sdf), removeHs=True)[0]
    if rec != q:
        src, _ = extract(raw / (q + ".pdb"), c, cc)
        dst, site = extract(raw / (rec + ".pdb"), rc, marker)
        keys = [
            k
            for k in sorted(src.keys() & dst.keys())
            if np.linalg.norm(site - dst[k], axis=1).min() <= 12
        ]
        a, b = (np.array([d[k] for k in keys]) for d in (src, dst))
        rot, _ = Rotation.align_vectors(b - b.mean(0), a - a.mean(0))
        cf = m.GetConformer()
        cf.SetPositions(rot.apply(cf.GetPositions() - a.mean(0)) + b.mean(0))
    return m


def verify(work, tables):
    bykey = defaultdict(dict)
    with (tables / "poses.csv").open() as f:
        for r in csv.DictReader(f):
            key = (r["target"], r["receptor"], r["seed"], r["slot"], r["state"])
            bykey[key][int(r["model"])] = float(r["rmsd"])
    max_rms = 0.0
    models = pools = 0
    hashes = []
    for t, cc, q, c, r, rc, marker in PANEL:
        for receiver, chain, lig in [(q, c, cc), (r, rc, marker)]:
            if receiver is None:
                continue
            files = sorted((work / "pools" / t / receiver).glob("*.pdbqt"))
            if not files:
                continue
            ref = independent_reference(work / "raw", q, c, cc, receiver, chain, lig)
            for path in files:
                seed, slot, state = path.stem.split("_")
                expected = bykey[t, receiver, seed, slot, state]
                pdbqt = PDBQTMolecule.from_file(str(path), skip_typing=True)
                m = Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(pdbqt)[0])
                assert m.GetNumAtoms() == ref.GetNumAtoms()
                if Chem.MolToSmiles(ref) != Chem.MolToSmiles(m):
                    raise ValueError(f"{q}: different reconstructed ligand graphs")
                errs = [
                    abs(rdMolAlign.CalcRMS(m, ref, prbId=conf.GetId()) - expected[j])
                    for j, conf in enumerate(m.GetConformers(), 1)
                ]
                max_rms = max(max_rms, max(errs))
                models += len(errs)
                pools += 1
                hashes.append(
                    f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(work)}"
                )
    if max_rms > 1e-4:
        raise ValueError(f"RMSD discrepancy {max_rms}")
    (tables / "pools.sha256").write_text("\n".join(hashes) + "\n")
    print(
        json.dumps({"pools": pools, "models": models, "max_RMSD_difference_A": max_rms})
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("work", type=Path)
    p.add_argument("tables", type=Path)
    args = p.parse_args()
    verify(args.work, args.tables)

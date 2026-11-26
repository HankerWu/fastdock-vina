"""Recompute every RMSD from RCSB ligand-instance SDF, Meeko and RDKit CalcRMS."""
import argparse
import csv
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from meeko import PDBQTMolecule, RDKitMolCreate
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from scipy.spatial.transform import Rotation

COMPONENT = {"181L": "BNZ", "182L": "BZF", "183L": "DEN", "184L": "I4B",
             "185L": "IND", "186L": "N4B", "187L": "PXY", "188L": "OXE", "1NHB": "PYJ",
             "4OGI": "R78", "4OGJ": "2TA", "4O76": "1M3", "4O77": "2RE"}


def chain(pdb):
    return "B" if pdb == "4OGJ" else "A"


def reference(work, pdb):
    path = work / "raw" / (pdb+"-instance.sdf")
    if not path.exists():
        line = next(x for x in (work/"raw"/(pdb+".pdb")).read_text().splitlines()
                    if x.startswith("HETATM") and x[17:20] == COMPONENT[pdb] and x[21] == chain(pdb))
        url = (f"https://models.rcsb.org/v1/{pdb.lower()}/ligand?auth_asym_id={chain(pdb)}"
               f"&auth_seq_id={int(line[22:26])}&encoding=sdf")
        with urllib.request.urlopen(url, timeout=60) as response:
            data = response.read()
        if b"M  END" not in data:
            raise ValueError(f"No ligand-instance molecule returned for {pdb}")
        path.write_bytes(data)
    crystal = Chem.SDMolSupplier(str(path), removeHs=True)[0]
    if pdb in ("4OGI", "4OGJ", "4O76"):
        # Match the predefined ligand preparation microstate, preserving the heavy graph.
        site, = crystal.GetSubstructMatches(Chem.MolFromSmarts("[N;X3;R;!$(N-C=O);!$(N-c)]"))
        amine = crystal.GetAtomWithIdx(site[0])
        amine.SetFormalCharge(1)
        amine.SetNumExplicitHs(1)
        Chem.SanitizeMol(crystal)
    return crystal


def ca(work, pdb):
    result = {}
    for x in (work/"raw"/(pdb+".pdb")).read_text().splitlines():
        if x.startswith("ATOM") and x[21] == chain(pdb) and x[16] in (" ", "A") and x[12:16].strip() == "CA":
            result[int(x[22:26])] = [float(x[k:k+8]) for k in (30, 38, 46)]
    return result


def verify(work):
    with (work/"cells.csv").open() as table:
        rows = list(csv.DictReader(table))
    ids = sorted({r["ligand"] for r in rows})
    with ThreadPoolExecutor(5) as pool:
        refs = dict(zip(ids, pool.map(lambda pdb: reference(work, pdb), ids)))
    deviations = []
    for row in rows:
        i, r, seed, arm = (row[k] for k in ("ligand", "receptor", "seed", "arm"))
        ref_ca, test_ca = ca(work, i), ca(work, r)
        if i in ("4OGI", "4OGJ", "4O76", "4O77"):
            keys = sorted(k for k in ref_ca.keys() & test_ca.keys() if 44 <= k <= 166)
        else:
            keys = [k for k in range(80, 161) if k not in range(107, 116)]
        a, b = (np.array([atoms[k] for k in keys]) for atoms in (ref_ca, test_ca))
        rotation, _ = Rotation.align_vectors(b-b.mean(0), a-a.mean(0))
        crystal = Chem.Mol(refs[i])
        conf = crystal.GetConformer()
        conf.SetPositions(rotation.apply(conf.GetPositions()-a.mean(0))+b.mean(0))
        pdbqt = PDBQTMolecule.from_file(str(work/"runs"/f"{i}_{r}_{seed}_{arm}.pdbqt"), skip_typing=True)
        probe = Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(pdbqt)[0])
        if Chem.MolToSmiles(probe) != Chem.MolToSmiles(crystal):
            raise ValueError(f"Independent graph mismatch in {i}: {Chem.MolToSmiles(crystal)} versus {Chem.MolToSmiles(probe)}")
        computed = [rdMolAlign.CalcRMS(probe, crystal, prbId=c.GetId()) for c in probe.GetConformers()]
        expected = [float(x) for x in row["pool_rmsds"].split(",")]
        deviations.extend(abs(np.asarray(computed)-expected))
    error = max(deviations)
    print(f"{len(rows)} runs, {len(deviations)} models, worst absolute RMSD difference {error:.8f} A")
    if error > 0.0001:
        raise ValueError("RMSD discrepancy exceeds coordinate conversion/rounding tolerance")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    verify(p.parse_args().work)

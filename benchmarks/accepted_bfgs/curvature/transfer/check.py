"""Independent raw-deposit C-alpha fit and Meeko/RDKit fixed-frame minimum RMSD."""

import argparse
import csv
import json
from pathlib import Path

import gemmi
import numpy as np
from meeko import PDBQTMolecule, RDKitMolCreate
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from scipy.spatial.transform import Rotation
from scipy.spatial.distance import cdist


def crystal(work, pdb):
    meta = json.loads((work / pdb / "metadata.json").read_text())
    path = work / "raw" / (pdb + ".pdb")
    structure = gemmi.read_structure(str(path))
    chain = structure[0]["C"]
    ligand = next(r for r in chain if r.name == meta["component"])
    points = np.array(
        [
            [a.pos.x, a.pos.y, a.pos.z]
            for a in ligand
            if a.element.name not in {"H", "D"}
        ]
    )
    ca = {}
    for ch in structure[0]:
        if ch.name not in "ABC":
            continue
        for res in ch:
            if res.name == meta["component"]:
                continue
            atoms = [a for a in res if a.name == "CA" and a.altloc in ("\0", "A")]
            if not atoms:
                continue
            pos = atoms[0].pos
            v = [pos.x, pos.y, pos.z]
            if cdist(points, [v]).min() <= 12:
                ca[(ch.name, str(res.seqid), res.name)] = v
    lines = path.read_text().splitlines()
    lig = [
        l
        for l in lines
        if l.startswith("HETATM") and l[21] == "C" and l[17:20] == meta["component"]
    ]
    ids = {int(l[6:11]) for l in lig}
    conect = [l for l in lines if l.startswith("CONECT") and int(l[6:11]) in ids]
    raw = Chem.MolFromPDBBlock(
        "\n".join([*lig, *conect, "END"]) + "\n",
        sanitize=False,
        removeHs=False,
        proximityBonding=False,
    )
    ref = AllChem.AssignBondOrdersFromTemplate(Chem.MolFromSmiles(meta["smiles"]), raw)
    Chem.AssignStereochemistryFrom3D(ref, replaceExistingTags=True)
    assert Chem.MolToInchiKey(ref) == meta["inchikey"]
    return ca, ref, meta["inchikey"]


def verify(work):
    with (work / "modes.csv").open() as f:
        rows = list(csv.DictReader(f))
    refdata = {k: crystal(work, k) for k in ("4W9H", "5LLI")}
    maxima = 0.0
    count = 0
    for path in sorted((work / "runs").glob("*.pdbqt")):
        q, r, se, arm = path.stem.split("_", 3)
        a, b = refdata[q][0], refdata[r][0]
        keys = sorted(a.keys() & b.keys())
        x, y = [np.array([c[k] for k in keys]) for c in (a, b)]
        rot, _ = Rotation.align_vectors(y - y.mean(0), x - x.mean(0))
        ref = Chem.Mol(refdata[q][1])
        conf = ref.GetConformer()
        moved = rot.apply(conf.GetPositions() - x.mean(0)) + y.mean(0)
        for i, pos in enumerate(moved):
            conf.SetAtomPosition(i, pos)
        pred = Chem.RemoveHs(
            RDKitMolCreate.from_pdbqt_mol(
                PDBQTMolecule.from_file(str(path), skip_typing=True)
            )[0]
        )
        assert Chem.MolToInchiKey(pred) == refdata[q][2]
        matches = ref.GetSubstructMatches(
            pred, uniquify=False, useChirality=True, maxMatches=10001
        )
        assert 0 < len(matches) < 10001
        maps = [list(enumerate(m)) for m in matches]
        expected = [
            v
            for v in rows
            if (v["query"], v["receptor"], v["seed"], v["arm"]) == (q, r, se, arm)
        ]
        assert len(expected) == pred.GetNumConformers()
        for rank, row in enumerate(expected):
            error = rdMolAlign.CalcRMS(pred, ref, prbId=rank, refId=0, map=maps)
            maxima = max(maxima, abs(error - float(row["rmsd"])))
            copy = Chem.Mol(pred)
            Chem.AssignStereochemistryFrom3D(
                copy, confId=rank, replaceExistingTags=True
            )
            assert Chem.MolToInchiKey(copy) == refdata[q][2]
            count += 1
    assert count == len(rows) and maxima < 1e-5, (count, len(rows), maxima)
    return {
        "modes": count,
        "max_rmsd_abs_difference_A": maxima,
        "pose_chiral_graph_checks": count,
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("work", type=Path)
    print(json.dumps(verify(ap.parse_args().work), indent=2))

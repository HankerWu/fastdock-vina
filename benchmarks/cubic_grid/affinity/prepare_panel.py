"""Freeze exact Davis Kd rows for WT SYK and FAK using verified compound identities."""

import csv
import hashlib
import json
import sys
from pathlib import Path
import xlrd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

repo, raw, out = [Path(p).resolve() for p in sys.argv[1:]]
sys.path.insert(0, str(repo / "benchmarks/kinase_state_affinity/source"))
import freeze
from prepare import ligand


def write_table(path, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0])
        w.writeheader()
        w.writerows(rows)


(out / "prepared").mkdir(parents=True, exist_ok=True)
src = json.loads(
    (repo / "benchmarks/kinase_state_affinity/source/sources.json").read_text()
)
for fname in ("table1.xls", "table4.xls", "identities.json", "ligands_iso.txt"):
    assert (
        hashlib.sha256((raw / fname).read_bytes()).hexdigest() == src[fname]["sha256"]
    )
ids = json.loads((raw / "identities.json").read_text())["PropertyTable"]["Properties"]
keys = {r["CID"]: r["InChIKey"] for r in ids}
smis = json.loads((raw / "ligands_iso.txt").read_text())
known = {n: (c, smis[str(c)], keys[c]) for n, c in freeze.COMPOUNDS.items()}
for name, rec in json.loads(
    (repo / "benchmarks/torsion_kinase_affinity/names.json").read_text()
).items():
    (r,) = rec["record"]
    known[name] = r["CID"], r["SMILES"], r["InChIKey"]
s1, s4 = [
    xlrd.open_workbook(raw / f).sheet_by_index(0) for f in ("table1.xls", "table4.xls")
]
panel, excluded = [], []
for target in ("SYK", "FAK"):
    (ident,) = [r for r in s1.get_rows() if r[2].value == target]
    (row,) = [r for r in s4.get_rows() if r[2].value == target]
    assert ident[0].value == row[0].value and ident[3].value == "NO"
    for j, name in enumerate(s4.row_values(0)[3:], 3):
        if not row[j].value:
            continue
        if name not in known:
            excluded.append(
                dict(
                    target=target,
                    name=name,
                    reason="outside fixed name-structure dictionary",
                )
            )
            continue
        cid, smi, key = known[name]
        m = Chem.MolFromSmiles(smi)
        assert Chem.MolToInchiKey(m) == key
        mol = max(Chem.GetMolFrags(m, asMols=True), key=lambda f: f.GetNumHeavyAtoms())
        if any(
            x.specified == Chem.StereoSpecified.Unspecified
            for x in Chem.FindPotentialStereo(mol)
        ):
            excluded.append(
                dict(target=target, name=name, reason="unspecified stereochemistry")
            )
            continue
        panel.append(
            dict(
                target=target,
                accession=ident[0].value,
                name=name,
                cid=str(cid),
                kd_nm=float(row[j].value),
                heavy=mol.GetNumHeavyAtoms(),
                inchi_key=key,
                source_smiles=smi,
                smiles=Chem.MolToSmiles(mol),
                scaffold=MurckoScaffold.MurckoScaffoldSmiles(
                    mol=mol, includeChirality=False
                ),
            )
        )
write_table(out / "panel.csv", panel)
write_table(out / "exclusions.csv", excluded)
for r in {r["cid"]: r for r in panel}.values():
    ligand(out, r)

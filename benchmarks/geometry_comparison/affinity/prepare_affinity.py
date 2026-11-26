"""Prepare exact same-assay BACE1 Ki and MAPK14 Kd rows for optimizer comparisons."""
import argparse
import csv
import json
import sys
import tarfile
from pathlib import Path

import xlrd
from scipy.spatial.distance import cdist
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "benchmarks/kinase_state_affinity/source"))
from prepare import ligand

sys.path.insert(0, str(ROOT / "benchmarks/torsion_kinase_affinity"))
from prepare_receptors import prepare as prepare_receptor

DAVIS_CIDS = {
    "AST-487": 11409972, "BIBW-2992": 10184653, "BIRB-796": 156422,
    "Dasatinib": 3062316, "EXEL-2880/GSK-1363089": 42642645, "Nilotinib": 644241,
    "PD-173955": 447077, "SB-203580": 176155, "Sorafenib": 216239, "VX-745": 3038525,
}
EXTERNAL = {"ALK": ("7BTT", "F8R"), "RET": ("7DUA", "HJ0"), "TNK2": ("7KP6", "WTP")}
ASSAYS = {"CHEMBL2010674": "CHEMBL2015090", "CHEMBL2150878": "CHEMBL2155063"}


def table(path, rows):
    with path.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def case(target, assay, endpoint, compound, name, nm, source_smiles, mol):
    return dict(target=target, assay=assay, endpoint=endpoint, compound=compound, name=name,
                value_nM=float(nm), source_smiles=source_smiles, smiles=Chem.MolToSmiles(mol),
                formal_charge=Chem.GetFormalCharge(mol), scaffold=MurckoScaffold.MurckoScaffoldSmiles(mol=mol))


def cases(work):
    raw = work / "raw"
    rows, exclusions = [], []
    tables = [xlrd.open_workbook(raw / f).sheet_by_index(0) for f in ("table1.xls", "table4.xls")]
    identities, data = [[r for r in sheet.get_rows() if r[2].value == "p38-alpha"][0] for sheet in tables]
    assert data[0].value == identities[0].value == "NP_620581.1" and identities[3].value == "NO"
    lookup = json.loads((raw / "ligands_iso.txt").read_text())
    for j, name in enumerate(tables[1].row_values(0)[3:], 3):
        if not data[j].value:
            continue
        cid = str(DAVIS_CIDS[name])
        smi = lookup[cid]
        mol = Chem.MolFromSmiles(smi)
        if any(x.specified == Chem.StereoSpecified.Unspecified for x in Chem.FindPotentialStereo(mol)):
            exclusions.append(dict(assay="Davis-p38-alpha", compound=cid, reason="unspecified stereochemistry"))
            continue
        rows.append(case("MAPK14", "Davis-p38-alpha", "Kd", cid, name, data[j].value, smi, mol))


    from freeze import COMPOUNDS

    all_names = {**COMPOUNDS, **DAVIS_CIDS}
    all_smiles = dict(lookup)
    for name, entry in json.loads((ROOT / "benchmarks/torsion_kinase_affinity/names.json").read_text()).items():
        obj, = entry["record"]
        assert Chem.MolToInchiKey(Chem.MolFromSmiles(obj["SMILES"])) == obj["InChIKey"]
        all_names[name] = obj["CID"]
        all_smiles[str(obj["CID"])] = obj["SMILES"]
    for target in EXTERNAL:
        ident, = [r for r in tables[0].get_rows() if r[2].value == target]
        data, = [r for r in tables[1].get_rows() if r[2].value == target]
        assert data[0].value == ident[0].value and ident[3].value == "NO"
        assay = "Davis-" + target
        eligible = []
        for j, name in enumerate(tables[1].row_values(0)[3:], 3):
            if not data[j].value:
                continue
            if name not in all_names:
                exclusions.append(dict(assay=assay, compound=name, reason="outside verified identity dictionary"))
                continue
            cid = str(all_names[name])
            smi = all_smiles[cid]
            mol = max(Chem.GetMolFrags(Chem.MolFromSmiles(smi), asMols=True), key=lambda m: m.GetNumHeavyAtoms())
            if any(x.specified == Chem.StereoSpecified.Unspecified for x in Chem.FindPotentialStereo(mol)):
                exclusions.append(dict(assay=assay, compound=cid, reason="unspecified stereochemistry"))
                continue
            eligible.append(case(target, assay, "Kd", cid, name, data[j].value, smi, mol))
        eligible.sort(key=lambda r: int(r["compound"]))
        rows.extend(eligible[:8])
        exclusions.extend(dict(assay=assay, compound=r["compound"], reason="beyond first eight eligible IDs")
                          for r in eligible[8:])

    # Predeclare a single cationic cyclic guanidinium microstate. This removes
    # the artificial exocyclic C=NH E/Z designation of a neutral drawing.
    guanidine = Chem.MolFromSmarts("[NX2;+0]=[C]([N])[N]")
    for doc, assay in ASSAYS.items():
        accepted = []
        archive = json.loads((raw / f"{doc}_Ki.json").read_text())
        for rec in archive["activities"]:
            if rec["assay_chembl_id"] != assay:
                continue
            cid = rec["molecule_chembl_id"]
            if (rec["standard_type"], rec["standard_units"], rec["standard_relation"]) != ("Ki", "nM", "="):
                exclusions.append(dict(assay=assay, compound=cid, reason="nonexact Ki"))
                continue
            assert rec["target_pref_name"] == "Beta-secretase 1" and rec["target_organism"] == "Homo sapiens"
            smi = rec["canonical_smiles"]
            mol = Chem.MolFromSmiles(smi)
            match, = mol.GetSubstructMatches(guanidine)
            atom = mol.GetAtomWithIdx(match[0])
            atom.SetNumExplicitHs(atom.GetTotalNumHs() + 1)
            atom.SetNoImplicit(True)
            atom.SetFormalCharge(1)
            Chem.SanitizeMol(mol)
            # In the one deposited amino tautomer, shift the positive double bond
            # to terminal NH2. All entries then use the same guanidinium resonance.
            if atom.IsInRing():
                central = match[1]
                terminal, = [a for a in (match[2], match[3]) if not mol.GetAtomWithIdx(a).IsInRing()]
                editable = Chem.RWMol(mol)
                editable.GetBondBetweenAtoms(match[0], central).SetBondType(Chem.BondType.SINGLE)
                editable.GetBondBetweenAtoms(central, terminal).SetBondType(Chem.BondType.DOUBLE)
                editable.GetAtomWithIdx(match[0]).SetFormalCharge(0)
                dest = editable.GetAtomWithIdx(terminal)
                dest.SetFormalCharge(1)
                dest.SetNumExplicitHs(2)
                mol = editable.GetMol()
                Chem.SanitizeMol(mol)
            if any(x.specified == Chem.StereoSpecified.Unspecified for x in Chem.FindPotentialStereo(mol)):
                exclusions.append(dict(assay=assay, compound=cid, reason="unspecified stereochemistry"))
                continue
            accepted.append(case("BACE1", assay, "Ki", cid, cid, rec["standard_value"], smi, mol))
        # Bound run time without selecting on potency or docking results.
        accepted.sort(key=lambda r: int(r["compound"].removeprefix("CHEMBL")))
        rows.extend(accepted[:8])
        exclusions.extend(dict(assay=assay, compound=r["compound"], reason="beyond first eight eligible IDs")
                          for r in accepted[8:])
    return rows, exclusions



def external_structure(work, pdb, ccd):
    # Select the ligand-bound protomer used in the independent PoseBusters
    # crystal data. Native PDB coordinates (including omitted chains and ions)
    # are the source of truth for the archived preparation.
    selected = {"7BTT": "A", "7DUA": "B", "7KP6": "A"}
    chain = selected[pdb]
    source = work / "raw" / (pdb + ".pdb")
    records = source.read_text().splitlines()
    xyz = lambda line: [float(line[k:k + 8]) for k in (30, 38, 46)]
    atoms = [l for l in records if l.startswith(("ATOM  ", "HETATM"))
             and l[76:78].strip() not in ("H", "D") and l[16] in " A"]
    lig = [l for l in atoms if l[17:20].strip() == ccd and l[21] == chain]
    omitted = [xyz(l) for l in atoms if l[17:20] not in ("DOD", "HOH")
               and not (l[21] == chain and (l.startswith("ATOM  ") or l[17:20].strip() == ccd))]
    closest = float(cdist([xyz(l) for l in lig], omitted).min()) if omitted else None
    if closest is not None and closest < 8:
        raise ValueError(f"{pdb}: omitted nonwater atom within 8 A ({closest})")
    # The shared preparer accepts chain A. Build an isolated copy with the
    # chosen chain aliased to A, preserving the unmodified original.
    import shutil
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        scratch = Path(tmp)
        (scratch / "raw").mkdir()
        view = [l[:21] + "A" + l[22:] for l in atoms if l[21] == chain]
        (scratch / "raw" / (pdb + ".pdb")).write_text("\n".join(view + ["END", ""]))
        prepare_receptor(scratch, pdb, ccd)
        for path in (scratch / "prepared").iterdir():
            shutil.copyfile(path, work / "prepared" / path.name)
    meta = json.loads((work / "prepared" / f"{pdb}.json").read_text())
    if meta["closest_repaired"] is not None and meta["closest_repaired"] < 6:
        raise ValueError(f"{pdb}: repaired atom within 6 A")
    meta.update(source_chain=chain, nearest_omitted_nonwater_A=closest)
    (work / "prepared" / f"{pdb}.json").write_text(json.dumps(meta, indent=2) + "\n")
    return meta["center"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    args = parser.parse_args()
    work = args.work.resolve()
    prepared = work / "prepared"
    prepared.mkdir(parents=True, exist_ok=True)
    rows, exclusions = cases(work)
    table(work / "panel.csv", rows)
    table(work / "exclusions.csv", exclusions)
    for row in {r["compound"]: r for r in rows}.values():
        ligand(work, {"cid": row["compound"], "smiles": row["smiles"]})
    for pdb, ccd in (("4DJW", "0KP"), ("4GID", "0GH")):
        prepare_receptor(work, pdb, ccd)
        center = json.loads((prepared / f"{pdb}.json").read_text())["center"]
        (prepared / f"{pdb}.box").write_text("".join(f"center_{a} = {v}\nsize_{a} = 26\n" for a, v in zip("xyz", center)))
    for pdb, ccd in EXTERNAL.values():
        center = external_structure(work, pdb, ccd)
        (prepared / f"{pdb}.box").write_text("".join(f"center_{a} = {v}\nsize_{a} = 26\n" for a, v in zip("xyz", center)))
    with tarfile.open(ROOT / "benchmarks/ring_ensemble/fixtures/structures-prepared.tar.gz") as tar:
        for pdb in ("3GCP", "1KV2"):
            for source, dest in (("receptor.pdbqt", f"{pdb}.pdbqt"), ("box.txt", f"{pdb}.box")):
                (prepared / dest).write_bytes(tar.extractfile(f"prepared/p38a/{pdb}/{source}").read())


if __name__ == "__main__":
    main()

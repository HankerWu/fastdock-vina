"""Rebuild same-assay FXa comparisons directly from ChEMBL JSON and docked coordinates."""

import argparse
import itertools
import json
import re
import statistics
from decimal import Decimal
from pathlib import Path

from rdkit import Chem


def audit(activity, folders):
    series = json.loads(activity.read_text())["activities"]
    labels = {}
    for row in series:
        if (
            row["assay_chembl_id"],
            row["standard_type"],
            row["standard_relation"],
            row["standard_units"],
        ) != ("CHEMBL661804", "Ki", "=", "nM"):
            continue
        mol = Chem.MolFromSmiles(row["canonical_smiles"])
        if any(
            tag == "?"
            for _, tag in Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        ):
            continue
        key = row["molecule_chembl_id"]
        if key in labels:
            raise ValueError(f"Duplicate experimental measurement for {key}")
        labels[key] = Decimal(row["standard_value"])
    pairs = [
        (a, b)
        for a, b in itertools.permutations(labels, 2)
        if labels[b] >= 3 * labels[a]
    ]
    predictions = {}
    chemistry = {}
    for folder in folders:
        for ligand in labels:
            for seed in (271828, 314159):
                source = folder / f"{ligand}_{seed}.pdbqt"
                text = source.read_text().split("ENDMDL", 1)[0]

                def energy(prefix):
                    return Decimal(
                        re.search(
                            "^REMARK " + re.escape(prefix) + r":\s+(\S+)",
                            text,
                            re.MULTILINE,
                        )[1]
                    )

                mol = Chem.MolFromSmiles(
                    re.search(r"^REMARK SMILES (?!IDX )(\S+)$", text, re.MULTILINE)[1]
                )
                state = (mol.GetNumHeavyAtoms(), Chem.GetFormalCharge(mol))
                if chemistry.setdefault(ligand, state) != state:
                    raise ValueError(f"Inconsistent inputs for {ligand}")
                predictions[folder, seed, ligand] = (
                    energy("VINA RESULT"),
                    energy("INTER + INTRA") - energy("UNBOUND"),
                )

    def credit(gap):
        return 1 if gap > 0 else 0.5 if gap == 0 else 0

    for subset in [[r] for r in folders] + ([folders] if len(folders) > 1 else []):
        groups = [(f, s) for f in subset for s in (271828, 314159)]
        vals = []
        for a, b in pairs:
            gaps = [
                statistics.median(
                    predictions[f, s, b][i] - predictions[f, s, a][i] for f, s in groups
                )
                for i in range(2)
            ]
            vals.append(
                (
                    chemistry[a] == chemistry[b],
                    credit(gaps[0]),
                    credit(gaps[1]),
                    credit(chemistry[a][0] - chemistry[b][0]),
                )
            )
        for matched in (False, True):
            chosen = [r[1:] for r in vals if not matched or r[0]]
            print(
                [f.name for f in subset],
                "size_charge_matched",
                matched,
                "pairs",
                len(chosen),
                "native/raw/size",
                list(map(sum, zip(*chosen, strict=True))),
            )
    print(
        f"{len(labels)} independent experimental labels; {len(predictions)} raw top models"
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("activity", type=Path)
    p.add_argument("outputs", nargs="+", type=Path)
    a = p.parse_args()
    audit(a.activity, a.outputs)

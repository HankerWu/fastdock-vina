"""Recompute order statistics from first MODEL records without reading primary score tables."""

import argparse
import csv
import itertools
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

# Transcribed independently from the discussion and compound characterizations,
# Holt et al., JACS 1993, 115, 9925–9938.
EXPECTED = {
    4: (186, 186, "C22H37NO4"),
    5: (110, 110, "C22H31NO4"),
    6: (12, 12, "C25H37NO7"),
    7: (250, 250, "C27H39NO4"),
    8: (10, 10, "C28H35NO4"),
    9: (7, 7, "C28H41NO4"),
    10: (300, 600, "C28H35NO4"),
    11: (300, 300, "C28H41NO4"),
}


def run(work):
    output = {}
    for tag in ("affinity", "affinity_confirmation"):
        root = work / tag
        rows = list(csv.DictReader((root / "ligands.csv").open()))
        for r in rows:
            ident = int(r["compound"])
            molecule = Chem.MolFromSmiles(r["smiles"])
            assert (
                int(r["ki_low_nm"]),
                int(r["ki_high_nm"]),
                CalcMolFormula(molecule),
            ) == EXPECTED[ident]
        samples = defaultdict(list)
        for fp in (root / "runs").glob("*.pdbqt"):
            comp, rec, _seed, arm = fp.stem.split("_")
            first = fp.read_text().split("ENDMDL")[0].splitlines()
            score = [
                float(l.split()[3])
                for l in first
                if l.startswith("REMARK VINA RESULT:")
            ]
            assert len(score) == 1
            samples[arm, int(comp), rec].append(score[0])
        rankings = {}
        for arm in ("native", "cubic"):
            mean = {
                i: np.mean([np.median(samples[arm, i, r]) for r in ("1FKG", "1FKH")])
                for i in EXPECTED
            }
            rankings[arm] = mean
        pairs = [
            (i, j)
            for i, j in itertools.permutations(EXPECTED, 2)
            if EXPECTED[j][0] / EXPECTED[i][1] >= 3
        ]
        credits = {
            arm: [int(mean[i] < mean[j]) for i, j in pairs]
            for arm, mean in rankings.items()
        }
        checked = list(csv.DictReader((root / "pairs.csv").open()))
        for r in checked:
            i, j = int(r["stronger"]), int(r["weaker"])
            for arm in ("native", "cubic"):
                assert int(r[arm + "_credit"]) == credits[arm][pairs.index((i, j))]
        output[tag] = {
            "n_unique_pairs": len(pairs),
            "n_unique_compounds": len(EXPECTED),
            "correct": {arm: sum(credit) for arm, credit in credits.items()},
            "gains": sum(a < b for a, b in zip(credits["native"], credits["cubic"])),
            "losses": sum(a > b for a, b in zip(credits["native"], credits["cubic"])),
            "matched_epimers": {
                f"{i}/{j}": {
                    arm: int(mean[i] < mean[j]) for arm, mean in rankings.items()
                }
                for i, j in [(8, 10), (9, 11)]
            },
        }
    return output


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    a = p.parse_args()
    result = run(a.work)
    (a.work / "affinity_independent.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print(json.dumps(result))

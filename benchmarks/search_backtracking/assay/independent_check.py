"""Recalculate Davis Kd labels, chemical identities and paired ranks from raw files."""

import argparse
import csv
import hashlib
import itertools
import json
from collections import defaultdict
from decimal import Decimal
from pathlib import Path

import numpy as np
import xlrd
from rdkit import Chem
from scipy.spatial.distance import cdist
from scipy.stats import kendalltau, spearmanr

HERE = Path(__file__).resolve().parent


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def xyz(text):
    return np.asarray(
        [
            [float(s[i : i + 8]) for i in (30, 38, 46)]
            for s in text.splitlines()
            if s[:6].strip() in {"ATOM", "HETATM"}
            and s[77:].strip() not in {"H", "HD", "HS"}
        ]
    )


def main(work):
    p = read(work / "panel.csv")
    labels = {}
    for name, src in json.loads((HERE / "sources.json").read_text()).items():
        assert (
            hashlib.sha256((work / "raw" / name).read_bytes()).hexdigest()
            == src["sha256"]
        )
    raw_graphs = json.loads((work / "raw/ligands_iso.txt").read_text())
    for target in ("EGFR", "p38-alpha"):

        def original_row(n, target):
            sheet = xlrd.open_workbook(work / "raw" / f"table{n}.xls").sheet_by_index(0)
            (row,) = [
                sheet.row_values(i)
                for i in range(1, sheet.nrows)
                if sheet.cell_value(i, 2) == target
            ]
            return dict(zip(sheet.row_values(0), row, strict=True))

        assay, prot = original_row(4, target), original_row(1, target)
        assert (
            prot["Mutant"] == "NO"
            and assay["Accession Number"] == prot["Accession Number"]
        )
        for r in (r for r in p if r["target"] == target):
            kd = Decimal(str(assay[r["compound"]]))
            assert (
                kd == Decimal(r["kd_nm"])
                and r["assay_accession"] == prot["Accession Number"]
            )
            cid = r["cid"]
            raw = (work / "prepared" / (cid + ".pdbqt")).read_bytes()
            assert hashlib.sha256(raw).hexdigest() == r["ligand_sha256"]
            (graph_line,) = [
                s
                for s in raw.decode().splitlines()
                if s.startswith("REMARK SMILES ")
                and not s.startswith("REMARK SMILES IDX ")
            ]
            mol = Chem.MolFromSmiles(graph_line.split()[2])
            fragments = Chem.GetMolFrags(
                Chem.MolFromSmiles(raw_graphs[cid]), asMols=True
            )
            biggest = max(fragments, key=lambda m: m.GetNumHeavyAtoms())
            assert (
                Chem.MolToInchiKey(biggest)
                == Chem.MolToInchiKey(mol)
                == Chem.MolToInchiKey(Chem.MolFromSmiles(r["smiles"]))
            )
            heavy, charge = mol.GetNumHeavyAtoms(), Chem.GetFormalCharge(mol)
            assert (heavy, charge) == (int(r["heavy"]), int(r["charge"]))
            labels[target, cid] = (kd, heavy, charge)
    report = {
        "target_ligand_rows": len(p),
        "distinct_compounds": len({r["cid"] for r in p}),
        "seeds": {},
    }
    for seed, study, tables in [
        (458921, "experiments", "results"),
        (703109, "cost_check", "cost_check"),
    ]:
        blocks = defaultdict(dict)
        poses = {}
        raw_runs = read(work / study / "runs.csv")
        for r in raw_runs:
            assert int(r["seed"]) == seed
            path = work / study / "raw" / (r["key"] + ".pdbqt")
            models = path.read_text().split("ENDMDL")[:-1]
            (score_line,) = [
                s for s in models[0].splitlines() if s.startswith("REMARK VINA RESULT:")
            ]
            score = Decimal(score_line.split()[3])
            assert len(models) == int(r["modes"]) and score == Decimal(r["first_score"])
            key = (r["panel"], r["cap"], r["arm"])
            assert r["compound"] not in blocks[key]
            blocks[key][r["compound"]] = score
            poses[r["key"]] = (score, models[0], r["receptor"])
        required = {
            (r["target"], r["cap"], r["arm"], r["subset"]): r
            for r in read(HERE / tables / "summary.csv")
        }
        ranks = {
            (r["target"], r["cap"], r["arm"]): r
            for r in read(HERE / tables / "cost_ranks.csv")
        }
        all_credits, matched_counts = {}, {}
        for (target, cap, arm), predictions in sorted(blocks.items()):
            keys = sorted(predictions)
            kd = [float(labels[target, c][0]) for c in keys]
            pred = [float(predictions[c]) for c in keys]
            ref = ranks[target, cap, arm]
            assert abs(kendalltau(kd, pred).statistic - float(ref["kendall_b"])) < 1e-12
            assert abs(spearmanr(kd, pred).statistic - float(ref["spearman"])) < 1e-12
            credits = defaultdict(list)
            for a, b in itertools.combinations(keys, 2):
                ka, ha, qa = labels[target, a]
                kb, hb, qb = labels[target, b]
                if max(ka, kb) < 3 * min(ka, kb):
                    continue
                sign = (ka - kb) * (predictions[a] - predictions[b])
                credit = 1 if sign > 0 else 0.5 if sign == 0 else 0
                credits["all"].append(credit)
                if (ha, qa) == (hb, qb):
                    credits["same_size_charge"].append(credit)
            for subset, values in credits.items():
                ref = required[target, cap, arm, subset]
                assert int(ref["pairs"]) == len(values) and float(
                    ref["correct"]
                ) == sum(values)
            nmatch = len(credits["same_size_charge"])
            actual = required[target, cap, arm, "same_size_charge"]
            assert nmatch == int(actual["pairs"])
            matched_counts[target] = nmatch
            all_credits.setdefault(target + "_" + cap, {})[arm] = sum(credits["all"])
        name = max(poses, key=lambda k: poses[k][0])
        score, pose, rec = poses[name]
        distances = cdist(
            xyz(pose), xyz((work / "prepared" / (rec + ".pdbqt")).read_text())
        ).min(axis=1)
        report["seeds"][str(seed)] = {
            "runs": len(raw_runs),
            "all_credits": all_credits,
            "matched_pairs": matched_counts,
            "maximum_first_score_key": name,
            "maximum_first_score": float(score),
            "nearest_heavy_A": round(float(distances.min()), 5),
            "ligand_heavy_atoms_within_2A": int((distances < 2).sum()),
        }
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    args.out.write_text(json.dumps(main(args.work), indent=2) + "\n")

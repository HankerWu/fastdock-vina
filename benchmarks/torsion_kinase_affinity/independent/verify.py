"""Independent confirmation of Davis Table 4 matched Kd and raw Vina pair signs."""

import argparse
import csv
import itertools
import json
import re
import statistics
from pathlib import Path

import numpy as np
import xlrd
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize


def audit(work, raw, dest):
    dest.mkdir(exist_ok=True, parents=True)
    expected = {"AURKA": ("2X6E", "2J4Z"), "CHEK1": ("2E9P", "2YWP")}
    arms = ("native", "balanced", "prepared")
    seeds = (60913, 60943, 61211, 61231)
    panel = list(csv.DictReader((work / "panel.csv").open()))
    t1, t4 = [
        xlrd.open_workbook(str(raw / f)).sheet_by_index(0)
        for f in ("table1.xls", "table4.xls")
    ]
    col = t4.row_values(0)
    secondary = json.loads((raw / "ligands_iso.txt").read_text())
    proteins = json.loads((raw / "proteins.txt").read_text())
    ext = np.loadtxt(
        raw / "drug-target_interaction_affinities_Kd__Davis_et_al.2011v1.txt"
    )
    check = []
    sign_map = {}
    details = {}
    for tgt, states in expected.items():
        (labels,) = [r for r in t4.get_rows() if r[2].value == tgt]
        (assay,) = [r for r in t1.get_rows() if r[2].value == tgt]
        assert labels[0].value == assay[0].value and assay[3].value == "NO"
        for r in [v for v in panel if v["target"] == tgt]:
            cid = r["cid"]
            kd = float(labels[col.index(r["compound"])].value)
            key = tgt, cid
            assert (
                key not in sign_map
                and kd == float(r["kd_nm"])
                and labels[0].value == r["accession"]
            )
            sign_map[key] = kd
            original = Chem.MolFromSmiles(r["source_smiles"])
            assert Chem.MolToInchiKey(original) == r["inchi_key"]
            parent = Chem.MolFromSmiles(r["smiles"])
            ligands = [
                Chem.MolToSmiles(x) for x in Chem.GetMolFrags(original, asMols=True)
            ]
            assert Chem.MolToSmiles(parent) in ligands
            pdbqt = (work / "prepared" / f"{cid}.pdbqt").read_text()
            (prepared,) = re.findall(
                r"^REMARK SMILES (?!IDX)(\S+)", pdbqt, re.MULTILINE
            )
            saved = Chem.MolFromSmiles(prepared)
            normalized = lambda m: Chem.MolToInchiKey(rdMolStandardize.ChargeParent(m))
            assert normalized(saved).split("-")[0] == normalized(parent).split("-")[0]
            unspecified = [
                str(v.centeredOn)
                for v in Chem.FindPotentialStereo(parent)
                if str(v.specified) == "Unspecified"
            ]
            if not unspecified:
                assert normalized(saved) == normalized(parent), (
                    cid,
                    prepared,
                    r["smiles"],
                )
            external = None
            if cid in secondary:
                assert normalized(Chem.MolFromSmiles(secondary[cid])) == normalized(
                    parent
                )
                external = ext[list(secondary).index(cid), list(proteins).index(tgt)]
                assert external == kd, (cid, external, kd)
            check.append(
                {
                    "target": tgt,
                    "compound": r["compound"],
                    "cid": cid,
                    "source_Kd_nM": kd,
                    "DeepDTA_CID_matched": external is not None,
                    "formal_charge": Chem.GetFormalCharge(parent),
                    "graph_unspecified_stereo": ",".join(unspecified),
                    "input_parent_key": normalized(parent),
                    "prepared_pose_key": normalized(saved),
                }
            )
            for st, sd, arm in itertools.product(states, seeds, arms):
                stem = work / "pools" / f"{st}_{cid}_{sd}_{arm}"
                text = stem.with_suffix(".pdbqt").read_text()
                assert text.endswith("ENDMDL\n")
                log = stem.with_suffix(".log").read_text()
                saved = float(
                    re.search(
                        r"^REMARK VINA RESULT:\s+(\S+)", text, re.MULTILINE
                    ).group(1)
                )
                rank = float(
                    re.search(r"^\s*1\s+(-?\d+(?:\.\d+)?)\s+", log, re.MULTILINE).group(
                        1
                    )
                )
                assert abs(saved - rank) <= 0.00501
                details[tgt, cid, st, sd, arm] = (saved, rank)
    assert len(details) == 816 and len(sign_map) == 34
    results = []
    for tgt, states in expected.items():
        labels = {cid: kd for (t, cid), kd in sign_map.items() if t == tgt}
        for scope in ["full", "stereodefined"] if tgt == "AURKA" else ["full"]:
            ids = [c for c in labels if scope != "stereodefined" or c != "9869779"]
            edges = [
                (i, j)
                for i, j in itertools.permutations(ids, 2)
                if labels[j] >= 3 * labels[i]
            ]
            for field, fid in (("serialized", 0), ("stdout", 1)):
                for cohort, sel in (
                    ("all", seeds),
                    ("first", seeds[:2]),
                    ("second", seeds[2:]),
                ):
                    outcomes = {}
                    for arm in arms:
                        energies = {
                            i: statistics.mean(
                                statistics.median(
                                    details[tgt, i, s, z, arm][fid] for z in sel
                                )
                                for s in states
                            )
                            for i in ids
                        }
                        outcomes[arm] = np.array(
                            [
                                1
                                if energies[i] + 1e-9 < energies[j]
                                else 0
                                if energies[j] + 1e-9 < energies[i]
                                else 0.5
                                for i, j in edges
                            ]
                        )
                    for arm in arms:
                        d = outcomes[arm] - outcomes["native"]
                        results.append(
                            {
                                "target": tgt,
                                "scope": scope,
                                "energy_source": field,
                                "cohort": cohort,
                                "arm": arm,
                                "compounds": len(ids),
                                "pairs": len(edges),
                                "correct": sum(outcomes[arm]),
                                "fixed": int(sum(d > 0)),
                                "lost": int(sum(d < 0)),
                                "net": float(d.sum()),
                            }
                        )
    for destfile, rows in (
        ("independent_ranking.csv", results),
        ("independent_chemistry.csv", check),
    ):
        with (dest / destfile).open("w") as out:
            w = csv.DictWriter(out, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("work", "raw", "out"):
        p.add_argument(key, type=Path)
    a = p.parse_args()
    audit(a.work, a.raw, a.out)

"""Independent reconstruction of two paired-crossdocking budgets from native PDB and PDBQT."""

import argparse
import csv
import itertools
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial.transform import Rotation

COMP = {
    "1A52": "EST",
    "3ERD": "DES",
    "2PRG": "BRL",
    "5Y2O": "8N6",
    "1XOQ": "ROF",
    "1XOM": "CIO",
}
TARGET = {
    s: g
    for g, q in [
        ("ESR1", ("1A52", "3ERD")),
        ("PPARG", ("2PRG", "5Y2O")),
        ("PDE4D", ("1XOQ", "1XOM")),
    ]
    for s in q
}


def xyz(r):
    return np.array([float(r[i : i + 8]) for i in (30, 38, 46)])


def complex_data(root, p):
    records = (root / "raw" / (p + ".pdb")).read_text().splitlines()
    pocket_ca = {
        (r[17:20].strip(), r[22:26], r[26]): xyz(r)
        for r in records
        if r.startswith("ATOM")
        and r[21] == "A"
        and r[16] in " A"
        and r[12:16].strip() == "CA"
    }
    dep = [
        r
        for r in records
        if r.startswith("HETATM")
        and r[21] == "A"
        and r[17:20].strip() == COMP[p]
        and r[76:78].strip() not in ("H", "D")
    ]
    alts = sorted({r[16] for r in dep if r[16] != " "}) or [" "]
    topology = (root / p / "ligand.pdbqt").read_text().splitlines()
    # Use canonical SMILES IDX and atom typing from molecular similarity, not original input-order mapping.
    canonical = Chem.MolFromSmiles(
        next(
            l.split(maxsplit=2)[2]
            for l in topology
            if l.startswith("REMARK SMILES ") and not l.startswith("REMARK SMILES IDX")
        )
    )
    serials = {}
    for l in topology:
        if l.startswith("REMARK SMILES IDX"):
            v = list(map(int, l.split()[3:]))
            serials.update(zip(v[::2], v[1::2]))
    serials = [serials[i + 1] for i in range(canonical.GetNumAtoms())]
    prep = {int(r[6:11]): xyz(r) for r in topology if r.startswith(("ATOM", "HETATM"))}
    template = Chem.Mol(canonical)
    conformer = Chem.Conformer(len(serials))
    for idx, serial in enumerate(serials):
        conformer.SetAtomPosition(idx, prep[serial])
    template.AddConformer(conformer)
    Chem.AssignStereochemistryFrom3D(template)
    assert Chem.MolToSmiles(template) == Chem.MolToSmiles(canonical)
    crystals = []
    expanded = []
    center = None
    for alt in alts:
        rl = [r[:16] + " " + r[17:] for r in dep if r[16] in (" ", alt)]
        if center is None:
            center = np.mean([xyz(r) for r in rl], axis=0)
        bound = Chem.MolFromPDBBlock(
            "\n".join(rl) + "\nEND\n",
            sanitize=False,
            removeHs=False,
            proximityBonding=True,
        )
        # An isolated neutral acid has two element-equivalent terminal oxygens before bond
        # assignment. Preserve CIO's CCD O25-C23 single / O24=C23 double bond explicitly.
        bound = AllChem.AssignBondOrdersFromTemplate(canonical, bound)
        if COMP[p] == "CIO":
            rw = Chem.RWMol(bound)
            ids = {
                a.GetPDBResidueInfo().GetName().strip(): a.GetIdx()
                for a in rw.GetAtoms()
            }
            for oxygen, order, nh in [
                ("O24", Chem.BondType.DOUBLE, 0),
                ("O25", Chem.BondType.SINGLE, 1),
            ]:
                rw.GetBondBetweenAtoms(ids["C23"], ids[oxygen]).SetBondType(order)
                atom = rw.GetAtomWithIdx(ids[oxygen])
                atom.SetNumExplicitHs(nh)
                atom.SetNoImplicit(True)
            bound = rw.GetMol()
            Chem.SanitizeMol(bound)
        Chem.AssignStereochemistryFrom3D(bound)
        matches = bound.GetSubstructMatches(
            template, uniquify=False, useChirality=True, maxMatches=10000
        )
        if not matches:
            raise ValueError((p, alt, "bond/stereochemistry mismatch"))
        altmaps = bound.GetSubstructMatches(
            canonical, uniquify=False, useChirality=True
        )
        for m in matches:
            crystals.append(bound.GetConformer().GetPositions()[list(m)])
        for m in set(matches) | set(altmaps):
            expanded.append(bound.GetConformer().GetPositions()[list(m)])
    return pocket_ca, center, crystals, serials, expanded


def fitted(a, b):
    ca1, c1 = a[:2]
    ca2, c2 = b[:2]
    ids = sorted(
        k
        for k in ca1.keys() & ca2.keys()
        if math.dist(ca1[k], c1) < 14 and math.dist(ca2[k], c2) < 14
    )
    x, y = [np.stack([mapping[k] for k in ids]) for mapping in (ca1, ca2)]
    rotation, _ = Rotation.align_vectors(y - y.mean(0), x - x.mean(0))

    def apply(coords):
        return rotation.apply(coords - x.mean(0)) + y.mean(0)

    return (
        np.array([apply(c) for c in a[2]]),
        np.array([apply(c) for c in a[4]]),
        len(ids),
        float(np.linalg.norm(apply(x) - y) / math.sqrt(len(ids))),
    )


def main(root, output):
    data = {p: complex_data(root, p) for p in COMP}
    rows = []
    max_error = max_align = 0.0
    model_counts = {}
    cost = []
    parity = defaultdict(dict)
    for budget in (3000, 12000):
        directory = root / f"runs_{budget}"
        with (root / f"poses_{budget}.csv").open() as h:
            reported = {
                (r["ligand"], r["receiver"], r["seed"], r["arm"], r["rank"]): r
                for r in csv.DictReader(h)
            }
        count = 0
        for path in sorted(directory.glob("*.pdbqt")):
            s, r, seed, arm = path.stem.split("_")
            refs, expanded, n, backbone = fitted(data[s], data[r])
            distances = []
            symdist = []
            log = path.with_suffix(".log").read_text()
            trajectories = [
                {k: int(v) for k, v in re.findall(r"(\w+)=(\d+)", l)}
                for l in log.splitlines()
                if l.startswith("Search trajectory:")
            ]
            assert [t["index"] for t in trajectories] == [1, 2, 3, 4]
            parity[budget, s, r, seed][arm] = [
                [t["hunt_evals"], t["second_evals"], t["steps"]] for t in trajectories
            ]
            final = {
                k: int(v)
                for k, v in re.findall(
                    r"(\w+)=(\d+)",
                    next(
                        l for l in log.splitlines() if l.startswith("Search counters:")
                    ),
                )
            }
            totals = {
                k: sum(t[k] for t in trajectories)
                for k in ("hunt_evals", "second_evals", "steps")
            }
            costs = {
                "search_evals": totals["hunt_evals"] + totals["second_evals"],
                "explicit_evals": final["explicit_evals"],
                "explicit_calls": final["explicit_calls"],
                "explicit_inputs": final["explicit_inputs"],
                "late_prescore": final["late_prescore"],
                "steps": totals["steps"],
            }
            cost.append(
                dict(budget=budget, source=s, receiver=r, seed=seed, arm=arm, **costs)
            )
            for model in path.read_text().split("ENDMDL"):
                if "MODEL" not in model:
                    continue
                rlines = model.splitlines()
                coords_by_serial = {
                    int(l[6:11]): xyz(l)
                    for l in rlines
                    if l.startswith(("ATOM", "HETATM"))
                }
                v = np.stack([coords_by_serial[k] for k in data[s][3]])
                rms = min(
                    float(np.linalg.norm(v - rr) / math.sqrt(len(v))) for rr in refs
                )
                distances.append(rms)
                symdist.append(
                    min(
                        float(np.linalg.norm(v - rr) / math.sqrt(len(v)))
                        for rr in expanded
                    )
                )
                count += 1
                saved = reported[s, r, seed, arm, str(len(distances))]
                assert n == int(saved["n_ca"])
                diff = abs(rms - float(saved["rmsd"]))
                max_error = max(max_error, diff)
                max_align = max(max_align, abs(backbone - float(saved["ca_rmsd"])))
            rows.append(
                {
                    "budget": budget,
                    "target": TARGET[s],
                    "source": s,
                    "receiver": r,
                    "kind": "cognate" if s == r else "cross",
                    "seed": seed,
                    "arm": arm,
                    "n_ca": n,
                    "ca_rmsd": backbone,
                    "top": distances[0],
                    "best": min(distances),
                    "top_hit": int(distances[0] <= 2),
                    "pool_hit": int(min(distances) <= 2),
                    "relative_stereo_top": symdist[0],
                    "relative_stereo_best": min(symdist),
                    "count": len(distances),
                }
            )
        assert count == len(reported)
        model_counts[budget] = count
    n_odd = 0
    n_even = 0
    for data_by_arm in parity.values():
        if set(data_by_arm) != {"native", "soft", "mixed"}:
            raise ValueError("Missing arm in paired evaluation")
        native = data_by_arm["native"]
        soft = data_by_arm["soft"]
        assert [native[0], native[2]] == [soft[0], soft[2]]
        n_odd += 2
        n_even += sum(native[k] != soft[k] for k in (1, 3))
    summary = []
    for budget, arm, kind in itertools.product(
        (3000, 12000), ("native", "soft", "mixed"), ("cognate", "cross")
    ):
        sel = [
            r for r in rows if (r["budget"], r["arm"], r["kind"]) == (budget, arm, kind)
        ]
        summary.append(
            {
                "budget": budget,
                "arm": arm,
                "kind": kind,
                "cells": len(sel),
                "top": sum(r["top_hit"] for r in sel),
                "pool": sum(r["pool_hit"] for r in sel),
                "relative_stereo_top": sum(r["relative_stereo_top"] <= 2 for r in sel),
                "relative_stereo_pool": sum(
                    r["relative_stereo_best"] <= 2 for r in sel
                ),
            }
        )
    cost_summary = []
    for budget, arm in itertools.product((3000, 12000), ("native", "soft", "mixed")):
        chosen = [r for r in cost if r["budget"] == budget and r["arm"] == arm]
        cost_summary.append(
            {
                "budget": budget,
                "arm": arm,
                "search_evals": sum(r["search_evals"] for r in chosen),
                "explicit_evals": sum(r["explicit_evals"] for r in chosen),
                "late_prescore": sum(r["late_prescore"] for r in chosen),
                "explicit_calls": sum(r["explicit_calls"] for r in chosen),
            }
        )
    assert max_error < 1e-7 and max_align < 1e-7

    def write(name, items):
        with (output / name).open("w") as f:
            w = csv.DictWriter(f, fieldnames=items[0])
            w.writeheader()
            w.writerows(items)

    write("review_cells.csv", rows)
    write("review_cost_by_run.csv", cost)
    out = {
        "pools": len(rows),
        "model_counts": model_counts,
        "largest_RMSD_difference_A": max_error,
        "largest_CA_alignment_difference_A": max_align,
        "identical_odd_trajectory_counters": n_odd,
        "changed_even_trajectory_counters": n_even,
        "summary": summary,
        "costs": cost_summary,
    }
    (output / "cross_audit.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("root", type=Path)
    ap.add_argument("output", type=Path)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    main(args.root, args.output)

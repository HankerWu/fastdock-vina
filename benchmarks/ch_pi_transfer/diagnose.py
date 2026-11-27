"""Independent carbon typing, atom-pair sum and localization on both lectins."""

import argparse
import json
import sys
from functools import cache
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from scipy.spatial import cKDTree
from scipy.stats import kendalltau, spearmanr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "galectin_affinity"))
from sugar import ligand_core

from model import CHPi, kernel

# Explicit conventional protein carbon classes independently enumerated by name.
RINGS = {
    "PHE": "CG CD1 CD2 CE1 CE2 CZ",
    "TYR": "CG CD1 CD2 CE1 CE2 CZ",
    "TRP": "CG CD1 CD2 CE2 CE3 CZ2 CZ3 CH2",
    "HIS": "CG CD2 CE1",
}
SIDE = {
    "ALA": "CB",
    "ARG": "CB CG CD",
    "ASN": "CB",
    "ASP": "CB",
    "CYS": "CB",
    "GLN": "CB CG",
    "GLU": "CB CG",
    "GLY": "",
    "HIS": "CB",
    "ILE": "CB CG1 CG2 CD1",
    "LEU": "CB CG CD1 CD2",
    "LYS": "CB CG CD CE",
    "MET": "CB CG CE",
    "PHE": "CB",
    "PRO": "CB CG CD",
    "SER": "CB",
    "THR": "CB CG2",
    "TRP": "CB",
    "TYR": "CB",
    "VAL": "CB CG1 CG2",
}
RENAMED = {
    "HID": "HIS",
    "HIE": "HIS",
    "HIP": "HIS",
    "HSD": "HIS",
    "HSE": "HIS",
    "HSP": "HIS",
}
FIELDS = ["ringCH_Trp", "ringCH_other", "otherCH_Trp", "otherCH_other", "reverse"]


def pos(l):
    return [float(l[i : i + 8]) for i in (30, 38, 46)]


@cache
def neighbors(path):
    points, typ, trp = [], [], []
    for l in path.read_text().splitlines():
        if not l.startswith(("ATOM  ", "HETATM")):
            continue
        res, atom = l[17:20], l[12:16].strip()
        res = RENAMED.get(res, res)
        if atom in RINGS.get(res, "").split():
            kind = "arom"
        elif res in SIDE and (atom == "CA" or atom in SIDE[res].split()):
            kind = "sp3"
        else:
            continue
        points.append(pos(l))
        typ.append(kind)
        trp.append(res == "TRP")
    return cKDTree(points), np.asarray(typ), np.asarray(trp)


@cache
def identities(path):
    lines = path.read_text().splitlines()
    (smi,) = [
        s.split()[2] for s in lines if s.startswith("REMARK SMILES ") and "IDX" not in s
    ]
    mol = Chem.MolFromSmiles(smi)
    idx = []
    for l in lines:
        if l.startswith("REMARK SMILES IDX "):
            vals = list(map(int, l.split()[3:]))
            idx.extend(zip(vals[::2], vals[1::2]))
    s = dict(idx)

    def matches(q):
        return [s[i[0] + 1] for i in mol.GetSubstructMatches(Chem.MolFromSmarts(q))]

    return matches("[C;X4;!H0]"), matches("[c]")


def decomposition(rpath, lpath, block, core):
    tree, kinds, trp = neighbors(rpath)
    ch, arom = identities(lpath)
    coords = {
        int(l[6:11]): pos(l)
        for l in block.splitlines()
        if l.startswith(("ATOM  ", "HETATM"))
    }
    values = np.zeros(5)
    for ids, partner_kind in [(ch, "arom"), (arom, "sp3")]:
        for i in ids:
            x = np.array(coords[i])
            nb = np.array(tree.query_ball_point(x, 8.0), dtype=int)
            nb = nb[kinds[nb] == partner_kind]
            dist = np.linalg.norm(x - tree.data[nb], axis=1)
            w = -0.26 * np.exp(-((dist - 4.49) ** 2) / (2 * 0.75**2))
            w[dist >= 8] = 0.0
            if partner_kind == "arom":
                k = 0 if i in core else 2
                values[k] += w[trp[nb]].sum()
                values[k + 1] += w[~trp[nb]].sum()
            else:
                values[-1] += w.sum()
    return values


def analytic_control():
    # One methane carbon at the origin, one phenyl carbon at the radial optimum.
    def line(i, name, res, typ, x):
        return (
            f"ATOM  {i:5d} {name:4s} {res:3s} A   1    "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00     0.000 {typ}"
        )

    template = "REMARK SMILES C\nREMARK SMILES IDX 1 1\n" + line(1, "C", "UNL", "C", 0)
    rec = line(1, "CZ", "PHE", "A", 4.49)
    expected = (-0.26, 0.0)
    if not np.allclose(CHPi(rec, template).components(template), expected, atol=1e-14):
        raise ValueError("one pair at optimal separation")
    mirror = line(1, "CZ", "PHE", "A", 20 - 4.49)
    ligand = template.replace(line(1, "C", "UNL", "C", 0), line(1, "C", "UNL", "C", 20))
    if not np.allclose(CHPi(mirror, ligand).components(ligand), expected, atol=1e-14):
        raise ValueError("rigid-coordinate invariance")
    if tuple(CHPi(line(1, "CZ", "ARG", "C", 4.49), template).components(template)) != (
        0.0,
        0.0,
    ):
        raise ValueError(
            "arginine guanidinium must not act as a CH donor or an aromatic carbon"
        )
    if (
        kernel(8.0) != 0.0
        or abs(kernel(4.49 + 0.75) / kernel(4.49) - np.exp(-0.5)) > 1e-14
    ):
        raise ValueError("radial width or boundary condition")


def analyze(gal, out):
    analytic_control()
    poses = pd.read_csv(out / "poses.csv")
    detail = []
    for (p, compound, seed), group in poses[poses.target.str.startswith("Gal")].groupby(
        ["pdb", "compound", "seed"]
    ):
        rpath = gal / "prepared" / f"{p}.pdbqt"
        lpath = gal / "prepared" / f"{compound}.pdbqt"
        core = set(ligand_core(gal, int(compound)))
        pool = (
            (gal / "pools" / f"{p}_{compound}_{seed}_native.pdbqt")
            .read_text()
            .split("ENDMDL")[:-1]
        )
        for i, (b, r) in enumerate(zip(pool, group.itertuples(), strict=True), 1):
            comp = decomposition(rpath, lpath, b, core)
            if (
                abs(sum(comp[:4]) - r.forward_raw) > 1e-7
                or abs(comp[-1] - r.reverse_raw) > 1e-7
            ):
                raise ValueError((p, compound, seed, i, comp, r))
            detail.append(
                dict(
                    target=r.target,
                    pdb=p,
                    compound=compound,
                    seed=seed,
                    mode=i,
                    native=r.native,
                    rmsd=r.rmsd,
                    factor=r.factor,
                    **dict(zip(FIELDS, comp)),
                )
            )
    raw = pd.DataFrame(detail)
    raw.to_csv(out / "localization.csv", index=False, float_format="%.9f")
    contrasts = []
    # Ablating groups is post hoc. Never use ligand or receptor identities in prediction.
    for key, cell in raw.groupby(["target", "pdb", "compound", "seed"]):
        for name, terms in [
            ("only_ring_CH_Trp", FIELDS[:1]),
            ("ring_CH_all", FIELDS[:2]),
            ("all_ligand_CH_except_ring_Trp", FIELDS[1:4]),
            ("all_ligand_CH", FIELDS[:4]),
            ("all_terms", FIELDS),
        ]:
            s = cell.native + 0.3 * cell.factor * cell[terms].sum(axis=1)
            winner = cell.loc[s.idxmin()]
            first = cell.iloc[0]
            contrasts.append(
                dict(
                    zip(["target", "pdb", "compound", "seed"], key),
                    ablation=name,
                    winner=winner["mode"],
                    native=first.rmsd,
                    selected=winner.rmsd,
                    delta=winner.rmsd - first.rmsd,
                    native_score_gap=winner.native - first.native,
                    **{k: 0.3 * winner.factor * (winner[k] - first[k]) for k in FIELDS},
                )
            )
    c = pd.DataFrame(contrasts)
    c.to_csv(out / "local_cells.csv", index=False, float_format="%.9f")
    s = c.groupby(["target", "ablation"]).agg(
        cells=("delta", "size"),
        changed=("winner", lambda x: (x != 1).sum()),
        better=("delta", lambda x: (x < 0).sum()),
        worse=("delta", lambda x: (x > 0).sum()),
        mean=("selected", "mean"),
        median=("selected", "median"),
        mean_delta=("delta", "mean"),
    )
    s.to_csv(out / "local_summary.csv", float_format="%.9f")

    lig = pd.read_csv(out / "ligand_scores.csv")
    ref = pd.read_csv(out / "affinity.csv")
    for (target, pdb), g in lig.groupby(["target", "pdb"]):
        for r in ref[(ref.target == target) & (ref.pdb == pdb)].itertuples():
            scores = g[r.arm]
            rho = spearmanr(-scores, g.pKd).statistic
            tau = kendalltau(-scores, g.pKd).statistic
            if max(abs(rho - r.rho), abs(tau - r.tau)) > 1e-7:
                raise ValueError("independent SciPy rank mismatch")
    # Table 1 compounds 5–16 share the aryltriazole series; 17 lacks this extension.
    descriptors = {}
    for compound in range(5, 18):
        lines = (gal / "prepared" / f"{compound}.pdbqt").read_text().splitlines()
        (smiles,) = [
            l.split()[2]
            for l in lines
            if l.startswith("REMARK SMILES ") and " IDX " not in l
        ]
        m = Chem.MolFromSmiles(smiles)
        descriptors[compound] = (
            m.GetNumHeavyAtoms(),
            Chem.GetFormalCharge(m),
            rdMolDescriptors.CalcMolFormula(m),
        )
    rows = []
    pairs = pd.read_csv(out / "pairs.csv")
    for (target, pdb), g in pairs[pairs.resolved].groupby(["target", "pdb"]):
        for subset in ("all", "triazoles", "same_size_and_charge", "same_formula"):
            for arm in (
                "native",
                "ligand_CH",
                "protein_CH",
                "bidirectional",
                "size_control",
            ):
                votes = []
                for r in g.itertuples():
                    a, b = descriptors[r.lig1], descriptors[r.lig2]
                    if (
                        (subset == "triazoles" and 17 in (r.lig1, r.lig2))
                        or (subset == "same_size_and_charge" and a[:2] != b[:2])
                        or (subset == "same_formula" and a[2] != b[2])
                    ):
                        continue
                    gap = (
                        a[0] - b[0]
                        if arm == "size_control"
                        else getattr(r, arm + "_weak_minus_strong")
                    )
                    votes.append((1 + np.sign(gap)) / 2)
                rows.append(
                    {
                        "target": target,
                        "pdb": pdb,
                        "subset": subset,
                        "arm": arm,
                        "pairs": len(votes),
                        "credit": sum(votes),
                    }
                )
    pd.DataFrame(rows).to_csv(out / "size_control.csv", index=False)
    print(
        json.dumps(
            {
                "galectin_modes_checked": len(raw),
                "affinity_cells_checked": len(ref),
                "kcal_tolerance": 1e-7,
            }
        )
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("gal", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    analyze(a.gal, a.out)

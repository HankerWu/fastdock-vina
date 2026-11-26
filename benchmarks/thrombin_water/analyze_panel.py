"""Score the frozen ensembles with and without the crystallographically observed S1 water prior."""

import csv
import itertools
import re
from pathlib import Path

import numpy as np
from rdkit import Chem

from run_panel import CORE, REC, SEEDS, SERIES, ca, lines, transform, xyz

RT = 0.00831446261815324 * 298.15


def write(name, rows):
    if not rows:
        return
    with Path(name).open("w") as h:
        w = csv.DictWriter(h, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)


def element_graph(mol):
    graph = Chem.RWMol()
    for a in mol.GetAtoms():
        graph.AddAtom(Chem.Atom(a.GetAtomicNum()))
    for b in mol.GetBonds():
        graph.AddBond(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), Chem.BondType.SINGLE)
    return graph.GetMol()


def crystal_refs(pdb, ccd, ring, rec):
    pdb_lines = [r for r in lines(pdb) if r[17:20] == ccd]
    ref = Chem.MolFromPDBBlock(
        "\n".join(pdb_lines) + "\nEND\n",
        sanitize=False,
        removeHs=False,
        proximityBonding=True,
    )
    refg = element_graph(ref)
    m = Chem.MolFromSmiles(CORE + ring)
    matches = refg.GetSubstructMatches(element_graph(m), uniquify=False)
    if (
        not matches
        or ref.GetNumBonds() != m.GetNumBonds()
        or ref.GetNumAtoms() != m.GetNumAtoms()
    ):
        raise ValueError((pdb, len(matches)))
    rot, t = transform(pdb, rec)
    return np.array(
        [ref.GetConformer().GetPositions()[list(map_)] @ rot + t for map_ in matches]
    )


def protein_site():
    # Reference water is chosen solely from independent 3D49 and two protein carbonyls
    # cited in Baum's S1 mechanism; other ligand positions never enter the scoring rule.
    src = lines("3D49")
    oxy = [
        xyz(
            next(
                r
                for r in src
                if r[21] == "H"
                and r[22:26].strip() == str(x)
                and r[12:16].strip() == "O"
            )
        )
        for x in (215, 227)
    ]
    site = min(
        (r for r in src if r[17:20] == "HOH"),
        key=lambda r: np.linalg.norm(xyz(r) - np.mean(oxy, axis=0)),
    )
    return xyz(site), site


def steric(d):
    t = np.clip((d - 2.0) / 0.6, 0, 1)
    a = t * t * (3 - 2 * t)
    # One crystal votes for the site out of one, using the prior workflow's (k+.5)/(n+1).
    return -RT / 4.184 * np.log(0.25 + 0.75 * a)


def main():
    q, water = protein_site()
    geometry = []
    poses = []
    cells = []
    pairs = []
    for i, (ring, dg, p, c) in SERIES.items():
        m = Chem.MolFromSmiles(CORE + ring)
        n = m.GetNumAtoms()
        top = Path("prepared", i + ".pdbqt").read_text()
        inds = {}
        for l in top.splitlines():
            if l.startswith("REMARK INDEX MAP"):
                ints = [int(x) for x in l.split()[3:]]
                inds.update({a - 1: b for a, b in zip(ints[::2], ints[1::2]) if a <= n})
        assert set(inds) == set(range(n))
        if p:
            r, t = transform("3D49", p)
            s = q @ r + t
            f = lines(p)
            ll = [l for l in f if l[17:20] == c]
            nearW = min(
                (l for l in f if l[17:20] == "HOH"),
                key=lambda l: np.linalg.norm(xyz(l) - s),
            )
            nearL = min(ll, key=lambda l: np.linalg.norm(xyz(l) - s))
            geometry.append(
                {
                    "ligand": i,
                    "PDB": p,
                    "source_water": water[17:26].strip(),
                    "pocket_CA_RMSD": float(
                        np.sqrt(
                            np.mean(np.sum((ca("3D49") @ r + t - ca(p)) ** 2, axis=-1))
                        )
                    ),
                    "water_to_nearest_ligand": float(np.linalg.norm(xyz(nearL) - s)),
                    "nearest_ligand_atom": nearL[12:16].strip(),
                    "water_to_nearest_crystal_water": float(
                        np.linalg.norm(xyz(nearW) - s)
                    ),
                    "nearest_crystal_water": nearW[21:26].strip(),
                }
            )
        for rec in REC:
            rot, t = transform("3D49", rec)
            site = q @ rot + t
            refs = crystal_refs(p, c, ring, rec) if p else None
            for seed in SEEDS:
                path = Path("results", f"{rec}_{i}_{seed}.pdbqt")
                if not path.exists():
                    raise FileNotFoundError(path)
                block = path.read_text()
                if "ENDMDL" not in block:
                    raise ValueError(f"Incomplete pool: {path}")
                pool = []
                for model, body in enumerate(
                    re.findall(r"MODEL\s+\d+\n(.*?)ENDMDL", block, flags=re.DOTALL), 1
                ):
                    raw_score = float(
                        re.search(r"REMARK VINA RESULT:\s+([-\d.]+)", body)[1]
                    )
                    coords = {
                        int(l[6:11]): xyz(l)
                        for l in body.splitlines()
                        if l.startswith(("ATOM", "HETATM"))
                    }
                    pos = np.array([coords[inds[j]] for j in range(n)])
                    rmsd = (
                        np.sqrt(np.sum((refs - pos) ** 2, axis=-1).mean(axis=-1)).min()
                        if p
                        else None
                    )
                    dist = np.linalg.norm(pos - site, axis=-1).min()
                    row = {
                        "rec": rec,
                        "ligand": i,
                        "seed": seed,
                        "model": model,
                        "native": raw_score,
                        "exclusion": raw_score + steric(dist),
                        "distance": dist,
                        "rmsd": rmsd,
                    }
                    pool.append(row)
                    poses.append(row)
                for score in ("native", "exclusion"):
                    lead = min(pool, key=lambda row: (row[score], row["model"]))
                    cells.append(
                        {
                            "rec": rec,
                            "ligand": i,
                            "seed": seed,
                            "score": score,
                            "model": lead["model"],
                            "energy": lead[score],
                            "rmsd": lead["rmsd"],
                            "site_distance": lead["distance"],
                            "best_available_rmsd": min(r["rmsd"] for r in pool)
                            if p
                            else None,
                        }
                    )
    for score in ("native", "exclusion"):
        by = {
            (r["ligand"], r["rec"], r["seed"]): r["energy"]
            for r in cells
            if r["score"] == score
        }
        for a, b in itertools.combinations(SERIES, 2):
            gap = abs(SERIES[a][1] - SERIES[b][1])
            strong, weak = sorted((a, b), key=lambda i: SERIES[i][1])
            if gap < RT * np.log(3):
                continue
            comparable = [
                (r, s)
                for r, s in itertools.product(REC, SEEDS)
                if (strong, r, s) in by and (weak, r, s) in by
            ]
            if len(comparable) != len(REC) * len(SEEDS):
                continue
            for subset in ("all", *REC):
                obs = [
                    by[weak, r, s] - by[strong, r, s]
                    for r, s in comparable
                    if subset == "all" or subset == r
                ]
                pred = float(np.median(obs))
                pairs.append(
                    {
                        "score": score,
                        "subset": subset,
                        "strong": strong,
                        "weak": weak,
                        "experimental_gap_kcal": gap / 4.184,
                        "predicted_gap_kcal": pred,
                        "correct": 1 if pred > 1e-9 else 0.5 if abs(pred) < 1e-9 else 0,
                    }
                )
    write("geometry.csv", geometry)
    write("poses.csv", poses)
    write("cells.csv", cells)
    write("pairs.csv", pairs)
    print(
        "geometry",
        [
            (
                r["PDB"],
                round(r["water_to_nearest_ligand"], 2),
                round(r["water_to_nearest_crystal_water"], 2),
                round(r["pocket_CA_RMSD"], 2),
            )
            for r in geometry
        ],
    )
    for score in ("native", "exclusion"):
        for subset in ("all", "noncognate", *REC):
            c = [
                r
                for r in cells
                if r["score"] == score
                and r["rmsd"] is not None
                and (
                    subset == "all"
                    or subset == "noncognate"
                    and r["rec"] != SERIES[r["ligand"]][2]
                    or r["rec"] == subset
                )
            ]
            ps = [r for r in pairs if r["score"] == score and r["subset"] == subset]
            print(
                score,
                subset,
                "pose",
                sum(x["rmsd"] <= 2 for x in c),
                len(c),
                "oracle",
                sum(x["best_available_rmsd"] <= 2 for x in c),
                "pairs",
                sum(x["correct"] for x in ps),
                len(ps),
            )


if __name__ == "__main__":
    main()

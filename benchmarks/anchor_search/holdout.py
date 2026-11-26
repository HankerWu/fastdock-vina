"""Reproduce independent small-sulfonamide pose tests on human CA I and CA XII."""

import argparse
import csv
import hashlib
import itertools
import json
import os
import re
import subprocess
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import gemmi
import numpy as np
from meeko import MoleculePreparation, PDBQTWriterLegacy
from openmm import Platform
from openmm.app import PDBFile
from pdbfixer import PDBFixer
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
from ionization import anion_pair_offset
from select_zinc_pose import assess_site, parse_atoms, read_models, select_pose

SOURCES = {
    "1AZM.pdb": "cfe92653bb7633f9fb87fd92f0c4a8e8ad3f9405449c0dc8ca3191d6c2c9c6ea",
    "1BZM.pdb": "e2e8a855618565bce08033ec4a07d75b87868aaf524559ea3ab4d747984a03a9",
    "1JD0.pdb": "325ef3ab8a8080af58bf6d4eb0a9f89cf1a8e3f8639e4bb3c071f9c7bc7bf9dc",
    "AZM.cif": "9afddb10dce8ffc0a481c46ee659eae0adfc75e42eb92388a5d4165d49ec7445",
    "MZM.cif": "1b3a2e38b5de3301512d566c10d2ab0cc0fd674ceef6f2239f7f38dfb0bd9acb",
}
COMPLEX = {"1AZM": "AZM", "1BZM": "MZM", "1JD0": "AZM"}
CASES = [
    ("1AZM", "1AZM"),
    ("1AZM", "1BZM"),
    ("1BZM", "1AZM"),
    ("1BZM", "1BZM"),
    ("1JD0", "1JD0"),
]
SEEDS = (314159, 271828)


def coords(line):
    return [float(line[k : k + 8]) for k in (30, 38, 46)]


def fetch(raw):
    for name, expected in SOURCES.items():
        path = raw / name
        if not path.exists():
            url = (
                "https://files.rcsb.org/"
                + ("download/" if name.endswith(".pdb") else "ligands/download/")
                + name
            )
            with urllib.request.urlopen(url, timeout=60) as response:
                path.write_bytes(response.read())
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f"Source checksum differs: {path}")


def receptor(prep, pdb, lines, box_size):
    protein = [
        l
        for l in lines
        if l.startswith("ATOM  ")
        and l[21] == "A"
        and l[16] in " A"
        and l[76:78].strip() not in ("H", "D")
    ]
    (zn,) = [
        l
        for l in lines
        if l.startswith("HETATM") and l[21] == "A" and l[17:20] == " ZN"
    ]
    protein_file = prep / f"{pdb}_protein.pdb"
    protein_file.write_text(
        "\n".join(l[:16] + " " + l[17:] for l in protein) + "\nEND\n"
    )
    os.environ["OPENMM_CPU_THREADS"] = "1"
    fixer = PDBFixer(
        filename=str(protein_file), platform=Platform.getPlatformByName("CPU")
    )
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    missing = {
        f"{r.name}:{r.id}": [a.name for a in atoms]
        for r, atoms in fixer.missingAtoms.items()
    }
    fixer.addMissingAtoms(seed=2026)
    with protein_file.open("w") as fp:
        PDBFile.writeFile(fixer.topology, fixer.positions, fp, keepIds=True)
    out = prep / f"{pdb}_receptor.pdbqt"
    with (prep / f"{pdb}_receptor.log").open("w") as fp:
        subprocess.run(
            [
                sys.executable,
                "-m",
                "meeko.cli.mk_prepare_receptor",
                "--read_pdb",
                str(protein_file),
                "--set_template",
                "A:94,96=HID,A:119=HIE",
                "-p",
                str(out),
            ],
            check=True,
            stdout=fp,
            stderr=subprocess.STDOUT,
        )
    qt = out.read_text().rstrip()
    zserial = max(parse_atoms(qt.splitlines())) + 1
    qt += "\n" + zn[:6] + f"{zserial:5d}" + zn[11:54] + "  1.00  0.00    +2.000 Zn\n"
    out.write_text(qt)
    site = assess_site(parse_atoms(qt.splitlines()), zserial)
    if not site["applicable"]:
        raise ValueError(site)
    zn_xyz, ray = np.array(site["zinc_xyz"]), np.array(site["vacant_direction"])
    center = zn_xyz + 5.0 * ray
    (prep / f"{pdb}.box").write_text(
        "".join(f"center_{c} = {v:.5f}\n" for c, v in zip("xyz", center))
        + "".join(f"size_{c} = {box_size:d}\n" for c in "xyz")
    )
    return {
        "zinc": zn_xyz.tolist(),
        "target": (zn_xyz + 2 * ray).tolist(),
        "center": center.tolist(),
        "site": site,
        "missing_sidechain_atoms": missing,
    }


def ligand(raw, prep, code):
    block = gemmi.cif.read_file(str(raw / f"{code}.cif")).sole_block()
    descriptors = block.get_mmcif_category("_pdbx_chem_comp_descriptor.")
    (smiles,) = [
        s
        for t, p, s in zip(
            descriptors["type"], descriptors["program"], descriptors["descriptor"]
        )
        if t == "SMILES_CANONICAL" and p == "OpenEye OEToolkits"
    ]
    neutral = Chem.MolFromSmiles(smiles)
    atoms = block.get_mmcif_category("_chem_comp_atom.")
    bonds = block.get_mmcif_category("_chem_comp_bond.")
    graph = Chem.RWMol()
    names = []
    idx = {}
    for name, el, charge in zip(
        atoms["atom_id"], atoms["type_symbol"], atoms["charge"]
    ):
        if el == "H":
            continue
        atom = Chem.Atom(el.capitalize())
        atom.SetFormalCharge(int(charge))
        idx[name] = graph.AddAtom(atom)
        names.append(name)
    orders = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE}
    for a, b, order in zip(
        bonds["atom_id_1"], bonds["atom_id_2"], bonds["value_order"]
    ):
        if a in idx and b in idx:
            graph.AddBond(idx[a], idx[b], orders[order])
    graph = graph.GetMol()
    Chem.SanitizeMol(graph)
    if Chem.MolToSmiles(neutral, isomericSmiles=False) != Chem.MolToSmiles(
        graph, isomericSmiles=False
    ):
        raise ValueError("CCD bond graph disagrees with canonical SMILES")
    maps = [
        [names[j] for j in match]
        for match in graph.GetSubstructMatches(neutral, uniquify=False)
    ]
    mol = Chem.Mol(neutral)
    (match,) = mol.GetSubstructMatches(Chem.MolFromSmarts("[N;H2;+0]S(=O)=O"))
    donor = match[0]
    mol.GetAtomWithIdx(donor).SetFormalCharge(-1)
    Chem.SanitizeMol(mol)
    h = Chem.AddHs(mol)
    opts = AllChem.ETKDGv3()
    opts.randomSeed = 2026
    if (
        AllChem.EmbedMolecule(h, opts) != 0
        or AllChem.MMFFOptimizeMolecule(h, maxIters=1000) != 0
    ):
        raise ValueError(f"Failed conformer for {code}")
    (setup,) = MoleculePreparation().prepare(h)
    qt, ok, reason = PDBQTWriterLegacy.write_string(setup, add_index_map=True)
    if not ok:
        raise ValueError(reason)
    (prep / f"{code}.pdbqt").write_text(qt)
    imap = {}
    for line in qt.splitlines():
        if line.startswith("REMARK INDEX MAP "):
            vals = list(map(int, line.split()[3:]))
            imap.update(zip(vals[0::2], vals[1::2]))
    return {
        "heavy_serials": [imap[i + 1] for i in range(mol.GetNumAtoms())],
        "atom_names_by_map": maps,
        "donor_serial": imap[donor + 1],
        "smiles": Chem.MolToSmiles(mol),
    }


def prepare(work, box_size):
    prep = work / "prepared"
    raw = work / "raw"
    prep.mkdir(parents=True, exist_ok=True)
    raw.mkdir(exist_ok=True)
    fetch(raw)
    lines = {p: (raw / f"{p}.pdb").read_text().splitlines() for p in COMPLEX}
    ligs = {c: ligand(raw, prep, c) for c in sorted(set(COMPLEX.values()))}
    sites = {p: receptor(prep, p, ll, box_size) for p, ll in lines.items()}
    refs = {}
    for rec, src in CASES:
        lig = ligs[COMPLEX[src]]
        target = sites[rec]["target"]
        (prep / f"{rec}_{src}.anchor").write_text(
            f"{lig['donor_serial']} " + " ".join(f"{x:.8f}" for x in target) + "\n"
        )
        (prep / f"{rec}_{src}.shell").write_text(
            f"{lig['donor_serial']} "
            + " ".join(f"{x:.8f}" for x in sites[rec]["zinc"])
            + " shell 2\n"
        )
        # Parse crystallographic ligand coordinates exclusively for retrospective evaluation.
        points = {
            l[12:16].strip(): coords(l)
            for l in lines[src]
            if l.startswith("HETATM")
            and l[21] == "A"
            and l[17:20] == COMPLEX[src]
            and l[16] in " A"
            and l[76:78].strip() not in ("H", "D")
        }
        alpha = [
            {
                l[22:27]: coords(l)
                for l in lines[p]
                if l.startswith("ATOM  ")
                and l[21] == "A"
                and l[12:16].strip() == "CA"
                and l[16] in " A"
            }
            for p in (src, rec)
        ]
        keys = sorted(alpha[0].keys() & alpha[1].keys())
        x, y = [np.array([a[k] for k in keys]) for a in alpha]
        u, _, v = np.linalg.svd((x - x.mean(0)).T @ (y - y.mean(0)))
        rot = u @ np.diag([1, 1, np.linalg.det(u @ v)]) @ v
        ref = np.array(
            [[points[name] for name in names] for names in lig["atom_names_by_map"]]
        )
        if len(points) != len(lig["heavy_serials"]):
            raise ValueError(f"Incomplete crystal ligand in {src}")
        ref = (ref - x.mean(0)) @ rot + y.mean(0)
        margin = box_size / 2 - abs(ref - sites[rec]["center"]).max()
        if margin < 0:
            raise ValueError(f"Reference does not fit box for {rec}/{src}")
        refs[f"{rec}_{src}"] = {
            "coords": ref.tolist(),
            "ca_count": len(keys),
            "ca_rmsd": np.sqrt(
                np.mean(np.sum(((x - x.mean(0)) @ rot + y.mean(0) - y) ** 2, axis=1))
            ),
            "box_margin": float(margin),
        }
    (prep / "manifest.json").write_text(
        json.dumps({"sites": sites, "ligands": ligs, "refs": refs}, indent=2) + "\n"
    )


def dock(work, binary, arm, workers):
    p, out = work / "prepared", work / "runs"
    out.mkdir(exist_ok=True)

    def one(case):
        (rec, src), seed = case
        name = f"{rec}_{src}_{seed}_{arm}"
        cmd = [
            str(binary.resolve()),
            "--receptor",
            str(p / f"{rec}_receptor.pdbqt"),
            "--ligand",
            str(p / f"{COMPLEX[src]}.pdbqt"),
            "--config",
            str(p / f"{rec}.box"),
            "--pose_selection",
            "late",
            "--cpu",
            "1",
            "--seed",
            str(seed),
            "--exhaustiveness",
            "4",
            "--max_evals",
            "50000",
            "--num_modes",
            "9",
            "--energy_range",
            "100",
            "--search_counters",
            "--out",
            str(out / f"{name}.pdbqt"),
        ]
        if arm in ("anchor", "shell"):
            cmd.extend(["--search_anchor", str(p / f"{rec}_{src}.{arm}")])
        with (out / f"{name}.log").open("w") as log:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)

    with ThreadPoolExecutor(workers) as pool:
        list(pool.map(one, itertools.product(CASES, SEEDS)))


def save_rows(path, rows):
    with path.open("w") as output:
        writer = csv.DictWriter(output, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def evaluate(work):
    data = json.loads((work / "prepared/manifest.json").read_text())
    rows = []
    cells = []
    reranked = []
    receptors = {
        p: parse_atoms(
            (work / "prepared" / f"{p}_receptor.pdbqt").read_text().splitlines()
        )
        for p in COMPLEX
    }
    for (rec, src), seed, arm in itertools.product(
        CASES, SEEDS, ("native", "anchor", "shell")
    ):
        entry = data["ligands"][COMPLEX[src]]
        site = data["sites"][rec]
        ref = np.array(data["refs"][f"{rec}_{src}"]["coords"])
        name = f"{rec}_{src}_{seed}_{arm}"
        pdbqt = (work / "runs" / f"{name}.pdbqt").read_text()
        models = re.findall(
            r"^MODEL\s+\d+\n.*?^ENDMDL$", pdbqt, re.MULTILINE | re.DOTALL
        )
        if not 1 <= len(models) <= 9:
            raise ValueError(f"Missing or unexpected models for {name}")
        start = len(rows)
        for model in models:
            coords_by_id = parse_atoms(model.splitlines())
            xyz = np.array([coords_by_id[n].xyz for n in entry["heavy_serials"]])
            donor = np.array(coords_by_id[entry["donor_serial"]].xyz)
            rmsd = np.sqrt(np.mean(np.sum((ref - xyz) ** 2, axis=2), axis=1)).min()
            score = re.search(r"^REMARK VINA RESULT:\s+(\S+)", model, re.MULTILINE)[1]
            rows.append(
                {
                    "receptor": rec,
                    "ligand": COMPLEX[src],
                    "crystal": src,
                    "seed": seed,
                    "arm": arm,
                    "model": int(model.splitlines()[0].split()[-1]),
                    "score": float(score),
                    "rmsd": round(rmsd, 6),
                    "donor_zinc_A": round(np.linalg.norm(donor - site["zinc"]), 6),
                    "donor_target_error_A": round(
                        np.linalg.norm(donor - site["target"]), 6
                    ),
                }
            )
        cellrows = rows[start:]
        by_number = {r["model"]: r for r in cellrows}
        for rule in ("hard", "continuous", "radial"):
            result = select_pose(
                receptors[rec],
                read_models(pdbqt.splitlines()),
                site["site"]["zinc_serial"],
                entry["donor_serial"],
                rule,
            )
            selected = result["selected_model"]
            reranked.append(
                {
                    "receptor": rec,
                    "crystal": src,
                    "ligand": COMPLEX[src],
                    "seed": seed,
                    "arm": arm,
                    "rule": rule,
                    "status": result["status"],
                    "selected_model": selected,
                    "selected_rmsd": by_number[selected]["rmsd"]
                    if selected is not None
                    else None,
                }
            )
        first = cellrows[0]
        log = (work / "runs" / f"{name}.log").read_text()
        cells.append(
            {
                "receptor": rec,
                "crystal": src,
                "ligand": COMPLEX[src],
                "seed": seed,
                "arm": arm,
                "score": first["score"],
                "top1_rmsd": first["rmsd"],
                "best_rmsd": min(r["rmsd"] for r in cellrows),
                "n_models": len(cellrows),
                "hunt_evals": sum(map(int, re.findall(r"hunt_evals=(\d+)", log))),
                "second_evals": sum(map(int, re.findall(r"second_evals=(\d+)", log))),
            }
        )
    pairs = []
    for rec, seed, arm in itertools.product(
        ("1AZM", "1BZM"), SEEDS, ("native", "anchor", "shell")
    ):
        scores = {
            c["ligand"]: c["score"]
            for c in cells
            if (c["receptor"], c["seed"], c["arm"]) == (rec, seed, arm)
        }
        delta = scores["AZM"] - scores["MZM"]
        pairs.append(
            {
                "receptor": rec,
                "seed": seed,
                "arm": arm,
                "stronger": "MZM",
                "weaker": "AZM",
                "MZM_DSF_Kd_uM": 1.2,
                "AZM_DSF_Kd_uM": 7.2,
                "AZM_score_minus_MZM": round(delta, 3),
                "concordance": (1 if delta > 0 else 0.5 if delta == 0 else 0),
            }
        )
    save_rows(work / "holdout_modes.csv", rows)
    save_rows(work / "holdout_cells.csv", cells)
    save_rows(work / "holdout_pairs.csv", pairs)
    # Morkunaite et al., Tables 1-2. This sensitivity is for the 25 C phosphate ITC condition only.
    # Experimental state-population penalties are added once to inter-ligand reported scores.
    offset = anion_pair_offset(first_pka=7.3, second_pka=7.1, ph=7.0, kelvin=298.15)
    linkage = []
    for pair in pairs:
        delta = pair["AZM_score_minus_MZM"]
        corrected = delta + offset
        linkage.append(
            {
                "receptor": pair["receptor"],
                "seed": pair["seed"],
                "arm": pair["arm"],
                "AZM_ITC_Kd_nM": round(1e9 / 1.19e6, 3),
                "MZM_ITC_Kd_nM": round(1e9 / 1.14e7, 3),
                "AZM_score_minus_MZM": delta,
                "differential_penalty_kcal_mol": round(offset, 6),
                "corrected_gap_kcal_mol": round(corrected, 6),
                "concordance_before": pair["concordance"],
                "concordance_after": 1
                if corrected > 0
                else (0.5 if corrected == 0 else 0),
            }
        )
    save_rows(work / "holdout_linkage.csv", linkage)
    save_rows(work / "holdout_rerank.csv", reranked)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "action", choices=("prepare", "native", "anchor", "shell", "evaluate")
    )
    parser.add_argument("work", type=Path)
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--box-size", type=int, choices=(20, 36), default=36)
    args = parser.parse_args()
    if args.action == "prepare":
        prepare(args.work, args.box_size)
    elif args.action == "evaluate":
        evaluate(args.work)
    else:
        if args.binary is None:
            parser.error("docking needs --binary")
        dock(args.work, args.binary, args.action, args.workers)

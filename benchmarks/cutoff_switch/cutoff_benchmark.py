"""Reproduce the T4 lysozyme and held-out BRD4 cutoff experiments in separate directories."""
import argparse
import csv
import hashlib
import itertools
import json
import re
import subprocess
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import gemmi
import numpy as np
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem

# Whole crystallographic series from Morton & Matthews, doi:10.1021/bi00027a007.
# Kd in uM: original ITC series, doi:10.1021/bi00027a006, as tabulated
# by Graves et al. doi:10.1021/jm0491187 (BindingDB's Endolysin [L99A] records).
# All complexes and both shared receptors have the C54T/C97A/L99A sequence.
T4_SERIES = {
    "181L": ("BNZ", "benzene", "c1ccccc1", 175),
    "182L": ("BZF", "benzofuran", "c1ccc2occc2c1", 112),
    "183L": ("DEN", "indene", "C1=Cc2ccccc2C1", 193),
    "184L": ("I4B", "isobutylbenzene", "CC(C)Cc1ccccc1", 19),
    "185L": ("IND", "indole", "c1ccc2[nH]ccc2c1", 290),
    "186L": ("N4B", "n-butylbenzene", "CCCCc1ccccc1", 14),
    "187L": ("PXY", "p-xylene", "Cc1ccc(C)cc1", 422),
    "188L": ("OXE", "o-xylene", "Cc1ccccc1C", 470),
    "1NHB": ("PYJ", "ethylbenzene", "CCc1ccccc1", 68),
}
BRD_SERIES = {
    "4OGI": ("R78", "BI-2536", "CC[C@H]1N(C2CCCC2)c3nc(Nc4ccc(cc4OC)C(=O)NC5CCN(C)CC5)ncc3N(C)C1=O", .037),
    "4OGJ": ("2TA", "TG-101348", "Cc1cnc(Nc2ccc(OCCN3CCCC3)cc2)nc1Nc4cccc(c4)S(=O)(=O)NC(C)(C)C", .164),
    "4O76": ("1M3", "TG-101209", "CN1CCN(CC1)c2ccc(Nc3ncc(C)c(Nc4cccc(c4)S(=O)(=O)NC(C)(C)C)n3)cc2", .123),
    "4O77": ("2RE", "SB-202190", "Oc1ccc(cc1)c2[nH]c(c3ccncc3)c(n2)c4ccc(F)cc4", 3.400),
}
# BRD4 Kd values are from Ciceri et al., doi:10.1038/nchembio.1471.
# Only the 4OGJ chain B was used because chain A lacks several side-chain atoms.
def chain(i):
    return "B" if i == "4OGJ" else "A"


def choose_panel(panel):
    global SERIES, SHARED, ANCHORS, SEEDS, SIDE, ARMS
    if panel == "brd4":
        SERIES, SHARED = BRD_SERIES, ("3MXF", "4OGI", "4OGJ")
        ANCHORS, SIDE = (81, 82, 83, 85, 97, 140, 146), 28
        SEEDS, ARMS = (314159, 271828), ("native", "mixed")
    else:
        SERIES, SHARED = T4_SERIES, ("1L90", "186L")
        ANCHORS, SIDE = (78, 99, 102, 118, 121, 133, 153), 16
        SEEDS, ARMS = (42, 1729), ("native", "smooth")


def coords(line):
    return [float(line[i:i+8]) for i in (30, 38, 46)]


def chain_atoms(lines, i):
    return [x for x in lines if x.startswith("ATOM  ") and x[21] == chain(i)
            and x[76:78].strip() not in ("H", "D") and x[16] in (" ", "A")]


def backbone(lines, i):
    return {int(x[22:26]): coords(x) for x in chain_atoms(lines, i)
            if x[12:16].strip() == "CA"}


def ccd_names(ccd, reference, mol):
    """CCD bond graph and atom names, independent of proximity in either pose."""
    block = gemmi.cif.read_file(str(ccd)).sole_block()
    atoms = block.get_mmcif_category("_chem_comp_atom.")
    bonds = block.get_mmcif_category("_chem_comp_bond.")
    rw = Chem.RWMol()
    names = []
    for name, element, charge in zip(atoms["atom_id"], atoms["type_symbol"], atoms["charge"]):
        atom = Chem.Atom(element.title())
        atom.SetFormalCharge(int(charge))
        atom.SetProp("ccd", name)
        rw.AddAtom(atom)
        names.append(name)
    orders = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE,
              "TRIP": Chem.BondType.TRIPLE}
    for a, b, order in zip(bonds["atom_id_1"], bonds["atom_id_2"], bonds["value_order"]):
        rw.AddBond(names.index(a), names.index(b), orders[order])
    m = rw.GetMol()
    Chem.SanitizeMol(m)
    m = Chem.RemoveHs(m)
    heavy_names = [a.GetProp("ccd") for a in m.GetAtoms()]
    if set(reference) != set(heavy_names):
        raise ValueError(f"CCD atom names disagree with deposited ligand {ccd}")
    if Chem.FindMolChiralCenters(mol):
        conf = Chem.Conformer(m.GetNumAtoms())
        for j, atom in enumerate(heavy_names):
            conf.SetAtomPosition(j, reference[atom])
        m.AddConformer(conf)
        Chem.AssignStereochemistryFrom3D(m)
    if Chem.MolToSmiles(m) != Chem.MolToSmiles(mol):
        raise ValueError(f"Incompatible ligand and crystal: {ccd}")
    return [[heavy_names[i] for i in match] for match in m.GetSubstructMatches(mol, uniquify=False, useChirality=True)]


def download(path, url):
    if path.exists():
        return
    with urllib.request.urlopen(url, timeout=60) as response:
        data = response.read()
    path.write_bytes(data)


def prepare(work):
    raw = work / "raw"
    prep = work / "prepared"
    for d in (raw, prep):
        d.mkdir(parents=True, exist_ok=True)
    ids = sorted(set(SHARED) | SERIES.keys())
    structures = {}
    digests = {}
    for i in ids:
        path = raw / f"{i}.pdb"
        download(path, f"https://files.rcsb.org/download/{i}.pdb")
        digests[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
        lines = path.read_text().splitlines()
        structures[i] = lines
        protein = chain_atoms(lines, i)
        if SERIES is T4_SERIES and {int(x[22:26]): x[17:20] for x in protein if int(x[22:26]) in (54, 97, 99)} != {54: "THR", 97: "ALA", 99: "ALA"}:
            raise ValueError(f"Unexpected construct in {i}")
        protein_file = prep / f"{i}_protein.pdb"
        protein_file.write_text("\n".join(protein) + "\n")
        if SERIES is BRD_SERIES:
            from openmm import Platform
            from openmm.app import PDBFile
            from pdbfixer import PDBFixer
            fixer = PDBFixer(filename=str(protein_file), platform=Platform.getPlatformByName("CPU"))
            fixer.findMissingResidues()
            fixer.missingResidues = {}
            fixer.findMissingAtoms()
            (prep / f"{i}_missing.json").write_text(json.dumps(
                {r.id: [atom.name for atom in atoms] for r, atoms in fixer.missingAtoms.items()}))
            fixer.addMissingAtoms(seed=2026)
            with protein_file.open("w") as fixed:
                PDBFile.writeFile(fixer.topology, fixer.positions, fixed, keepIds=True)
        with (prep / f"{i}_receptor.log").open("w") as log:
            subprocess.run([sys.executable, "-m", "meeko.cli.mk_prepare_receptor",
                            "--read_pdb", str(protein_file), "--default_altloc", "A",
                            "-p", str(prep / f"{i}_receptor.pdbqt")], check=True,
                           stdout=log, stderr=subprocess.STDOUT)
        ca = backbone(lines, i)
        center = np.mean([ca[k] for k in ANCHORS], axis=0)
        (prep / f"{i}_box.txt").write_text("".join(f"center_{k} = {v:.4f}\nsize_{k} = {SIDE}\n" for k, v in zip("xyz", center)))
    metadata = {}
    for i, (component, name, smiles, kd) in SERIES.items():
        mol = Chem.MolFromSmiles(smiles)
        docking_mol = Chem.RWMol(mol)
        if SERIES is BRD_SERIES and component in ("R78", "2TA", "1M3"):
            # Protonate the non-aryl, non-amide tertiary aliphatic cyclic amine.
            pattern = Chem.MolFromSmarts("[N;X3;R;!$(N-C=O);!$(N-c)]")
            sites = docking_mol.GetSubstructMatches(pattern)
            if len(sites) != 1:
                raise ValueError(f"Ambiguous protonation site {component}: {sites}")
            amine = docking_mol.GetAtomWithIdx(sites[0][0])
            amine.SetFormalCharge(1)
            amine.SetNumExplicitHs(1)
        Chem.SanitizeMol(docking_mol)
        supplied = Chem.AddHs(docking_mol)
        if AllChem.EmbedMolecule(supplied, randomSeed=2026) != 0 or AllChem.MMFFOptimizeMolecule(supplied, maxIters=1000) != 0:
            raise RuntimeError(f"Could not generate independent starting conformer for {i}")
        setups = MoleculePreparation().prepare(supplied)
        assert len(setups) == 1
        data, ok, error = PDBQTWriterLegacy.write_string(setups[0], add_index_map=True)
        if not ok:
            raise RuntimeError(error)
        (prep / f"{i}_ligand.pdbqt").write_text(data)
        serial_to_idx = {}
        for line in data.splitlines():
            if line.startswith("REMARK INDEX MAP "):
                vals = list(map(int, line.split()[3:]))
                serial_to_idx.update((serial, source-1) for source, serial in zip(vals[::2], vals[1::2])
                                     if source <= mol.GetNumAtoms())
        if set(serial_to_idx.values()) != set(range(mol.GetNumAtoms())):
            raise ValueError(f"Incomplete ligand atom map for {i}")
        ccd = raw / f"{component}.cif"
        download(ccd, f"https://files.rcsb.org/ligands/download/{component}.cif")
        digests[ccd.name] = hashlib.sha256(ccd.read_bytes()).hexdigest()
        crystal = [x for x in structures[i] if x.startswith("HETATM") and x[21] == chain(i)
                   and x[17:20] == component and x[76:78].strip() != "H"]
        reference = {x[12:16].strip(): coords(x) for x in crystal}
        assert len(crystal) == len(reference) and all(x[16] == " " for x in crystal)
        matches = ccd_names(ccd, reference, mol)
        alignments = {}
        for r in set(SHARED) | {i}:
            moving, fixed = backbone(structures[i], i), backbone(structures[r], r)
            # C-terminal cavity domain. Exclude the known mobile F helix (107-115)
            # so the fit uses receptor framework independently of ligand positions.
            keys = ([k for k in range(80, 161) if not 107 <= k <= 115]
                    if SERIES is T4_SERIES else sorted(k for k in moving.keys() & fixed.keys() if 44 <= k <= 166))
            x, y = (np.array([table[k] for k in keys]) for table in (moving, fixed))
            xc, yc = x.mean(0), y.mean(0)
            u, _, vt = np.linalg.svd((x-xc).T @ (y-yc))
            rotation = u @ np.diag([1, 1, np.linalg.det(u @ vt)]) @ vt
            shift = yc - xc @ rotation
            reference_by_match = [[(np.array(reference[names[serial_to_idx[s]]]) @ rotation + shift).tolist()
                                   for s in sorted(serial_to_idx)] for names in matches]
            center = np.mean([fixed[k] for k in ANCHORS], axis=0)
            margin = SIDE / 2 - np.abs(np.asarray(reference_by_match[0]) - center).max()
            alignments[r] = {"box_margin_A": float(margin), "ca_count": len(keys), "ca_rmsd": float(np.sqrt(np.mean(np.sum((x@rotation+shift-y)**2, axis=1)))),
                             "rotation": rotation.tolist(), "translation": shift.tolist(),
                             "reference_by_match": reference_by_match}
        metadata[i] = {"name": name, "component": component, "kd_uM": kd, "smiles": smiles, "docking_smiles": Chem.MolToSmiles(docking_mol),
                       "heavy_serials": sorted(serial_to_idx), "alignments": alignments}
    (prep / "manifest.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (prep / "sources.sha256").write_text("".join(f"{v}  {k}\n" for k, v in sorted(digests.items())))


SHARED_ONLY = False


def jobs():
    for i in SERIES:
        for r in (SHARED if SHARED_ONLY else dict.fromkeys((*SHARED, i))):
            yield from itertools.product([i], [r], SEEDS, ARMS)


def dock(work, binary, workers):
    runs = work / "runs"
    runs.mkdir(exist_ok=True)
    prep = work / "prepared"
    def run(case):
        i, r, seed, arm = case
        basename = f"{i}_{r}_{seed}_{arm}"
        out = runs / (basename + ".pdbqt")
        cmd = [str(binary.resolve()), "--receptor", str(prep / f"{r}_receptor.pdbqt"),
               "--ligand", str(prep / f"{i}_ligand.pdbqt"), "--config", str(prep / f"{r}_box.txt"),
               "--cpu", "1", "--seed", str(seed), "--exhaustiveness", "4", "--max_evals", "50000",
               "--num_modes", "9", "--energy_range", "100", "--search_counters", "--out", str(out)]
        if ARMS[1] == "mixed":
            cmd += ["--pose_selection", "late"]
        if arm == "mixed":
            cmd.append("--mixed_cutoff_search")
        if arm == "smooth":
            cmd.append("--smooth_cutoff")
        start = time.monotonic()
        with (runs / (basename+".log")).open("w") as log:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
        return basename, time.monotonic()-start
    with ThreadPoolExecutor(workers) as executor:
        for case, wall in executor.map(run, jobs()):
            print(f"{case}: {wall:.2f} s", flush=True)


def evaluate(work):
    meta = json.loads((work / "prepared/manifest.json").read_text())
    rows = []
    for i, r, seed, arm in jobs():
        base = work / "runs" / f"{i}_{r}_{seed}_{arm}"
        models = re.findall(r"^MODEL .*?^ENDMDL$", base.with_suffix(".pdbqt").read_text(), re.M | re.S)
        assert models, base
        entry = meta[i]
        refs = np.array(entry["alignments"][r]["reference_by_match"])
        results = []
        for model in models:
            score = float(re.search(r"REMARK VINA RESULT:\s+(\S+)", model).group(1))
            atoms = {int(x[6:11]): coords(x) for x in model.splitlines() if x.startswith(("ATOM  ", "HETATM"))}
            xyz = np.array([atoms[s] for s in entry["heavy_serials"]])
            rmsd = float(np.sqrt(np.mean(np.sum((refs - xyz)**2, axis=2), axis=1)).min())
            results.append((score, rmsd))
        log = base.with_suffix(".log").read_text()
        counters = next(line for line in log.splitlines() if line.startswith("Search counters:"))
        rows.append(dict(ligand=i, receptor=r, seed=seed, arm=arm, kd_uM=entry["kd_uM"],
                         score=results[0][0], top1_rmsd=results[0][1], best_rmsd=min(x[1] for x in results),
                         n_models=len(models), pool_rmsds=",".join(f"{x[1]:.5f}" for x in results),
                         hunt_evals=sum(map(int, re.findall(r"hunt_evals=(\d+)", log))),
                         second_evals=sum(map(int, re.findall(r"second_evals=(\d+)", log))),
                         counters=counters))
    with (work / "cells.csv").open("w") as out:
        writer = csv.DictWriter(out, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    pairs = []
    for r, seed, arm in itertools.product(SHARED, SEEDS, ARMS):
        block = [row for row in rows if (row["receptor"], row["seed"], row["arm"]) == (r, seed, arm)]
        for a, b in itertools.combinations(sorted(block, key=lambda z: z["kd_uM"]), 2):
            delta = b["score"] - a["score"]
            pairs.append(dict(receptor=r, seed=seed, arm=arm, stronger=a["ligand"], weaker=b["ligand"],
                              kd_ratio=b["kd_uM"]/a["kd_uM"], score_gap=delta,
                              concordance=1 if delta>0 else (0.5 if delta==0 else 0)))
    with (work / "pairs.csv").open("w") as out:
        writer = csv.DictWriter(out, fieldnames=pairs[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(pairs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["prepare", "dock", "evaluate"])
    parser.add_argument("work", type=Path)
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--panel", choices=["t4l","brd4"], default="t4l")
    parser.add_argument("--box-size", type=int, help="Sensitivity box side length, in Angstrom")
    parser.add_argument("--shared-only", action="store_true", help="Only run the fixed set of shared receptors")
    args = parser.parse_args()
    global SIDE, SHARED_ONLY
    choose_panel(args.panel)
    if args.box_size is not None:
        SIDE = args.box_size
    SHARED_ONLY = args.shared_only
    if args.stage == "prepare":
        prepare(args.work)
    elif args.stage == "dock":
        if not args.binary:
            parser.error("dock requires --binary")
        dock(args.work, args.binary, args.workers)
    else:
        evaluate(args.work)


if __name__ == "__main__":
    main()

"""Prepare and search the registered human GBA1 inhibitor series in receptor coordinates."""

import argparse
import hashlib
import itertools
import json
import re
import sys
import time
from importlib.metadata import version
from pathlib import Path

import numpy as np
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "hydration_transfer"))
from panel import cross, timed_process, write_csv

CASES = {"2NSX": ("B", "IFM"), "2V3D": ("A", "NBV"), "2V3E": ("A", "NND")}


def protonate(smiles, dest):
    # Preserve the neutral template's heavy-atom order, including throughout AddHs.
    mol = Chem.MolFromSmiles(smiles)
    (nitrogen,) = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 7]
    nitrogen.SetNumExplicitHs(nitrogen.GetTotalNumHs() + 1)
    nitrogen.SetFormalCharge(1)
    Chem.SanitizeMol(mol)
    full = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 43103
    if AllChem.EmbedMolecule(full, params) or AllChem.MMFFOptimizeMolecule(
        full, maxIters=1000
    ):
        raise ValueError("Protonated graph embedding/relaxation failed")
    (setup,) = MoleculePreparation().prepare(full)
    text, ok, message = PDBQTWriterLegacy.write_string(setup, add_index_map=True)
    if not ok:
        raise ValueError(message)
    dest.write_text(text)
    serial = {}
    for line in text.splitlines():
        if line.startswith("REMARK INDEX MAP"):
            fields = list(map(int, line.split()[3:]))
            serial.update(zip(fields[::2], fields[1::2]))
    return [serial[i + 1] for i in range(mol.GetNumAtoms())], Chem.MolToSmiles(mol)


def select_missing_primary_altloc(lines):
    """Retain a deposited protein conformer when a residue has no A conformer."""
    occupancies = {}
    for line in lines:
        if line.startswith("ATOM  ") and line[16] != " ":
            groups = occupancies.setdefault(line[17:27], {})
            groups.setdefault(line[16], []).append(float(line[54:60]))
    chosen = {
        residue: min(groups, key=lambda label: (-np.mean(groups[label]), label))
        for residue, groups in occupancies.items()
        if "A" not in groups
    }
    adjusted = [
        line[:16] + "A" + line[17:]
        if line.startswith("ATOM  ") and line[16] == chosen.get(line[17:27])
        else line
        for line in lines
    ]
    return adjusted, chosen


def prepare(w, cases=CASES):
    raw = w / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    meta, costs = {}, []
    for code, (chain, comp) in cases.items():
        d = w / code
        d.mkdir(exist_ok=True)
        pdb = raw / (code + ".pdb")
        ccd = raw / (comp + ".cif")
        urls = {
            pdb: f"https://files.rcsb.org/download/{code}.pdb",
            ccd: f"https://files.rcsb.org/ligands/download/{comp}.cif",
        }
        hashes = {u: cross.fetch(p, u) for p, u in urls.items()}
        # The shared prep routine operates on chain A. Relabel only the selected
        # receptor instance; keep source identifiers and all coordinates intact.
        lines = [
            l[:21] + "A" + l[22:]
            for l in pdb.read_text().splitlines()
            if l.startswith(("ATOM  ", "HETATM")) and l[21] == chain
        ]
        lines, altlocs = select_missing_primary_altloc(lines)
        ref = cross.dep_records(lines, comp)
        deposited = np.asarray(list(ref[0].values()))
        # 2V3E's terminal Gln497 ends at CB in the deposition. Building its
        # missing side chain creates a 1.77 Å C--C clash with residue 34.
        # Retain all observed atoms as terminal Ala, far outside the scoring box.
        truncation = None
        if code == "2V3E":
            stub = [
                l for l in lines if l.startswith("ATOM") and l[22:26].strip() == "497"
            ]
            assert {l[12:16].strip() for l in stub} == {"N", "CA", "C", "O", "CB"}
            nearest = min(
                np.linalg.norm(deposited - cross.coords(l), axis=1).min() for l in stub
            )
            assert nearest > 40
            truncation = dict(
                residue="GLN497", representation="ALA", nearest_ligand_A=nearest
            )
            lines = [l[:17] + "ALA" + l[20:] if l in stub else l for l in lines]
        then = time.perf_counter()
        rec = cross.receptor(lines, d, code, deposited)
        costs.append(
            dict(pdb=code, stage="receptor", wall_s=time.perf_counter() - then)
        )
        then = time.perf_counter()
        lig = cross.ligand(ccd, ref, d / "neutral.pdbqt")
        serial, charged = protonate(lig["canonical_smiles"], d / "protonated.pdbqt")
        costs.append(
            dict(pdb=code, stage="both_ligands", wall_s=time.perf_counter() - then)
        )
        ca = {
            l[22:27].strip() + ":" + l[17:20]: cross.coords(l)
            for l in lines
            if l.startswith("ATOM") and l[12:16].strip() == "CA" and l[16] in " A"
        }
        meta[code] = dict(
            chain=chain,
            comp=comp,
            center=deposited.mean(axis=0).tolist(),
            ca=ca,
            sources=hashes,
            remote_terminal_stub=truncation,
            **({"alternate_protein_conformers": altlocs} if altlocs else {}),
            protonated_serial=serial,
            protonated_smiles=charged,
            **rec,
            **lig,
        )
    (w / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    (w / "software.json").write_text(
        json.dumps(
            {
                p: version(p)
                for p in (
                    "pdbfixer",
                    "openmm",
                    "numpy",
                    "scipy",
                    "meeko",
                    "gemmi",
                    "rdkit",
                )
            },
            indent=2,
        )
        + "\n"
    )
    write_csv(w / "preparation.csv", costs)


def dock(w, pro, official):
    directory = w / "runs"
    directory.mkdir(exist_ok=True)
    engine = {"pro": pro.resolve(), "official": official.resolve()}
    (w / "binaries.json").write_text(
        json.dumps(
            {k: hashlib.sha256(p.read_bytes()).hexdigest() for k, p in engine.items()},
            indent=2,
        )
        + "\n"
    )
    meta = json.loads((w / "metadata.json").read_text())
    costs = []
    for src, dst, seed in itertools.product(CASES, CASES, cross.SEEDS):
        for state, arm in (
            ("protonated", "pro"),
            ("protonated", "official"),
            ("neutral", "pro"),
        ):
            tag = f"{src}_{dst}_{state}_{arm}_{seed}"
            output = f"runs/{tag}.pdbqt"
            log = directory / (tag + ".log")
            args = [
                "--ligand",
                f"{src}/{state}.pdbqt",
                "--receptor",
                f"{dst}/receptor.pdbqt",
                "--cpu",
                "1",
                "--exhaustiveness",
                "4",
                "--max_evals",
                "12000",
                "--num_modes",
                "9",
                "--energy_range",
                "100",
                "--seed",
                str(seed),
                "--out",
                output,
            ]
            if arm == "pro":
                args.extend(
                    [
                        "--pose_selection",
                        "late",
                        "--refine_candidates",
                        "40",
                        "--search_counters",
                    ]
                )
            for axis, value in zip("xyz", meta[dst]["center"]):
                args.extend(["--center_" + axis, str(value), "--size_" + axis, "24"])
            try:
                times = timed_process([str(engine[arm]), *args], w, log)
                status = "ok" if (w / output).is_file() else "missing"
            except RuntimeError:
                times = dict(wall_s="", cpu_s="", peak_rss_kib="")
                status = "failed"
            trace = log.read_text()
            search = re.findall(r"hunt_evals=(\d+).*?second_evals=(\d+)", trace)
            explicit = re.findall(r"explicit_evals=(\d+)", trace)
            costs.append(
                dict(
                    src=src,
                    dst=dst,
                    state=state,
                    arm=arm,
                    seed=seed,
                    status=status,
                    **times,
                    search_evals=sum(int(a) + int(b) for a, b in search)
                    if search
                    else "",
                    explicit_evals=sum(map(int, explicit)) if explicit else "",
                    args=" ".join(args),
                )
            )
            write_csv(w / "compute.csv", costs)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    sub = p.add_subparsers(dest="action", required=True)
    sub.add_parser("prepare")
    d = sub.add_parser("dock")
    for k in ("pro", "official"):
        d.add_argument("--" + k, required=True, type=Path)
    a = p.parse_args()
    if a.action == "prepare":
        prepare(a.work.resolve())
    else:
        dock(a.work.resolve(), a.pro, a.official)

"""Prepare Merski's 4W57 F-helix substates and independent alkylbenzene queries."""

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [
    str(ROOT / "benchmarks/anchor_search"),
    str(ROOT / "benchmarks/cutoff_switch"),
    str(ROOT / "tools"),
]
from cutoff_benchmark import backbone, ccd_names, coords, download
from make_scaffold_anchors import all_correspondences, format_assignments
from transfer_mapping import choose_reference

SERIES = {
    "4W52": ("BNZ", 0, -5.19),
    "4W53": ("MBN", 1, -5.52),
    "4W54": ("PYJ", 2, -5.76),
    "4W55": ("3H0", 3, -6.55),
    "4W57": ("N4B", 4, None),
    "4W58": ("3H2", 5, None),
    "4W59": ("3GZ", 6, None),
}
DONOR = "4W57"
QUERIES = tuple(k for k in SERIES if k != DONOR)
STATES = ("A", "B", "C")


def receptor(prep, lines, state):
    kept = {}
    for x in lines:
        if not x.startswith("ATOM  ") or x[21] != "A" or x[76:78].strip() in ("H", "D"):
            continue
        n = int(x[22:26])
        desired = state if 109 <= n <= 113 else "A"
        if x[16] not in (" ", desired):
            continue
        key = n, x[12:16].strip()
        if key in kept:
            raise ValueError(f"Repeated residue/atom {key}")
        kept[key] = x[:16] + " " + x[17:54] + "  1.00" + x[60:]
    path = prep / f"4W57_{state}.pdb"
    path.write_text("\n".join(kept.values()) + "\nEND\n")
    fixed = prep / "4W57_common.pdb"
    if state == "A":
        from openmm import Platform
        from openmm.app import PDBFile
        from pdbfixer import PDBFixer

        fixer = PDBFixer(filename=str(path), platform=Platform.getPlatformByName("CPU"))
        fixer.findMissingResidues()
        fixer.missingResidues = {}
        fixer.findMissingAtoms()
        (prep / "missing_atoms.json").write_text(
            json.dumps(
                {
                    r.id: [a.name for a in atoms]
                    for r, atoms in fixer.missingAtoms.items()
                }
            )
        )
        if any(109 <= int(r.id) <= 113 for r in fixer.missingAtoms):
            raise ValueError("Incomplete helix")
        fixer.addMissingAtoms(seed=2026)
        with fixed.open("w") as output:
            PDBFile.writeFile(fixer.topology, fixer.positions, output, keepIds=True)
    output = []
    for line in fixed.read_text().splitlines():
        if line.startswith("ATOM  ") and 109 <= int(line[22:26]) <= 113:
            line = (
                line[:30]
                + kept[int(line[22:26]), line[12:16].strip()][30:54]
                + line[54:]
            )
        if (
            line.startswith("ATOM  ")
            and line[12:16].strip() == "OXT"
            and (164, "OXT") not in kept
        ):
            # Place the unobserved terminal oxygen in the carboxylate plane.
            c, ca, o = (np.array(coords(kept[164, name])) for name in ("C", "CA", "O"))
            u = (c - ca) / np.linalg.norm(c - ca) + (c - o) / np.linalg.norm(c - o)
            xyz = c + 1.25 * u / np.linalg.norm(u)
            line = line[:30] + "".join(f"{v:8.3f}" for v in xyz) + line[54:]
        output.append(line)
    path.write_text("\n".join(output) + "\n")
    with (prep / f"{state}_prepare.log").open("w") as log:
        subprocess.run(
            [
                sys.executable,
                "-m",
                "meeko.cli.mk_prepare_receptor",
                "--read_pdb",
                str(path),
                "-p",
                str(prep / f"{state}_receptor.pdbqt"),
            ],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    return kept


def ligand(prep, raw, code):
    component, length, _ = SERIES[code]
    mol = Chem.MolFromSmiles("C" * length + "c1ccccc1")
    supplied = Chem.AddHs(mol)
    if (
        AllChem.EmbedMolecule(supplied, randomSeed=2026) != 0
        or AllChem.MMFFOptimizeMolecule(supplied, maxIters=1000) != 0
    ):
        raise ValueError(f"Independent initial structure failed for {code}")
    (setup,) = MoleculePreparation().prepare(supplied)
    text, ok, error = PDBQTWriterLegacy.write_string(setup, add_index_map=True)
    if not ok:
        raise ValueError(error)
    (prep / f"{code}_ligand.pdbqt").write_text(text)
    index_to_serial = {}
    for line in text.splitlines():
        if line.startswith("REMARK INDEX MAP "):
            vals = list(map(int, line.split()[3:]))
            index_to_serial.update(
                (idx - 1, s)
                for idx, s in zip(vals[::2], vals[1::2], strict=True)
                if idx <= mol.GetNumAtoms()
            )
    if set(index_to_serial) != set(range(mol.GetNumAtoms())):
        raise ValueError("Missing ligand atom")
    lines = (raw / f"{code}.pdb").read_text().splitlines()
    deposits = [
        x
        for x in lines
        if x.startswith("HETATM")
        and x[21:26] == "A 200"
        and x[17:20].strip() == component
        and x[76:78].strip() != "H"
    ]
    alts = sorted({x[16] for x in deposits} - {" "}) or [" "]
    refs = {}
    for alt in alts:
        selected = [x for x in deposits if x[16] in (" ", alt)]
        xyz = {x[12:16].strip(): coords(x) for x in selected}
        if len(xyz) != len(selected) or len(xyz) != mol.GetNumAtoms():
            raise ValueError(f"Incomplete {code}/{alt}")
        names = ccd_names(raw / f"{component}.cif", xyz, mol)
        refs[alt.strip() or "."] = {
            "occupancy": float(np.mean([float(x[54:60]) for x in selected])),
            "mapped": [np.array([xyz[n] for n in match]) for match in names],
        }
    return mol, index_to_serial, refs


def prepare(work):
    raw = work / "raw"
    prep = work / "prepared"
    raw.mkdir(exist_ok=True, parents=True)
    prep.mkdir(exist_ok=True)
    for code, (component, _, _) in SERIES.items():
        download(raw / f"{code}.pdb", f"https://files.rcsb.org/download/{code}.pdb")
        download(
            raw / f"{component}.cif",
            f"https://files.rcsb.org/ligands/download/{component}.cif",
        )
    proteins = {code: (raw / f"{code}.pdb").read_text().splitlines() for code in SERIES}
    recs = [receptor(prep, proteins[DONOR], state) for state in STATES]
    for rec in recs[1:]:
        if rec.keys() != recs[0].keys() or any(
            rec[k] != recs[0][k] for k in rec if k[0] not in range(109, 114)
        ):
            raise ValueError("Receptor states differ outside declared F-helix")
    ca = backbone(proteins[DONOR], DONOR)
    center = np.mean([ca[k] for k in (78, 99, 102, 118, 121, 133, 153)], axis=0)
    (prep / "box.txt").write_text(
        "".join(f"center_{k} = {v:.4f}\nsize_{k} = 20\n" for k, v in zip("xyz", center))
    )
    ligands = {code: ligand(prep, raw, code) for code in SERIES}
    donor, _, donor_refs = ligands[DONOR]
    donor_xyz = donor_refs["A"]["mapped"][0]
    manifest = {}
    for code in QUERIES:
        mol, idx, refs = ligands[code]
        _, qa, da = choose_reference(code, {code: mol, DONOR: donor})
        single = tuple(sorted(zip(qa, da, strict=True)))
        mappings = all_correspondences(mol, donor)
        if single not in mappings:
            raise ValueError("Unrecoverable baseline assignment")
        (prep / f"{code}.single").write_text(
            format_assignments([zip(qa, da, strict=True)], idx, donor_xyz)
        )
        (prep / f"{code}.multi").write_text(
            format_assignments(mappings, idx, donor_xyz)
        )
        frame = [k for k in range(80, 161) if not 107 <= k <= 115]
        moving = backbone(proteins[code], code)
        a, b = (np.array([v[k] for k in frame]) for v in (moving, ca))
        u, _, vt = np.linalg.svd((a - a.mean(0)).T @ (b - b.mean(0)))
        rot = u @ np.diag([1, 1, np.linalg.det(u @ vt)]) @ vt
        trans = b.mean(0) - a.mean(0) @ rot
        for data in refs.values():
            transformed = np.array(data["mapped"]) @ rot + trans
            data["mapped"] = transformed.tolist()
            data["box_margin"] = float(np.min(10 - np.abs(transformed - center)))
        manifest[code] = {
            "refs": refs,
            "serials": [idx[j] for j in range(len(idx))],
            "dg": SERIES[code][2],
            "single": single,
            "mapping_count": len(mappings),
            "rotation": rot.tolist(),
            "translation": trans.tolist(),
        }
    (prep / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
    (prep / "sources.sha256").write_text(
        "".join(
            hashlib.sha256(f.read_bytes()).hexdigest() + "  " + f.name + "\n"
            for f in sorted(raw.iterdir())
            if f.is_file()
        )
    )
    print([(k, v["mapping_count"], list(v["refs"])) for k, v in manifest.items()])


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    prepare(p.parse_args().work)

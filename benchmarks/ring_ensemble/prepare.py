"""Build a locked new-human-target small-ring pose panel from wwPDB deposits."""

import argparse
import json
import os
import subprocess
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import gemmi
import numpy as np
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
from ring_conformers import fit, prepare_pair

# target, query component, crystal/chain, other receptor/chain/local marker (if available).
PANEL = (
    ("HSP90a", "2GJ", "2VCI", "A", "2XAB", "A", "VHD"),
    ("DPP4", "715", "1X70", "A", "2RGU", "A", "356"),
    ("p38a", "B96", "1KV2", "A", "3GCP", "A", "SB2"),
    ("CDK6", "LQQ", "5L2I", "A", "2EUF", "B", "LQQ"),
    ("GR", "DEX", "1M2Z", "A", None, None, None),
    ("PR", "STR", "1A28", "A", None, None, None),
    ("PARP1", "09L", "5DS3", "A", "4ZZZ", "A", "FSU"),
    ("PDE5", "VIA", "1UDT", "A", "2H42", "A", "VIA"),
)


def get(raw, name, url):
    path = raw / name
    if not path.exists():
        with urllib.request.urlopen(url, timeout=60) as response:
            data = response.read()
        path.write_bytes(data)
    return path


def ccd_graph(path):
    atoms, bonds = (
        gemmi.cif.read_file(str(path)).sole_block().get_mmcif_category(cat)
        for cat in ("_chem_comp_atom.", "_chem_comp_bond.")
    )
    rw = Chem.RWMol()
    names = atoms["atom_id"]
    ideal = all(
        v is not None
        for axis in "xyz"
        for v in atoms["pdbx_model_Cartn_" + axis + "_ideal"]
    )
    conf = Chem.Conformer(len(names)) if ideal else None
    for i, (symbol, charge, name) in enumerate(
        zip(atoms["type_symbol"], atoms["charge"], names, strict=True)
    ):
        a = Chem.Atom(symbol.title())
        a.SetFormalCharge(int(charge))
        a.SetProp("ccd_name", name)
        rw.AddAtom(a)
        if ideal:
            conf.SetAtomPosition(
                i,
                [
                    float(atoms["pdbx_model_Cartn_" + axis + "_ideal"][i])
                    for axis in "xyz"
                ],
            )
    orders = {
        "SING": Chem.BondType.SINGLE,
        "DOUB": Chem.BondType.DOUBLE,
        "TRIP": Chem.BondType.TRIPLE,
    }
    for a, b, order in zip(
        bonds["atom_id_1"], bonds["atom_id_2"], bonds["value_order"], strict=True
    ):
        rw.AddBond(names.index(a), names.index(b), orders[order])
    m = rw.GetMol()
    Chem.SanitizeMol(m)
    if ideal:
        m.AddConformer(conf)
        Chem.AssignStereochemistryFrom3D(m)
    m = Chem.RemoveHs(m)
    m.RemoveAllConformers()
    descriptors = (
        gemmi.cif.read_file(str(path))
        .sole_block()
        .get_mmcif_category("_pdbx_chem_comp_descriptor.")
    )
    stereo = {
        Chem.MolToSmiles(Chem.MolFromSmiles(v))
        for k, v in zip(descriptors["type"], descriptors["descriptor"])
        if k == "SMILES_CANONICAL"
    }
    if stereo != {Chem.MolToSmiles(m)}:
        raise ValueError(
            f"CCD graph/stereochemistry differs from the deposited isomeric descriptors in {path}"
        )
    return m


def xyz(line):
    return [float(line[i : i + 8]) for i in (30, 38, 46)]


def records(path, chain, marker):
    lines = path.read_text().splitlines()
    heavy = [
        l
        for l in lines
        if l[21:22] == chain
        and l[16:17] in (" ", "A")
        and l[76:78].strip() not in ("H", "D")
    ]
    protein = [l[:16] + " " + l[17:] for l in heavy if l.startswith("ATOM  ")]
    lig = {
        l[12:16].strip(): xyz(l)
        for l in heavy
        if l.startswith("HETATM") and l[17:20].strip() == marker
    }
    if not lig:
        raise ValueError((path, chain, marker))
    return lines, protein, lig


class IncompletePocketError(ValueError):
    pass


def receptor(path, folder, pdb, chain, marker):
    lines, protein, marker_atoms = records(path, chain, marker)
    folder.mkdir(exist_ok=True, parents=True)
    site = np.array(list(marker_atoms.values()))
    ca = {(l[22:27], l[17:20]): xyz(l) for l in protein if l[12:16].strip() == "CA"}
    site_ca = {
        k: v for k, v in ca.items() if np.linalg.norm(site - v, axis=1).min() <= 12
    }
    protein_path = folder / "protein.pdb"
    protein_path.write_text("\n".join(protein) + "\nEND\n")
    qt_path = folder / "receptor.pdbqt"
    with (folder / "prepare.log").open("w") as log:
        subprocess.run(
            [
                sys.executable,
                "-m",
                "meeko.cli.mk_prepare_receptor",
                "--read_pdb",
                str(protein_path),
                "--allow_bad_res",
                "-p",
                str(qt_path),
            ],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    qt = qt_path.read_text()
    kept = {l[22:27] for l in qt.splitlines() if l.startswith(("ATOM", "HETATM"))}
    omitted = [l for l in protein if l[22:27] not in kept]
    near = min(
        (np.linalg.norm(site - xyz(l), axis=1).min() for l in omitted), default=None
    )
    if near is not None and near < 8:
        # Restore missing atoms in incomplete observed sidechains; leave unobserved loops unresolved.
        from openmm import Platform
        from openmm.app import PDBFile
        from pdbfixer import PDBFixer

        fixer = PDBFixer(
            filename=str(protein_path), platform=Platform.getPlatformByName("CPU")
        )
        fixer.findMissingResidues()
        fixer.missingResidues = {}
        fixer.findMissingAtoms()
        fixer.addMissingAtoms(seed=20261126)
        with (folder / "repaired.pdb").open("w") as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
        # Add modelled atoms only to observed residues around the pocket. A distant incomplete
        # terminus can acquire implausible rebuilt bonds and is irrelevant to these grids.
        repaired = (folder / "repaired.pdb").read_text().splitlines()
        nearby = {
            l[22:27] for l in omitted if np.linalg.norm(site - xyz(l), axis=1).min() < 8
        }
        missing_keys = {(l[22:27], l[12:16]) for l in protein}
        added = [
            l
            for l in repaired
            if l.startswith("ATOM")
            and l[22:27] in nearby
            and (l[22:27], l[12:16]) not in missing_keys
        ]
        combined = []
        for i, l in enumerate(protein):
            combined.append(l)
            if i + 1 == len(protein) or protein[i + 1][22:27] != l[22:27]:
                combined.extend(a for a in added if a[22:27] == l[22:27])
        (folder / "repaired.pdb").write_text("\n".join(combined) + "\nEND\n")
        with (folder / "repair.log").open("w") as log:
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "meeko.cli.mk_prepare_receptor",
                    "--read_pdb",
                    str(folder / "repaired.pdb"),
                    "--allow_bad_res",
                    "-p",
                    str(qt_path),
                ],
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
            )
        qt = qt_path.read_text()
        kept = {l[22:27] for l in qt.splitlines() if l.startswith(("ATOM", "HETATM"))}
        omitted = [l for l in protein if l[22:27] not in kept]
        near = min(
            (np.linalg.norm(site - xyz(l), axis=1).min() for l in omitted), default=None
        )
        if near is not None and near < 8:
            raise IncompletePocketError(
                f"{pdb}: omitted pocket residue at {near:.2f} A"
            )
    metals = [
        l
        for l in lines
        if l.startswith("HETATM")
        and l[21] == chain
        and l[16] in (" ", "A")
        and l[76:78].strip() in ("MG", "ZN")
    ]
    # Keep native catalytic metal centers. Vina uses atom types for its standard potential.
    for l in metals:
        qt += l[:54] + "  1.00  0.00     2.000 " + l[76:78].strip().title() + "\n"
    qt_path.write_text(qt)
    (folder / "box.txt").write_text(
        "".join(
            f"center_{k} = {v:.5f}\nsize_{k} = 26\n"
            for k, v in zip("xyz", site.mean(0))
        )
    )
    return (
        ca,
        site_ca,
        {
            "pdb": pdb,
            "chain": chain,
            "box_center": site.mean(0).tolist(),
            "ignored_residue_min_distance_A": near,
            "metals": [l[17:20].strip() for l in metals],
        },
    )


def prepare(work, item):
    target, cc, q, ch, r, rc, marker = item
    raw = work / "raw"
    raw.mkdir(parents=True, exist_ok=True)

    def pp(pdb):
        return get(raw, f"{pdb}.pdb", f"https://files.rcsb.org/download/{pdb}.pdb")

    cpath = get(raw, f"{cc}.cif", f"https://files.rcsb.org/ligands/download/{cc}.cif")
    mol = ccd_graph(cpath)
    folder = work / "prepared" / target
    pair = prepare_pair(mol, folder)
    names = [a.GetProp("ccd_name") for a in mol.GetAtoms()]
    _, _, lig = records(pp(q), ch, cc)
    if set(names) != set(lig):
        raise ValueError(
            f"{target}: incomplete heavy-atom correspondence: {set(names) ^ set(lig)}"
        )
    x = np.array([lig[n] for n in names])
    perms = mol.GetSubstructMatches(mol, uniquify=False, useChirality=True)
    qca, _, record = receptor(pp(q), folder / q, q, ch, cc)
    record.update(reference=x.tolist(), alignment_n=0, alignment_rmsd=0, kind="cognate")
    receptors = {q: record}
    exclusions = {}
    if r:
        try:
            _, site_ca, record = receptor(pp(r), folder / r, r, rc, marker)
        except IncompletePocketError as exc:
            exclusions[r] = str(exc)
    if r and r not in exclusions:
        keys = sorted(set(site_ca) & set(qca))
        a = np.array([qca[k] for k in keys])
        b = np.array([site_ca[k] for k in keys])
        if len(keys) < 12:
            raise ValueError(f"Too few shared pocket CA {target}: {len(keys)}")
        rot, shift = fit(a, b)
        error = np.mean(np.sum((a @ rot + shift - b) ** 2, axis=1)) ** 0.5
        record.update(
            reference=(x @ rot + shift).tolist(),
            alignment_n=len(keys),
            alignment_rmsd=error,
            kind="other_deposit_same_ligand" if cc == marker else "cross_ligand",
        )
        receptors[r] = record
    data = {
        "target": target,
        "ccd": cc,
        "crystal": q,
        "names": names,
        "permutations": perms,
        "states": pair,
        "receptors": receptors,
        "exclusions": exclusions,
    }
    (folder / "panel.json").write_text(json.dumps(data, indent=2) + "\n")
    print(target, pair["alternate_distinct"], pair["ring_distance_A"], flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("work", type=Path)
    p.add_argument("--workers", type=int, default=2)
    p.add_argument("--target", nargs="*")
    args = p.parse_args()
    os.environ.setdefault("OPENMM_CPU_THREADS", "1")
    selected = [row for row in PANEL if not args.target or row[0] in args.target]
    with ThreadPoolExecutor(args.workers) as pool:
        for future in [pool.submit(prepare, args.work, row) for row in selected]:
            future.result()

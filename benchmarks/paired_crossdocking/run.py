"""Compare two Late search engines in reciprocal cross-docking on human drug targets."""

import argparse
import csv
import hashlib
import itertools
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen

import gemmi
import numpy as np
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem

PANELS = {
    "ESR1": {"1A52": "EST", "3ERD": "DES"},
    "PPARG": {"2PRG": "BRL", "5Y2O": "8N6"},
    "PDE4D": {"1XOQ": "ROF", "1XOM": "CIO"},
}
SEEDS = (43103, 43117, 43133)
SIDE = 24
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "benchmarks/cutoff_switch"))
from cutoff_benchmark import coords


def fetch(path, url):
    if not path.is_file():
        data = urlopen(url, timeout=60).read()
        path.write_bytes(data)
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dep_records(lines, code):
    atoms = [
        l
        for l in lines
        if l.startswith("HETATM")
        and l[21] == "A"
        and l[17:20].strip() == code
        and l[76:78].strip() not in {"H", "D"}
    ]
    alts = sorted({l[16] for l in atoms} - {" "}) or [" "]
    return [
        {l[12:16].strip(): coords(l) for l in atoms if l[16] in (" ", alt)}
        for alt in alts
    ]


def ligand(ccd, refs, out):
    block = gemmi.cif.read_file(str(ccd)).sole_block()
    desc = block.get_mmcif_category("_pdbx_chem_comp_descriptor.")
    smi = next(
        s for t, s in zip(desc["type"], desc["descriptor"]) if t == "SMILES_CANONICAL"
    )
    identity = Chem.MolFromSmiles(smi)
    mol = Chem.MolFromSmiles(Chem.MolToSmiles(identity))
    # Generate the starting coordinates from the chemical graph alone.
    start = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 43103
    if AllChem.EmbedMolecule(start, params) or AllChem.MMFFOptimizeMolecule(
        start, maxIters=1000
    ):
        raise ValueError(f"Embedding/relaxation failed: {ccd}")
    (setup,) = MoleculePreparation().prepare(start)
    pdbqt, ok, msg = PDBQTWriterLegacy.write_string(setup, add_index_map=True)
    if not ok:
        raise ValueError(msg)
    out.write_text(pdbqt)
    serial = {}
    for l in pdbqt.splitlines():
        if l.startswith("REMARK INDEX MAP"):
            v = list(map(int, l.split()[3:]))
            serial.update(
                (i - 1, j) for i, j in zip(v[::2], v[1::2]) if i <= mol.GetNumAtoms()
            )
    if set(serial) != set(range(mol.GetNumAtoms())):
        raise ValueError("Incomplete PDBQT atom map")
    atoms = block.get_mmcif_category("_chem_comp_atom.")
    bonds = block.get_mmcif_category("_chem_comp_bond.")
    rw = Chem.RWMol()
    ids = {}
    for name, element, charge in zip(
        atoms["atom_id"], atoms["type_symbol"], atoms["charge"]
    ):
        atom = Chem.Atom(element.title())
        atom.SetFormalCharge(int(charge))
        atom.SetProp("ccd", name)
        ids[name] = rw.AddAtom(atom)
    for a, b, o in zip(bonds["atom_id_1"], bonds["atom_id_2"], bonds["value_order"]):
        rw.AddBond(
            ids[a],
            ids[b],
            {
                "SING": Chem.BondType.SINGLE,
                "DOUB": Chem.BondType.DOUBLE,
                "TRIP": Chem.BondType.TRIPLE,
            }[o],
        )
    graph = rw.GetMol()
    Chem.SanitizeMol(graph)
    graph = Chem.RemoveHs(graph)
    names = [a.GetProp("ccd") for a in graph.GetAtoms()]
    maps = []
    for ref in refs:
        if set(names) != set(ref):
            raise ValueError(f"Incomplete crystal ligand {ccd}")
        geom = Chem.Mol(graph)
        conf = Chem.Conformer(len(names))
        for i, n in enumerate(names):
            conf.SetAtomPosition(i, ref[n])
        geom.AddConformer(conf)
        Chem.AssignStereochemistryFrom3D(geom)
        matches = geom.GetSubstructMatches(
            mol, uniquify=False, useChirality=True, maxMatches=10001
        )
        if not matches or len(matches) > 10000:
            raise ValueError(f"Incompatible reference chemistry {ccd}")
        maps.extend([[ref[names[j]] for j in match] for match in matches])
    return {
        "canonical_smiles": Chem.MolToSmiles(mol),
        "serial": [serial[i] for i in range(len(serial))],
        "mapped": maps,
        "n_altloc": len(refs),
    }


def receptor(lines, out, code, center):
    from openmm import Platform
    from openmm.app import PDBFile
    from pdbfixer import PDBFixer

    protein = [
        l[:16] + " " + l[17:]
        for l in lines
        if l.startswith("ATOM  ")
        and l[21] == "A"
        and l[16] in " A"
        and l[76:78].strip() not in ("H", "D")
    ]
    # A lone terminal atom cannot define a peptide backbone. Remove only such
    # terminal stubs well outside the pocket, recording every exclusion.
    by_res = {}
    for l in protein:
        by_res.setdefault(l[22:27], []).append(l)
    trimmed = []
    for key in dict.fromkeys((next(iter(by_res)), next(reversed(by_res)))):
        ls = by_res[key]
        if len({"N", "CA", "C"} & {l[12:16].strip() for l in ls}) < 2:
            nearest = min(
                np.linalg.norm(np.array(center) - coords(l), axis=1).min() for l in ls
            )
            if nearest < 10:
                raise ValueError("Unresolved terminal residue near pocket")
            trimmed.append(key.strip())
            protein = [l for l in protein if l[22:27] != key]
    pdb = out / "chain-A.pdb"
    pdb.write_text("\n".join(protein) + "\nEND\n")
    f = PDBFixer(filename=str(pdb), platform=Platform.getPlatformByName("CPU"))
    f.findMissingResidues()
    f.missingResidues = {}
    f.findMissingAtoms()
    missing = {r.id: [a.name for a in atoms] for r, atoms in f.missingAtoms.items()}
    f.addMissingAtoms(seed=43103)
    with pdb.open("w") as h:
        PDBFile.writeFile(f.topology, f.positions, h, keepIds=True)
    lines_after = pdb.read_text().splitlines()
    for i, l in enumerate(lines_after):
        if (
            l.startswith("ATOM")
            and l[12:16].strip() == "OXT"
            and not any(
                orig[12:16].strip() == "OXT" and orig[22:27] == l[22:27]
                for orig in protein
            )
        ):
            geom = {
                orig[12:16].strip(): np.array(coords(orig))
                for orig in lines_after
                if orig.startswith("ATOM") and orig[22:27] == l[22:27]
            }
            c = geom["C"]
            u, v = c - geom["CA"], c - geom["O"]
            vec = u / np.linalg.norm(u) + v / np.linalg.norm(v)
            pos = c + 1.25 * vec / np.linalg.norm(vec)
            lines_after[i] = l[:30] + "".join(f"{v:8.3f}" for v in pos) + l[54:]
    pdb.write_text("\n".join(lines_after) + "\n")
    repaired = [
        l
        for l in pdb.read_text().splitlines()
        if l.startswith("ATOM")
        and l[12:16].strip() in missing.get(l[22:26].strip(), [])
    ]
    mindist = min(
        (
            float(np.linalg.norm(np.array(coords(l)) - v, axis=1).min())
            for l in repaired
            for v in [np.array(center)]
        ),
        default=None,
    )
    with (out / "prep.log").open("w") as log:
        subprocess.run(
            [
                sys.executable,
                "-m",
                "meeko.cli.mk_prepare_receptor",
                "--read_pdb",
                str(pdb),
                "-p",
                str(out / "receptor.pdbqt"),
            ],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    # Retain both measured catalytic metal coordinates for PDE4D.
    metals = [
        l
        for l in lines
        if l.startswith("HETATM") and l[21] == "A" and l[17:20].strip() in {"ZN", "MG"}
    ]
    if code in PANELS["PDE4D"] and {l[17:20].strip() for l in metals} != {"MG", "ZN"}:
        raise ValueError("Metal missing")
    with (out / "receptor.pdbqt").open("a") as f:
        for i, l in enumerate(metals, 1):
            f.write(
                l[:6]
                + f"{9000 + i:5d}"
                + l[11:66]
                + f"    {2.0:6.3f} "
                + l[17:20].strip().title()
                + "\n"
            )
    return {
        "trimmed_remote_terminal_stubs": trimmed,
        "missing_atoms": missing,
        "closest_repaired_ligand_distance": mindist,
        "metals": [l[17:20].strip() for l in metals],
    }


def prepare(work):
    raw = work / "raw"
    raw.mkdir(exist_ok=True, parents=True)
    meta = {}
    for target, entries in PANELS.items():
        for code, comp in entries.items():
            out = work / code
            out.mkdir(exist_ok=True)
            cif, pdb = raw / (comp + ".cif"), raw / (code + ".pdb")
            digest = [
                fetch(pdb, f"https://files.rcsb.org/download/{code}.pdb"),
                fetch(cif, f"https://files.rcsb.org/ligands/download/{comp}.cif"),
            ]
            lines = pdb.read_text().splitlines()
            refs = dep_records(lines, comp)
            lig = ligand(cif, refs, out / "ligand.pdbqt")
            center = np.mean(list(refs[0].values()), axis=0)
            rec = receptor(lines, out, code, list(refs[0].values()))
            ca = {
                l[22:27].strip() + ":" + l[17:20]: coords(l)
                for l in lines
                if l.startswith("ATOM  ")
                and l[21] == "A"
                and l[12:16].strip() == "CA"
                and l[16] in " A"
            }
            meta[code] = dict(
                target=target,
                ligand=comp,
                sha256=digest,
                **lig,
                **rec,
                ca=ca,
                center=center.tolist(),
            )
            (work / "metadata.json").write_text(json.dumps(meta))


def alignment(source, receiver):
    # Only protein atoms determine the transform. The receiver's known ligand
    # defines its pocket, including for reciprocal cross-docking.
    keys = sorted(
        k
        for k in source["ca"].keys() & receiver["ca"].keys()
        if np.linalg.norm(np.array(receiver["ca"][k]) - receiver["center"]) < 14
        and np.linalg.norm(np.array(source["ca"][k]) - source["center"]) < 14
    )
    if len(keys) < 12:
        raise ValueError("Insufficient pocket alignment")
    x, y = (np.array([m["ca"][k] for k in keys]) for m in (source, receiver))
    u, _s, vt = np.linalg.svd((x - x.mean(0)).T @ (y - y.mean(0)))
    u[:, -1] *= np.sign(np.linalg.det(u @ vt))
    rot = u @ vt
    trans = y.mean(0) - x.mean(0) @ rot
    rms = float(np.sqrt(np.mean(np.sum((x @ rot + trans - y) ** 2, axis=1))))
    if rms > 2.0:
        raise ValueError("Pocket backbone alignment RMSD > 2 Å")
    return rot, trans, rms, len(keys)


def run(work, binary, mixed_binary, jobs, cap, seeds):
    meta = json.loads((work / "metadata.json").read_text())
    todo = [
        (s, r, seed, arm)
        for entries in PANELS.values()
        for s, r in itertools.product(entries, repeat=2)
        for seed, arm in itertools.product(
            seeds, ("native", "soft", "mixed") if mixed_binary else ("native", "soft")
        )
    ]

    def one(item):
        s, r, seed, arm = item
        out = work / f"runs_{cap}" / f"{s}_{r}_{seed}_{arm}"
        out.parent.mkdir(exist_ok=True)
        rot, trans, backbone, n = alignment(meta[s], meta[r])
        ref = np.array(meta[s]["mapped"]) @ rot + trans
        cmd = [
            str(mixed_binary if arm == "mixed" else binary),
            "--receptor",
            str(work / r / "receptor.pdbqt"),
            "--ligand",
            str(work / s / "ligand.pdbqt"),
            "--pose_selection",
            "late",
            "--exhaustiveness",
            "4",
            "--max_evals",
            str(cap),
            "--num_modes",
            "10",
            "--refine_candidates",
            "40",
            "--energy_range",
            "100",
            "--cpu",
            "1",
            "--search_counters",
            "--seed",
            str(seed),
            "--out",
            str(out.with_suffix(".pdbqt")),
        ]
        for axis, pos in zip("xyz", meta[r]["center"]):
            cmd.extend(["--center_" + axis, str(pos), "--size_" + axis, str(SIDE)])
        if arm == "mixed":
            cmd.append("--mixed_vinardo_search")
        if arm == "soft":
            cmd.append("--soft_receptor_search")
        res = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
            text=True,
            timeout=180,
        )
        out.with_suffix(".log").write_text(res.stdout)
        res.check_returncode()
        rows = []
        for block in out.with_suffix(".pdbqt").read_text().split("ENDMDL"):
            atom = {
                int(l[6:11]): coords(l)
                for l in block.splitlines()
                if l.startswith(("ATOM  ", "HETATM"))
            }
            if not atom:
                continue
            xyz = np.array([atom[k] for k in meta[s]["serial"]])
            rms = float(np.sqrt(((ref - xyz) ** 2).sum(axis=2).mean(axis=1)).min())
            score = float(
                next(
                    l.split()[3]
                    for l in block.splitlines()
                    if l.startswith("REMARK VINA RESULT:")
                )
            )
            rows.append(
                {
                    "target": meta[s]["target"],
                    "ligand": s,
                    "receiver": r,
                    "kind": "redock" if r == s else "cross",
                    "seed": seed,
                    "arm": arm,
                    "rank": len(rows) + 1,
                    "rmsd": rms,
                    "score": score,
                    "ca_rmsd": backbone,
                    "n_ca": n,
                }
            )
        if not rows:
            raise ValueError("No modes")
        return rows

    with ThreadPoolExecutor(jobs) as pool:
        rows = [r for group in pool.map(one, todo) for r in group]
    with (work / f"poses_{cap}.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("--prepare", action="store_true")
    p.add_argument("--vina", type=Path)
    p.add_argument("--mixed-vina", type=Path)
    p.add_argument("--cap", type=int, default=3000)
    p.add_argument("--seeds", type=int, nargs="+", default=SEEDS)
    p.add_argument("--jobs", type=int, default=2)
    a = p.parse_args()
    if a.prepare:
        prepare(a.work)
    if a.vina:
        run(
            a.work,
            a.vina.resolve(),
            a.mixed_vina.resolve() if a.mixed_vina else None,
            a.jobs,
            a.cap,
            a.seeds,
        )

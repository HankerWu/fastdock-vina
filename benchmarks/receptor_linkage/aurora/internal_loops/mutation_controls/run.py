"""Run same-receptor reference engines and independently measure all published poses."""

import argparse
import json
import re
import sys
from pathlib import Path

import gemmi
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import check_predictions as audit
import replay
from frame_check import positions, reference_graph, untyped_graph

SEEDS = (71261, 71263)
OFFICIAL = "f31f774f723bba7bbe6e9d1c47577020eea9a8da16424284c043d22593570644"


def brd_input(brd, out):
    path = brd / "6YQW/ligand.pdbqt"
    lines = path.read_text().splitlines()
    smiles = next(
        line.split()[2] for line in lines if line.startswith("REMARK SMILES ")
    )
    molecule = Chem.MolFromSmiles(smiles)
    by_atom = {}
    edges = []
    for row in lines:
        if row.startswith("REMARK SMILES IDX "):
            v = list(map(int, row.split()[3:]))
            by_atom.update(zip(v[::2], v[1::2], strict=True))
        elif row.startswith("REMARK H PARENT "):
            v = list(map(int, row.split()[3:]))
            edges += [(by_atom[i], j) for i, j in zip(v[::2], v[1::2], strict=True)]
    for b in molecule.GetBonds():
        edges.append((by_atom[b.GetBeginAtomIdx() + 1], by_atom[b.GetEndAtomIdx() + 1]))
    sidecar = out / "82I.bonds"
    sidecar.write_text(
        "".join(f"{i} {j}\n" for i, j in sorted(tuple(sorted(e)) for e in edges))
    )
    mol = untyped_graph(molecule)
    serials = [by_atom[i + 1] for i in range(mol.GetNumAtoms())]
    mol.AddConformer(Chem.Conformer(len(serials)))
    # Gemmi fit uses only observed pocket backbone atoms. Chain-A numbering of
    # the common BRD9 domain differs by 116 (verified by residue name).
    protein = {}
    for name, ligand in (("6YQW", "82I"), ("6Y7L", "QMG")):
        struct = gemmi.read_structure(str(brd / "raw" / f"{name}.pdb"))
        struct.merge_chain_parts()
        ch = struct[0]["A"]
        xyz = [
            list(a.pos)
            for r in ch
            if r.name == ligand
            for a in r
            if not a.element.is_hydrogen and a.altloc in ("\0", "A")
        ]
        center = np.mean(xyz, axis=0)
        nearby = {
            (r.seqid.num - (116 if name == "6Y7L" else 0), r.name): a.pos
            for r in ch
            if r.het_flag == "A"
            for a in r
            if a.name == "CA"
            and a.altloc in ("\0", "A")
            and np.linalg.norm(np.array(list(a.pos)) - center) < 14
        }
        protein[name] = nearby, center
    moving = protein["6YQW"][0]
    fixed, center = protein["6Y7L"]
    keys = sorted(fixed.keys() & moving.keys())
    fit = gemmi.superpose_positions([fixed[i] for i in keys], [moving[i] for i in keys])
    assert len(keys) == 37 and 0.391 < fit.rmsd < 0.393
    ref = reference_graph(brd / "raw/6YQW.pdb", "82I", "A")
    conf = ref.GetConformer()
    for i, pos in enumerate(conf.GetPositions()):
        conf.SetAtomPosition(i, tuple(fit.transform.apply(gemmi.Position(*pos))))
    meta = json.loads((brd / "metadata.json").read_text())["6Y7L"]
    np.testing.assert_allclose(center, meta["center"], atol=1e-8, rtol=0)
    return (
        brd / "6Y7L/receptor.pdbqt",
        path,
        sidecar,
        meta["center"],
        24,
        mol,
        serials,
        ref,
    )


def inputs(aurora, brd, out):
    systems = {"82I": brd_input(brd, out)}
    for name in ("at9283", "danusertib"):
        path = aurora / "prepared" / (name + ".pdbqt")
        graph = (
            aurora
            / "raw"
            / (
                "at9283_chemical.bonds"
                if name == "at9283"
                else "danusertib_native.bonds"
            )
        )
        receptor = aurora / replay.PROTEINS["1OL7_released"]
        center = json.loads((aurora / "prepared/2X6E.json").read_text())["center"]
        ref, _ = audit.refmol(aurora, name)
        mol, serials = audit.frame.read_query(aurora, name)
        systems[name] = (receptor, path, graph, center, 26, mol, serials, ref)
    return systems


def measure(path, ligand, mol, serials, ref):
    # Independent coordinate check against prepared heavy bonds on the original
    # PDBQT precision scale; only existing graph edges are used in this test.
    template = {
        int(l[6:11]): np.array([float(l[k : k + 8]) for k in (30, 38, 46)])
        for l in ligand.read_text().splitlines()
        if l.startswith(("ATOM  ", "HETATM"))
    }
    bonds = [
        (serials[b.GetBeginAtomIdx()], serials[b.GetEndAtomIdx()])
        for b in mol.GetBonds()
    ]
    scores = [
        float(l.split()[3])
        for l in path.read_text().splitlines()
        if l.startswith("REMARK VINA RESULT:")
    ]
    result = []
    for mode, score in zip(positions(path), scores, strict=True):
        for i, s in enumerate(serials):
            mol.GetConformer().SetAtomPosition(i, mode[s])
        err = rdMolAlign.CalcRMS(mol, ref, symmetrizeConjugatedTerminalGroups=False)
        drift = max(
            abs(
                np.linalg.norm(np.subtract(mode[a], mode[b]))
                - np.linalg.norm(template[a] - template[b])
            )
            for a, b in bonds
        )
        assert np.isfinite(err) and drift < 0.004
        result.append((err, score, drift))
    return result


def run(aurora, brd, out, pro, official):
    out.mkdir(exist_ok=False)
    assert replay.digest(official) == OFFICIAL
    (out / "pools").mkdir()
    systems = inputs(aurora, brd, out)
    cells, modes, manifests = [], [], []
    for name, (rec, ligand, graph, center, size, mol, serials, ref) in systems.items():
        for seed in (*SEEDS, 61231) if name == "at9283" else SEEDS:
            arms = (
                ("official", "pro_double")
                if seed % 4 == 1
                else ("pro_double", "official")
            )
            for arm in arms:
                raw = out / "pools" / f"{name}_{seed}_{arm}.pdbqt"
                cap = 25000 if arm == "official" else 50000
                exe = official if arm == "official" else pro
                cmd = [
                    exe,
                    "--receptor",
                    rec,
                    "--ligand",
                    ligand,
                    "--seed",
                    seed,
                    "--exhaustiveness",
                    2,
                    "--cpu",
                    1,
                    "--num_modes",
                    9,
                    "--energy_range",
                    100,
                    "--max_evals",
                    cap,
                    "--out",
                    raw,
                ]
                for k, c in zip("xyz", center, strict=True):
                    cmd.extend([f"--center_{k}", c, f"--size_{k}", size])
                if arm == "pro_double":
                    cmd.extend(
                        [
                            "--ligand_bonds",
                            graph,
                            "--pose_selection",
                            "late",
                            "--search_counters",
                        ]
                    )
                log = raw.with_suffix(".log")
                costs = replay.invoke(cmd, log)
                counts = re.findall(
                    r"\b(?:hunt|second|explicit)_evals=(\d+)", log.read_text()
                )
                result = measure(raw, ligand, mol, serials, ref)
                rms = [r[0] for r in result]
                assert result
                key = dict(ligand=name, seed=seed, arm=arm)
                modes.extend(
                    dict(**key, rank=i, rmsd_A=r, score=s, max_bond_drift_A=d)
                    for i, (r, s, d) in enumerate(result, 1)
                )
                cells.append(
                    dict(
                        **key,
                        n=len(rms),
                        top1=rms[0],
                        best5=min(rms[:5]),
                        oracle=min(rms),
                        regret=rms[0] - min(rms),
                        evals=sum(map(int, counts)) if counts else "",
                        **costs,
                        pool_sha256=replay.digest(raw),
                    )
                )
                manifests.append(
                    dict(
                        **key,
                        receptor=replay.digest(rec),
                        ligand_pdbqt=replay.digest(ligand),
                        bonds=replay.digest(graph) if arm == "pro_double" else None,
                        exe=replay.digest(exe),
                    )
                )
    replay.geom.write_csv(out, "cells.csv", cells)
    replay.geom.write_csv(out, "modes.csv", modes)
    (out / "manifest.json").write_text(json.dumps(manifests, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("aurora", "brd", "out", "pro", "official"):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    run(*(v.resolve() for v in vars(a).values()))

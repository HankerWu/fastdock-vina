"""Read FKBP12 output chemistry, actual search work and the two deposited ligands."""

import argparse
import csv
import hashlib
import json
import re
import sys
import tarfile
from pathlib import Path

import numpy as np
from rdkit import Chem
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "benchmarks/paired_crossdocking"))
from run import alignment


def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def write_csv(path, entries):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=entries[0], lineterminator="\n")
        w.writeheader()
        w.writerows(entries)


def xyz(block, serials):
    atoms = {
        int(l[6:11]): [float(l[k : k + 8]) for k in (30, 38, 46)]
        for l in block.splitlines()
        if l.startswith(("ATOM  ", "HETATM"))
    }
    return np.array([atoms[s] for s in serials])


def molecule(text):
    smiles = next(
        l.split()[2] for l in text.splitlines() if l.startswith("REMARK SMILES ")
    )
    mol = Chem.MolFromSmiles(smiles)
    serial = {}
    for line in text.splitlines():
        if line.startswith("REMARK SMILES IDX"):
            vals = list(map(int, line.split()[3:]))
            serial.update(zip(vals[::2], vals[1::2], strict=True))
    serials = [serial[i + 1] for i in range(mol.GetNumAtoms())]
    bonds = np.array([(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()])
    return mol, serials, bonds


def stereo(mol, coords):
    m = Chem.Mol(mol)
    c = Chem.Conformer(m.GetNumAtoms())
    for i, p in enumerate(coords):
        c.SetAtomPosition(i, p)
    m.AddConformer(c)
    Chem.AssignStereochemistryFrom3D(m, replaceExistingTags=True)
    return Chem.MolToSmiles(m)


def main(work, out):
    out.mkdir(parents=True, exist_ok=True)
    bench = ROOT / "benchmarks/cubic_grid"
    with tarfile.open(bench / "human_panel/prepared.tar.xz") as t:
        refs = {
            p: json.load(t.extractfile(f"{p}/preparation.json"))
            for p in ("1FKG", "1FKH")
        }
    ligand_refs = {d["canonical_smiles"]: d for d in refs.values()}
    with tarfile.open(bench / "human_panel/affinity-pools.tar.xz") as t:
        ligands = {
            i: t.extractfile(f"affinity/{i}.pdbqt").read().decode()
            for i in range(4, 12)
        }
    history = ROOT / "benchmarks/accepted_bfgs/curvature/fkbp/results/runs.csv"
    historic = {
        (r["receiver"], r["compound"], r["seed"]): r["sha256"]
        for r in rows(history)
        if r["arm"] == "fixed"
    }
    detailed, measured, panels, descriptors = [], [], [], []
    runs = rows(work / "runs.csv")
    equal_baselines = n_stereo = n_total = 0
    max_bond_change = 0.0
    for cid, text in ligands.items():
        mol, serial, bonds = molecule(text)
        typed = Chem.MolToSmiles(mol)
        initial = xyz(text, serial)
        assert typed == stereo(mol, initial)
        start_lengths = np.linalg.norm(
            initial[bonds[:, 0]] - initial[bonds[:, 1]], axis=1
        )
        nrot = int(
            next(l.split()[1] for l in text.splitlines() if l.startswith("TORSDOF"))
        )
        descriptors.append(
            dict(
                compound=cid,
                heavy_atoms=len(serial),
                charge=Chem.GetFormalCharge(mol),
                torsdof=nrot,
            )
        )
        for run in (r for r in runs if int(r["compound"]) == cid):
            arm, rec, seed = (run[k] for k in ("arm", "receiver", "seed"))
            name = f"{cid}_{rec}_{seed}_ad4"
            stem = work / arm / "pools" / name
            pdb = stem.with_suffix(".pdbqt").read_bytes()
            digest = hashlib.sha256(pdb).hexdigest()
            assert digest == run["sha256"]
            if arm == "parent":
                assert digest == historic[rec, str(cid), seed]
                equal_baselines += 1
            log = stem.with_suffix(".log").read_text()
            iters = [
                dict(re.findall(r"(\w+)=(\d+)", l))
                for l in log.splitlines()
                if l.startswith("Search trajectory:")
            ]
            assert len(iters) == 3 and "explicit_evals=0" in log
            measured.append(
                dict(
                    compound=cid,
                    receiver=rec,
                    seed=seed,
                    arm=arm,
                    **{
                        k: sum(int(i[k]) for i in iters)
                        for k in ("steps", "hunt_evals", "second_evals")
                    },
                )
            )
            ref = ligand_refs.get(typed)
            if ref:
                # The preparation orders crystal coordinates by its canonical SMILES.
                assert ref["canonical_smiles"] == next(
                    l.split()[2]
                    for l in text.splitlines()
                    if l.startswith("REMARK SMILES ")
                )
                rot, shift, rms_ca, n_ca = alignment(ref, refs[rec])
                targets = np.array(ref["mapped"]) @ rot + shift
            errors = []
            for rank, block in enumerate(pdb.decode().strip().split("ENDMDL")[:-1], 1):
                positions = xyz(block, serial)
                lengths = np.linalg.norm(
                    positions[bonds[:, 0]] - positions[bonds[:, 1]], axis=1
                )
                max_bond_change = max(
                    max_bond_change, float(np.max(np.abs(lengths - start_lengths)))
                )
                n_stereo += typed != stereo(mol, positions)
                n_total += 1
                if ref:
                    err = float(
                        np.sqrt(
                            np.sum((positions - targets) ** 2, axis=2).mean(axis=1)
                        ).min()
                    )
                    detailed.append(
                        dict(
                            compound=cid,
                            crystal=ref["pdb"],
                            receiver=rec,
                            seed=seed,
                            arm=arm,
                            rank=rank,
                            rmsd_A=err,
                            n_ca=n_ca,
                            ca_rmsd_A=rms_ca,
                        )
                    )
                    errors.append(err)
            if ref:
                panels.append(
                    dict(
                        compound=cid,
                        crystal=ref["pdb"],
                        receiver=rec,
                        seed=seed,
                        arm=arm,
                        top1=errors[0],
                        best5=min(errors[:5]),
                        oracle=min(errors),
                        regret=errors[0] - min(errors),
                        severe=int(errors[0] >= 8),
                    )
                )
    assert (
        equal_baselines == 48
        and n_total == sum(int(r["modes"]) for r in runs)
        and len(panels) == 24
    )
    write_csv(out / "geometry_modes.csv", detailed)
    write_csv(out / "geometry_cells.csv", panels)
    distributions = []
    for arm in ("parent", "fixed"):
        selected = [r for r in panels if r["arm"] == arm]
        values = [r["top1"] for r in selected]
        qs = np.quantile(values, [0.25, 0.5, 0.75, 0.9, 0.95])
        distributions.append(
            dict(
                arm=arm,
                cells=len(selected),
                top1_mean=np.mean(values),
                **dict(zip(("q25", "median", "q75", "p90", "p95"), qs, strict=True)),
                **{
                    k: np.mean([r[k] for r in selected])
                    for k in ("best5", "oracle", "regret")
                },
                severe=sum(r["severe"] for r in selected),
            )
        )
    write_csv(out / "geometry_distributions.csv", distributions)
    write_csv(out / "counts.csv", measured)
    scores = {
        (r["arm"], int(r["compound"])): float(r["predicted"])
        for r in rows(work / "scores.csv")
        if r["receiver"] == "receptor_mean"
    }
    for d in descriptors:
        cid = d["compound"]
        d["score_change"] = scores["fixed", cid] - scores["parent", cid]
    write_csv(out / "descriptors.csv", descriptors)
    summary = dict(
        modes=n_total,
        changed_stereochemical_identities=n_stereo,
        max_heavy_bond_change_A=max_bond_change,
        baseline_byte_matches=equal_baselines,
        n_finite_reference_compounds=len(ligand_refs),
        score_change_spearman={
            k: float(
                spearmanr(
                    [d[k] for d in descriptors],
                    [d["score_change"] for d in descriptors],
                ).statistic
            )
            for k in ("heavy_atoms", "torsdof")
        },
    )
    (out / "chemistry.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("out", type=Path)
    args = p.parse_args()
    main(args.work, args.out)

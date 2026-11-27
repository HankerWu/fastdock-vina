"""Energy decomposition at the exact archived first poses and full-potential free minima."""

import argparse
import csv
import lzma
import sys
import tarfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from tempfile import TemporaryDirectory
from time import perf_counter

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold

HERE = Path(__file__).resolve().parent
FROZEN = HERE.parent
BENCH = FROZEN.parent
sys.path.insert(0, str(BENCH.parent / "tools"))
import score_physical_pool as physical


def read(path):
    with lzma.open(path, "rt") if path.suffix == ".xz" else path.open() as stream:
        return list(csv.DictReader(stream))


def write(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def electrostatic(mol, cid):
    # At fixed geometry MMFF's contribution is additive. Verify the difference
    # implementation by evaluating the electrostatic term alone as well.
    props = physical.mmff_properties(mol)

    def value():
        ff = AllChem.MMFFGetMoleculeForceField(
            mol, props, confId=cid, nonBondedThresh=100.0
        )
        return ff.CalcEnergy()

    full = value()
    props.SetMMFFEleTerm(False)
    difference = full - value()
    props.SetMMFFEleTerm(True)
    props.SetMMFFBondTerm(False)
    props.SetMMFFAngleTerm(False)
    props.SetMMFFStretchBendTerm(False)
    props.SetMMFFOopTerm(False)
    props.SetMMFFTorsionTerm(False)
    props.SetMMFFVdWTerm(False)
    direct = value()
    if abs(direct - difference) > 1e-6:
        raise ValueError(f"Nonadditive component: {direct} vs {difference}")
    return direct


def nearest_ion_centers(mol, cid):
    # Post-hoc structural diagnostic in the zwitterionic Factor Xa microstates.
    # Use the central carbons to avoid a choice of charge-localized resonance form.
    carboxylate = mol.GetSubstructMatches(Chem.MolFromSmarts("[O-]-[C;X3]=[O;X1]"))
    amidinium = mol.GetSubstructMatches(Chem.MolFromSmarts("[C;X3](=[N+])-[N]"))
    if not (len(carboxylate) == 1 and len(amidinium) == 2):
        return ""
    xyz = mol.GetConformer(cid).GetPositions()
    ac = carboxylate[0][1]
    return min(float(np.linalg.norm(xyz[ac] - xyz[am[0]])) for am in amidinium)


def free_minimum(ligand):
    # Same generation, full-potential minimization, and selection as gas_reference.
    # Retain the geometry of the winning conformer to evaluate its components.
    mol = Chem.Mol(ligand)
    mol.RemoveAllConformers()
    params = AllChem.ETKDGv3()
    params.randomSeed = physical.GAS_EMBED_SEED
    params.numThreads = 1
    ids = AllChem.EmbedMultipleConfs(
        mol, numConfs=physical.GAS_CONFORMERS, params=params
    )
    props = physical.mmff_properties(mol)
    energies = []
    for cid in ids:
        try:
            energies.append((physical.minimized_energy(mol, props, cid, ()), cid))
        except ValueError:
            continue
    energy, cid = min(energies)
    return energy, electrostatic(mol, cid), len(energies), nearest_ion_centers(mol, cid)


def read_first(path):
    return next(s for s in physical.read_poses(path) if s.native_rank == 1)


def compute(job):
    meta, states, scale, archived_free = job
    ligand, identity = physical.decode_ligand(read_first(states[0][0]))
    query = Chem.MolFromSmiles(meta["smiles"])
    if Chem.MolToSmiles(query) != identity:
        raise ValueError(f"Assay and pose microstates differ for {meta['compound']}")
    charge = Chem.GetFormalCharge(query)
    heavy = query.GetNumHeavyAtoms()
    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=query)
    free, free_elec, nfree, free_centers = free_minimum(ligand)
    rows = []
    for path, receptor, seed, archived_native, archived_full in states:
        pose = read_first(path)
        mol, current = physical.decode_ligand(pose)
        if current != identity or archived_native != pose.native_score:
            raise ValueError(f"Wrong frozen mode: {path}")
        fixed = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
        conf = mol.GetConformer()
        before = conf.GetPositions()[fixed]
        bound = physical.minimized_energy(
            mol, physical.mmff_properties(mol), conf.GetId(), fixed
        )
        if not np.array_equal(before, conf.GetPositions()[fixed]):
            raise ValueError("Heavy coordinate changed")
        bound_elec = electrostatic(mol, conf.GetId())
        base = pose.native_score - scale * (pose.native_intra - pose.native_unbound)
        full = physical.mmff_affinity_score(pose, scale, bound, free)
        delta = scale * (bound_elec - free_elec)
        rows.append(
            {
                "assay": meta["assay"],
                "compound": meta["compound"],
                "receptor": receptor,
                "seed": seed,
                "ki_nm": meta["ki_nm"],
                "heavy": heavy,
                "charge": charge,
                "scaffold": scaffold,
                "scale": scale,
                "free": free,
                "nfree": nfree,
                "free_elec": free_elec,
                "bound": bound,
                "bound_elec": bound_elec,
                "free_centers_A": free_centers,
                "bound_centers_A": nearest_ion_centers(mol, conf.GetId()),
                "vina_intra": pose.native_intra,
                "vina_unbound": pose.native_unbound,
                "free_replay_error": free - archived_free,
                "full_replay_error": full - archived_full,
                "native": pose.native_score,
                "full": full,
                "without_elec": full - delta,
                "elec_only": base + delta,
            }
        )
    return rows


def jobs(work):
    for prefix, folder in [
        ("native/outputs.tar.xz", "bace"),
        ("affinity/fxa_frozen.tar.xz", "fxa"),
    ]:
        with tarfile.open(FROZEN / prefix) as archive:
            archive.extractall(work / folder, filter="data")
    strain = FROZEN / "strain"
    scales = {
        r["compound"]: float(r["flexibility_scale"])
        for r in read(strain / "affinity_modes.csv.xz")
    }
    free = {
        r["compound"]: float(r["free_8"]) for r in read(strain / "free_references.csv")
    }
    labels = {
        r["compound"]: {
            "assay": r["assay"],
            "compound": r["compound"],
            "smiles": r["smiles"],
            "ki_nm": r["value_nM"],
        }
        for r in read(BENCH / "geometry_comparison/affinity/panel.csv")
        if r["target"] == "BACE1"
    }
    states = defaultdict(list)
    for r in read(strain / "fresh_mmff_modes.csv.xz"):
        if r["rank"] != "1" or r["target"] != "BACE1":
            continue
        c, rec, seed = r["compound"], r["receptor"], r["seed"]
        states[c].append(
            (
                work / "bace/pools" / f"{rec}_{c}_{seed}.pdbqt",
                rec,
                seed,
                float(r["native"]),
                float(r["unit_8"]),
            )
        )
    for r in read(BENCH / "affinity_partition/fxa_series.csv"):
        if r["included"] == "1":
            c = r["ligand"]
            labels[c] = {
                "assay": "CHEMBL661804",
                "compound": c,
                "smiles": r["docked_smiles"],
                "ki_nm": r["ki_nm"],
            }
    for r in read(strain / "fxa_modes.csv.xz"):
        if r["rank"] != "1":
            continue
        c, rec, seed = r["compound"], r["receptor"], r["seed"]
        states[c].append(
            (
                work / "fxa" / rec / f"{c}_{seed}.pdbqt",
                rec,
                seed,
                float(r["native"]),
                float(r["full_unit"]),
            )
        )
        scales[c] = float(r["flexibility_scale"])
        free[c] = float(r["free_mmff"])
    if set(states) != set(labels) or any(
        len({s[1:3] for s in v}) != 4 or len(v) != 4 for v in states.values()
    ):
        raise ValueError("Unpaired receptor/seed cell")
    return [(labels[c], states[c], scales[c], free[c]) for c in sorted(states)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    parser.add_argument("--jobs", type=int, default=1)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    start = perf_counter()
    with TemporaryDirectory() as tmp, ProcessPoolExecutor(args.jobs) as pool:
        rows = [row for group in pool.map(compute, jobs(Path(tmp))) for row in group]
    write(args.out / "scores.csv", rows)
    print(f"{len(rows)} fixed native first poses, {perf_counter() - start:.2f}s")
    for field in ["free_replay_error", "full_replay_error"]:
        err = max(abs(r[field]) for r in rows)
        print(f"Maximum absolute {field}: {err:.9g}")
        if err > 1e-4:
            raise ValueError("Reference replay differs from archived ledger")


if __name__ == "__main__":
    main()

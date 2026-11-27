"""Frozen ligand-only OBC-II single points and same-assay counterfactuals."""

import argparse
import hashlib
import itertools
import json
import sys
import tarfile
from collections import defaultdict
from pathlib import Path
from statistics import median
from tempfile import TemporaryDirectory
from time import perf_counter
from types import SimpleNamespace

import numpy as np
from openmm import app
from openmm.app.internal.customgbforces import GBSAOBC2Force
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.stats import kendalltau, spearmanr

FROZEN = Path(__file__).resolve().parents[1]
sys.path[:0] = [
    str(FROZEN / "strain_components"),
    str(FROZEN / "transfer/galectin_mmff"),
    str(FROZEN.parents[1] / "tools"),
]
import decompose
import polar_endpoint
import replay as galectin
import score_physical_pool as physical
from summarize import assignments, pair_rows, sensitivity, write

ARMS = ("native", "full", "gb_fixed", "gb_reselected")


def ligand_context(mol, platform):
    topology = app.Topology()
    residue = topology.addResidue("LIG", topology.addChain())
    atoms = [
        topology.addAtom(
            a.GetSymbol(), app.Element.getByAtomicNumber(a.GetAtomicNum()), residue
        )
        for a in mol.GetAtoms()
    ]
    for b in mol.GetBonds():
        topology.addBond(atoms[b.GetBeginAtomIdx()], atoms[b.GetEndAtomIdx()])
    radii_nm, scales = np.transpose(GBSAOBC2Force.getStandardParameters(topology))
    props = physical.mmff_properties(mol)
    q = [props.GetMMFFPartialCharge(i) for i in range(len(atoms))]
    if abs(sum(q) - Chem.GetFormalCharge(mol)) > 1e-6:
        raise ValueError("MMFF charges do not conserve the assigned charge")
    params = SimpleNamespace(charges=q, radii=radii_nm * 10, scales=scales)
    return polar_endpoint.PolarContext([params], platform, solute_dielectric=1)


def free_conformers(ligand, context):
    mol = Chem.Mol(ligand)
    mol.RemoveAllConformers()
    params = AllChem.ETKDGv3()
    params.randomSeed = physical.GAS_EMBED_SEED
    params.numThreads = 1
    ids = AllChem.EmbedMultipleConfs(
        mol, numConfs=physical.GAS_CONFORMERS, params=params
    )
    props = physical.mmff_properties(mol)
    rows = []
    for cid in ids:
        try:
            energy = physical.minimized_energy(mol, props, cid, ())
        except ValueError:
            continue
        gb = context.energy(mol.GetConformer(cid).GetPositions())
        rows.append(
            {
                "cid": cid,
                "vacuum": energy,
                "gb": gb,
                "solution": energy + gb,
                "ion_centers_A": decompose.nearest_ion_centers(mol, cid),
            }
        )
    return rows


def score(meta, initial, blocks, platform):
    ligand, identity = physical.decode_ligand(initial)
    ctx = ligand_context(ligand, platform)
    conformers = free_conformers(ligand, ctx)
    f = min(conformers, key=lambda x: x["vacuum"])
    fw = min(conformers, key=lambda x: x["solution"])
    for c in conformers:
        c.update(compound=meta["compound"], vacuum_min=c is f, solution_min=c is fw)

    # OBC parameters refer to the original explicit atom order.
    def order(m):
        return [
            (
                a.GetAtomicNum(),
                a.GetFormalCharge(),
                tuple(n.GetIdx() for n in a.GetNeighbors()),
            )
            for a in m.GetAtoms()
        ]

    expected_order = order(ligand)
    outputs = []
    for pose, row in blocks:
        mol, ident = physical.decode_ligand(pose)
        if ident != identity:
            raise ValueError("Mismatched prepared molecular graph")
        if order(mol) != expected_order:
            raise ValueError("Mismatched MMFF/OBC atom order")
        heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
        conf = mol.GetConformer()
        before = conf.GetPositions()[heavy]
        eb = physical.minimized_energy(
            mol, physical.mmff_properties(mol), conf.GetId(), heavy
        )
        xyz = conf.GetPositions()
        if not np.array_equal(before, xyz[heavy]):
            raise ValueError("Frozen heavy coordinates moved")
        scale = float(row["scale"])
        native = pose.native_score
        full = physical.mmff_affinity_score(pose, scale, eb, f["vacuum"])
        if (
            abs(native - float(row["native"])) > 1e-10
            or max(
                abs(full - float(row["full"])), abs(f["vacuum"] - float(row["free"]))
            )
            > 1e-5
        ):
            raise ValueError(f"Archived native or full-MMFF score diverges for {meta}")
        gb = ctx.energy(xyz)
        fixed = full + scale * (gb - f["gb"])
        reselected = fixed + scale * (f["solution"] - fw["solution"])
        direct = physical.mmff_affinity_score(pose, scale, eb + gb, fw["solution"])
        if abs(direct - reselected) > 1e-9:
            raise ValueError("Independent score reconstruction failed")
        outputs.append(
            dict(
                **meta,
                pdb=row["receptor"],
                seed=row["seed"],
                model=pose.model,
                core_rmsd=row.get("core_rmsd", ""),
                charge=Chem.GetFormalCharge(mol),
                scale=scale,
                bound_mmff=eb,
                bound_gb=gb,
                free_vac_id=f["cid"],
                free_solution_id=fw["cid"],
                native=native,
                full=full,
                gb_fixed=fixed,
                gb_reselected=reselected,
            )
        )
    return outputs, conformers


def proteases(work):
    metadata = {
        r["compound"]: r
        for r in decompose.read(FROZEN / "strain_components/scores.csv")
    }
    for meta, states, scale, free in decompose.jobs(work):
        c = meta["compound"]
        cases = [
            (
                decompose.read_first(p),
                {
                    "receptor": rec,
                    "seed": seed,
                    "native": native,
                    "full": full,
                    "free": free,
                    "scale": scale,
                },
            )
            for p, rec, seed, native, full in states
        ]
        yield {k: metadata[c][k] for k in ("compound", "assay")}, cases[0][0], cases


def lectins(work):
    with tarfile.open(galectin.GAL / "evidence.tar.xz") as arc:
        arc.extractall(work, filter="data")
    ledgers = {
        (r["pdb"], int(r["compound"]), int(r["seed"]), int(r["model"])): r
        for r in decompose.read(FROZEN / "transfer/galectin_mmff/modes.csv.xz")
    }
    for c in range(5, 18):
        ligand = work / "prepared" / f"{c}.pdbqt"
        first = physical.Pose(1, 0, 1, ligand.read_text(), 0, 0)
        cases = []
        for pdb, seed in itertools.product(
            sum(galectin.TARGETS.values(), ()), galectin.SEEDS
        ):
            pool = work / "pools" / f"{pdb}_{c}_{seed}_native.pdbqt"
            if (
                hashlib.sha256(pool.read_bytes()).hexdigest()
                != json.loads(pool.with_suffix(".json").read_text())["out_sha256"]
            ):
                raise ValueError("Changed native pool")
            poses = physical.read_poses(pool)
            physical.fixed_coordinate_inputs(poses, ligand, None)
            for pose in poses:
                r = ledgers.pop((pdb, c, seed, pose.model))
                cases.append(
                    (
                        pose,
                        {
                            **r,
                            "receptor": pdb,
                            "seed": seed,
                            "core_rmsd": float(r["core_rmsd"]),
                        },
                    )
                )
        yield {"compound": c, "assay": "lectin"}, first, cases
    if ledgers:
        raise ValueError("Incomplete lectin pools")


def summarize_proteases(rows, out):
    metadata = {
        r["compound"]: r
        for r in decompose.read(FROZEN / "strain_components/scores.csv")
    }
    by_rec = defaultdict(lambda: defaultdict(dict))
    for r in rows:
        for rec in [r["pdb"], "both"] if r["assay"] == "CHEMBL661804" else [r["pdb"]]:
            by_rec[r["assay"], rec][r["compound"]][r["pdb"], r["seed"]] = r
    pairs, summary, correlations, influence = [], [], [], []
    for (assay, rec), states in sorted(by_rec.items()):
        edges = list(pair_rows(assay, rec, states, metadata, ARMS))
        pairs.extend(edges)
        subsets = defaultdict(list)
        for e in edges:
            a, b = e["stronger"], e["weaker"]
            for group in assignments(metadata[a], metadata[b]):
                subsets[group].append(e)
            if e["fold"] >= 10:
                subsets["tenfold"].append(e)
            if assay == "CHEMBL661804":
                subsets[
                    "distorted_incident" if "CHEMBL369480" in (a, b) else "others"
                ].append(e)
        for subset, group in sorted(subsets.items()):
            summary.append(
                dict(
                    assay=assay,
                    receptor=rec,
                    subset=subset,
                    pairs=len(group),
                    **{arm: sum(e[arm] for e in group) for arm in ARMS},
                )
            )
        for arm in ARMS[1:]:
            low, high = sensitivity(edges, states, arm, "native")
            influence.append(
                {
                    "assay": assay,
                    "receptor": rec,
                    "arm": arm,
                    "low_pp": low,
                    "high_pp": high,
                }
            )
        if rec != "both":
            kd = [float(metadata[c]["ki_nm"]) for c in states]
            for arm in ARMS:
                scores = [
                    median(r[arm] for r in blocks.values())
                    for blocks in states.values()
                ]
                correlations.append(
                    {
                        "assay": assay,
                        "receptor": rec,
                        "arm": arm,
                        "spearman": spearmanr(scores, kd).statistic,
                        "kendall": kendalltau(scores, kd).statistic,
                    }
                )
    write(out / "protease_pairs.csv.xz", pairs)
    write(out / "protease_summary.csv", summary)
    write(out / "protease_correlations.csv", correlations)
    write(out / "protease_loo.csv", influence)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("out", type=Path)
    p.add_argument("--platform", choices=("Reference", "CPU"), default="Reference")
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    modes, refs = [], []
    with TemporaryDirectory() as tmp:
        work = Path(tmp)
        for loader in (proteases, lectins):
            start = perf_counter()
            results = []
            for meta, first, cases in loader(work):
                m, f = score(meta, first, cases, args.platform)
                results.extend(m)
                refs.extend(f)
            print(
                f"{loader.__name__}: {len(results)} fixed modes in {perf_counter() - start:.2f}s",
                flush=True,
            )
            if loader is proteases:
                summarize_proteases(results, args.out)
            else:
                galectin.evaluate(results, work, args.out, ARMS, write_table=write)
                pools = defaultdict(list)
                for row in results:
                    pools[row["pdb"], row["compound"], row["seed"]].append(row)
                for pool in pools.values():
                    orders = [
                        [
                            r["model"]
                            for r in sorted(pool, key=lambda r: (r[arm], r["model"]))
                        ]
                        for arm in ("gb_fixed", "gb_reselected")
                    ]
                    if orders[0] != orders[1]:
                        raise ValueError(
                            "Changing free references modified a pose ordering"
                        )
            modes.extend(results)
    write(args.out / "modes.csv.xz", modes)
    write(args.out / "references.csv.xz", refs)


if __name__ == "__main__":
    main()

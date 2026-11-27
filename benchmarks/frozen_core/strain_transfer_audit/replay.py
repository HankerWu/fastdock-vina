"""Unit-weight MMFF score diagnostics on archived input, bound and free geometries."""

import argparse
import itertools
import math
import sys
import tarfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from decimal import Decimal
from pathlib import Path
from statistics import median
from tempfile import TemporaryDirectory
from time import perf_counter

from scipy.stats import kendalltau, spearmanr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "strain_components"))
sys.path.insert(0, str(HERE.parents[1] / "strain_transfer"))

import decompose as previous
from alkylbenzene_diagnostic import TERMS, term_values
from evaluate import fixed_energy, physical
from rdkit import Chem
from rdkit.Chem import AllChem
from summarize import assignments, write

LOCAL = ("Bond", "Angle", "StretchBend", "Oop")
ARMS = ("native", "full", "without_local", "local_only", "input_reference")


def free_geometry(ligand):
    # The production helper returns only a scalar. Retain the very same
    # converged geometry here so all seven single-point terms can be evaluated.
    mol = Chem.Mol(ligand)
    mol.RemoveAllConformers()
    settings = AllChem.ETKDGv3()
    settings.randomSeed = physical.GAS_EMBED_SEED
    settings.numThreads = 1
    props = physical.mmff_properties(mol)
    energies = []
    for cid in AllChem.EmbedMultipleConfs(
        mol, numConfs=physical.GAS_CONFORMERS, params=settings
    ):
        try:
            energies.append((physical.minimized_energy(mol, props, cid, ()), cid))
        except ValueError:
            continue
    energy, cid = min(energies)
    return energy, checked_terms(mol, cid, energy)


def checked_terms(mol, cid, energy):
    terms = term_values(mol, cid)
    rebuilt = AllChem.MMFFGetMoleculeForceField(
        mol, physical.mmff_properties(mol), confId=cid
    ).CalcEnergy()
    if not math.isclose(sum(terms.values()), rebuilt, abs_tol=1e-6, rel_tol=0):
        raise ValueError(
            f"Seven single-point terms do not reconstruct the full single point: {rebuilt} vs {sum(terms.values())}"
        )
    return terms


def fixed_geometry(pose):
    mol, identity = physical.decode_ligand(pose)
    energy = fixed_energy(mol, physical.mmff_properties(mol))
    return energy, checked_terms(mol, mol.GetConformer().GetId(), energy), identity


def components(job):
    (meta, states, scale, archived_free), prepared_path = job
    prepared_text = prepared_path.read_text()
    prep = physical.Pose(0, 0, 0, prepared_text, None, None)
    input_energy, input_terms, identity = fixed_geometry(prep)
    ligand, _ = physical.decode_ligand(prep)
    free, free_terms = free_geometry(ligand)
    if not math.isclose(free, archived_free, abs_tol=1e-6, rel_tol=0):
        raise ValueError(
            "Independent reference replay changed the selected free energy"
        )
    rows = []
    for path, receptor, seed, native, full in states:
        pose = previous.read_first(path)
        physical.fixed_coordinate_inputs([pose], prepared_path, None)
        bound, bound_terms, current = fixed_geometry(pose)
        if current != identity or native != pose.native_score:
            raise ValueError("Native pose or input identity mismatch")
        base = native - scale * (pose.native_intra - pose.native_unbound)
        delta = {t: bound_terms[t] - free_terms[t] for t in TERMS}
        local_delta = scale * sum(delta[t] for t in LOCAL)
        if not math.isclose(
            base + scale * (bound - free), full, abs_tol=1e-6, rel_tol=0
        ):
            raise ValueError("Independent full score differs from published full score")
        rows.append(
            dict(
                assay=meta["assay"],
                compound=meta["compound"],
                receptor=receptor,
                seed=seed,
                ki_nm=meta["ki_nm"],
                scale=scale,
                free=free,
                prepared=input_energy,
                bound=bound,
                bound_rebuild_delta=sum(bound_terms.values()) - bound,
                free_rebuild_delta=sum(free_terms.values()) - free,
                input_rebuild_delta=sum(input_terms.values()) - input_energy,
                **{"bound_minus_free_" + t: delta[t] for t in TERMS},
                **{
                    "prepared_minus_free_" + t: input_terms[t] - free_terms[t]
                    for t in TERMS
                },
                native=native,
                full=full,
                without_local=full - local_delta,
                local_only=base + local_delta,
                input_reference=base + scale * (bound - input_energy),
            )
        )
    return rows


def tabulate(rows, out):
    labels = {
        r["compound"]: r
        for r in previous.read(HERE.parent / "strain_components/scores.csv")
    }
    data = defaultdict(lambda: defaultdict(dict))
    for row in rows:
        assay, c, rec, seed = (
            row[k] for k in ("assay", "compound", "receptor", "seed")
        )
        for block in (rec, "both") if assay == "CHEMBL661804" else (rec,):
            data[assay, block][c][rec, seed] = row
    pairs, totals, corrs = [], [], []
    for (assay, block), group in sorted(data.items()):
        subsets = defaultdict(list)
        compounds = sorted(group, key=lambda c: Decimal(labels[c]["ki_nm"]))
        for a, b in itertools.combinations(compounds, 2):
            fold = Decimal(labels[b]["ki_nm"]) / Decimal(labels[a]["ki_nm"])
            if fold < 3:
                continue
            if group[a].keys() != group[b].keys():
                raise ValueError("Unbalanced permutation blocks")
            pair = {
                "assay": assay,
                "receptor": block,
                "stronger": a,
                "weaker": b,
                "fold": float(fold),
            }
            for arm in ARMS:
                if assay == "CHEMBL661804":
                    gap = median(
                        group[b][k][arm] - row[arm] for k, row in group[a].items()
                    )
                else:
                    gap = median(r[arm] for r in group[b].values()) - median(
                        r[arm] for r in group[a].values()
                    )
                pair[arm + "_gap"] = gap
                pair[arm] = 0.5 if abs(gap) < 1e-10 else int(gap > 0)
            pairs.append(pair)
            for subset in assignments(labels[a], labels[b]):
                subsets[subset].append(pair)
            if fold >= 10:
                subsets["tenfold"].append(pair)
        for subset, records in sorted(subsets.items()):
            totals.append(
                dict(
                    assay=assay,
                    receptor=block,
                    subset=subset,
                    pairs=len(records),
                    **{arm: sum(r[arm] for r in records) for arm in ARMS},
                )
            )
        if block != "both":
            ki = [float(labels[c]["ki_nm"]) for c in compounds]
            for arm in ARMS:
                score = [median(r[arm] for r in group[c].values()) for c in compounds]
                corrs.append(
                    {
                        "assay": assay,
                        "receptor": block,
                        "arm": arm,
                        "ligands": len(ki),
                        "spearman": spearmanr(score, ki).statistic,
                        "kendall_tau_b": kendalltau(score, ki).statistic,
                    }
                )
    baseline = previous.read(HERE.parent / "strain_components/summary.csv")
    key = lambda r: (r["assay"], r["receptor"], r["subset"])

    def original_counts(records):
        return {
            key(r): tuple(float(r[f]) for f in ("pairs", "native", "full"))
            for r in records
        }

    if original_counts(totals) != original_counts(baseline):
        raise ValueError("Historical pair counts do not reproduce")
    write(out / "pairs.csv.xz", pairs)
    write(out / "summary.csv", totals)
    write(out / "correlations.csv", corrs)


def check_panel(rows, out):
    """Independently eliminate the gas reference using two published scores."""
    strain = HERE.parent / "strain"
    firsts = [
        r for r in previous.read(strain / "fresh_mmff_modes.csv.xz") if r["rank"] == "1"
    ]
    key_fields = ("assay", "receptor", "compound", "seed")
    direct = {
        tuple(r[k] for k in key_fields): r["input_reference"]
        for r in rows
        if r["assay"] != "CHEMBL661804"
    }
    by_ligand = defaultdict(list)
    checked = set()
    for r in firsts:
        key = tuple(r[k] for k in key_fields)
        score = float(r["unit_8"]) - float(r["prepared_prior_8"]) + float(r["native"])
        if key in direct:
            if abs(direct[key] - score) > 1e-6:
                raise ValueError(
                    "Coordinate replay differs from prepared-prior identity"
                )
            checked.add(key)
        by_ligand[key[:3]].append(score)
    if checked != direct.keys() or any(len(v) != 2 for v in by_ligand.values()):
        raise ValueError("Incomplete native first-pose reference replay")
    totals = defaultdict(
        lambda: {"pairs": 0, "native": 0, "full": 0, "input_reference": 0}
    )
    for r in previous.read(strain / "fresh_mmff_pairs.csv.xz"):
        if r["policy"] != "first" or r["arm"] not in ("native", "unit_8"):
            continue
        cell = totals[r["target"], r["assay"], r["receptor"]]
        if r["arm"] == "unit_8":
            cell["full"] += float(r["credit"])
            continue
        group = r["assay"], r["receptor"]
        gap = median(by_ligand[*group, r["weaker"]]) - median(
            by_ligand[*group, r["stronger"]]
        )
        cell["pairs"] += 1
        cell["native"] += float(r["credit"])
        cell["input_reference"] += 0.5 if abs(gap) < 1e-10 else int(gap > 0)
    write(
        out / "panel.csv",
        [
            dict(target=k[0], assay=k[1], receptor=k[2], **v)
            for k, v in sorted(totals.items())
        ],
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("out", type=Path)
    p.add_argument("--jobs", type=int, default=1)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    start = perf_counter()
    with TemporaryDirectory() as tmp, ProcessPoolExecutor(args.jobs) as pool:
        work = Path(tmp)
        with tarfile.open(
            previous.BENCH / "geometry_comparison/affinity/inputs.tar.xz"
        ) as archive:
            archive.extractall(work, filter="data")
        tasks = []
        for job in previous.jobs(work):
            meta = job[0]
            name = meta["compound"] + ".pdbqt"
            if meta["assay"] == "CHEMBL661804":
                path = work / "fxa/1FJS" / name
                if path.read_bytes() != (work / "fxa/2P16" / name).read_bytes():
                    raise ValueError("Factor Xa receptors used different ligand inputs")
            else:
                path = work / "prepared" / name
            tasks.append((job, path))
        rows = [r for part in pool.map(components, tasks) for r in part]
    write(args.out / "components.csv", rows)
    tabulate(rows, args.out)
    check_panel(rows, args.out)
    print(
        f"{len(tasks)} molecules, {len(rows)} frozen poses, {perf_counter() - start:.2f} s"
    )


if __name__ == "__main__":
    main()

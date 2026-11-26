"""Probe a protein-derived zinc vacancy against cognate and cross-docked hCA II pools."""

import argparse
import csv
import itertools
import json
import math
import re
import subprocess
import sys
from functools import cache, partial
from pathlib import Path

import numpy as np
from rdkit import Chem

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "tools"))
from select_zinc_pose import (
    COORDINATION_DISTANCES,
    assess_pose,
    assess_site,
    parse_atoms,
    read_models,
    read_score_header,
    score_sulfonamide_pose,
)


@cache
def ligand_mapping(prepared, ligand):
    """Return graph, heavy-atom PDBQT serials and the unique phenylsulfonamide head."""
    mol = Chem.SDMolSupplier(str(prepared / f"{ligand}_graph.sdf"))[0]
    query = Chem.MolFromSmarts("[N-]S(=O)(=O)c1ccccc1")
    matches = mol.GetSubstructMatches(query)
    if len(matches) != 1:
        raise ValueError(f"{ligand}: expected one deprotonated phenylsulfonamide")
    mapping = {}
    for line in (prepared / f"{ligand}_ligand.pdbqt").read_text().splitlines():
        if line.startswith("REMARK INDEX MAP "):
            pairs = list(map(int, line.split()[3:]))
            mapping.update(zip(pairs[::2], pairs[1::2], strict=True))
    serials = tuple(mapping[i + 1] for i in range(mol.GetNumAtoms()))
    if len(set(serials)) != len(serials):
        raise ValueError(f"{ligand}: duplicate mapped serial")
    return mol, serials, matches[0]


@cache
def receptor_site(prepared, receptor):
    lines = (prepared / f"{receptor}_receptor.pdbqt").read_text().splitlines()
    atoms = parse_atoms(lines)
    (zinc,) = (s for s, a in atoms.items() if a.atom_type == "Zn")
    site = assess_site(atoms, zinc)
    if not site["applicable"]:
        raise ValueError(f"{receptor}: {site['reason']}")
    return site


def load_anchor(prepared, receptor, ligand):
    """Construct the donor serial and vacancy in the receptor frame without ligand labels."""
    site = receptor_site(prepared, receptor)
    _, serials, head = ligand_mapping(prepared, ligand)
    ideal = COORDINATION_DISTANCES["N"][0]
    point = np.array(site["zinc_xyz"]) + ideal * np.array(site["vacant_direction"])
    return serials[head[0]], point


@cache
def symmetry_mappings(prepared, ligand):
    mol, _, _ = ligand_mapping(prepared, ligand)
    permutations = list(mol.GetSubstructMatches(mol, uniquify=False, maxMatches=100000))
    if len(permutations) >= 100000:
        raise ValueError("Graph symmetry enumeration exhausted")
    # RDKit's receptor-frame metric permits carboxylate oxygen exchanges.
    for _, left, right in mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)[O-]")):
        exchanged = []
        for permutation in permutations:
            new = list(permutation)
            new[left], new[right] = new[right], new[left]
            exchanged.append(tuple(new))
        permutations.extend(exchanged)
    return np.asarray(sorted(set(permutations)))


def error_components(prepared, ligand, coordinates, reference):
    """Use one full-graph mapping for the whole pose and both disjoint fragments."""
    _, serials, head = ligand_mapping(prepared, ligand)
    query = np.array([coordinates[s].xyz for s in serials])
    ref = np.array(reference)
    permutations = symmetry_mappings(prepared, ligand)
    squared = np.square(query[permutations] - ref).sum(axis=2)
    errors = squared[np.argmin(squared.sum(axis=1))]
    head_mask = np.zeros(len(errors), dtype=bool)
    head_mask[list(head)] = True
    return {
        "rmsd": float(np.sqrt(errors.mean())),
        "head_rmsd": float(np.sqrt(errors[head_mask].mean())),
        "tail_rmsd": float(np.sqrt(errors[~head_mask].mean())),
    }


@cache
def thr_amide_position(prepared, receptor):
    lines = (prepared / f"{receptor}_receptor.pdbqt").read_text().splitlines()
    (amide,) = [
        line
        for line in lines
        if line.startswith("ATOM")
        and line[12:16].strip() == "N"
        and line[17:26].split() == ["THR", "A", "199"]
    ]
    return next(iter(parse_atoms([amide]).values())).xyz


def candidate(prepared, rec, lig, lines, native_score, model):
    atoms = parse_atoms(lines)
    _, serials, head = ligand_mapping(prepared, lig)
    donor, _ = load_anchor(prepared, rec, lig)
    site = receptor_site(prepared, rec)
    pose = dict(native_score=native_score, **assess_pose(atoms, donor, site))
    negative_control = dict(
        site, vacant_direction=[-x for x in site["vacant_direction"]]
    )
    reverse = dict(
        native_score=native_score, **assess_pose(atoms, donor, negative_control)
    )
    oxygen_dist = min(
        math.dist(atoms[serials[i]].xyz, thr_amide_position(prepared, rec))
        for i in head[2:4]
    )
    directional = score_sulfonamide_pose(pose, site, "continuous")
    scores = {
        "native": native_score,
        "hard": native_score if pose["passes_geometry"] else None,
        "eligible": native_score if directional["scoring_eligible"] else None,
        "radial": score_sulfonamide_pose(pose, site, "radial")["ranking_score"],
        "directional": directional["ranking_score"],
        "reversed_axis": score_sulfonamide_pose(
            reverse, negative_control, "continuous"
        )["ranking_score"],
    }
    return atoms, dict(
        model=model,
        zinc_N=pose["zinc_distance_angstrom"],
        vacancy_error=pose["virtual_point_error_angstrom"],
        Thr199_N_to_sulfonyl_O=oxygen_dist,
        **scores,
    )


def pool_rows(prepared, pools, metadata):
    for panel, folder in pools:
        for path in sorted(folder.glob("*_native.pdbqt")):
            rec, lig, seed, _ = path.stem.split("_")
            reference = metadata[lig]["references"][rec]["coords"]
            for model, lines in read_models(path.read_text().splitlines()).items():
                native = read_score_header(model, lines)["native_score"]
                atoms, result = candidate(prepared, rec, lig, lines, native, model)
                yield dict(
                    panel=panel,
                    receptor=rec,
                    ligand=lig,
                    seed=int(seed),
                    **result,
                    **error_components(prepared, lig, atoms, reference),
                )


def normalized_raw(binary, work, lig, pose):
    cmd = [
        str(binary),
        "--cpu",
        "1",
        "--seed",
        "42",
        "--score_only",
        "--net_polar",
        "native",
        "--verbosity",
        "0",
        "--config",
        str(work / "prepared" / f"{lig}_box.txt"),
        "--receptor",
        str(work / "prepared" / f"{lig}_receptor.pdbqt"),
        "--ligand",
        str(work / "prepared" / f"{lig}_ligand.pdbqt"),
        "--ligand_bonds",
        str(work / "crystal_start" / f"{lig}.bonds"),
        "--pose_coordinates",
        str(pose),
    ]
    output = subprocess.check_output(cmd, text=True)
    terms = re.search(r"^Net polar native: (.+)$", output, re.MULTILINE)
    return float(terms[1].split()[2])


METHODS = ("native", "hard", "eligible", "radial", "directional", "reversed_axis")


def choices(rows):
    for method in METHODS:
        eligible = [r for r in rows if r[method] is not None]
        best = min(eligible, key=lambda r: (r[method], r["model"]), default=None)
        yield method, eligible, best


def measured_pose(binary, work, lig, reference, path, model):
    prepared = work / "prepared"
    atoms, row = candidate(
        prepared,
        lig,
        lig,
        path.read_text().splitlines(),
        normalized_raw(binary, work, lig, path),
        model,
    )
    return dict(**row, **error_components(prepared, lig, atoms, reference))


def spiked_pools(binary, work, metadata, symmetry_poses, reflected_poses):
    for lig, meta in metadata.items():
        reference = meta["references"][lig]["coords"]
        base = work / "crystal_start"
        measure = partial(measured_pose, binary, work, lig, reference)

        crystal_poses = {
            "projected": measure(base / f"{lig}_start.pdbqt", 0),
            "native_refined": measure(base / f"{lig}_native_refined.pdbqt", 0),
        }
        for directory, prefix in (
            (symmetry_poses, "symmetric"),
            (reflected_poses, "reflected_symmetric"),
        ):
            if directory:
                crystal_poses[prefix + "_projected"] = measure(
                    directory / f"{lig}_start.pdbqt", 0
                )
                crystal_poses[prefix + "_refined"] = measure(
                    directory / f"{lig}_native_refined.pdbqt", 0
                )
        for seed in (42, 1729):
            paths = sorted(base.glob(f"{lig}_{seed}_*.pdbqt"))
            if not paths:
                raise ValueError(f"Missing crystal-start decoys for {lig} {seed}")
            pool = [measure(path, int(path.stem.split("_")[-1])) for path in paths]
            for spike, crystal in crystal_poses.items():
                for method, eligible, best in choices([*pool, crystal]):
                    decoys = [r[method] for r in eligible if r["model"] != 0]
                    yield {
                        "ligand": lig,
                        "seed": seed,
                        "spike": spike,
                        "method": method,
                        "crystal_eligible": crystal[method] is not None,
                        "crystal_rmsd": crystal["rmsd"],
                        "crystal_zinc_N": crystal["zinc_N"],
                        "crystal_vacancy_error": crystal["vacancy_error"],
                        "crystal_Thr199_N_to_sulfonyl_O": crystal[
                            "Thr199_N_to_sulfonyl_O"
                        ],
                        "eligible_decoys": len(decoys),
                        "crystal_minus_best_decoy": (
                            crystal[method] - min(decoys)
                            if crystal[method] is not None and decoys
                            else None
                        ),
                        "winner_model": best["model"] if best else None,
                        "winner_rmsd": best["rmsd"] if best else None,
                    }


def selected_cells(rows):
    key = lambda r: (r["panel"], r["receptor"], r["ligand"], r["seed"])
    for (panel, rec, lig, seed), group in itertools.groupby(sorted(rows, key=key), key):
        group = list(group)
        for method, eligible, best in choices(group):
            yield dict(
                panel=panel,
                receptor=rec,
                ligand=lig,
                seed=seed,
                method=method,
                eligible=len(eligible),
                n_modes=len(group),
                oracle_rmsd=min(r["rmsd"] for r in group),
                winner_model=best["model"] if best else None,
                selected_score=best[method] if best else None,
                **{
                    field: best[field] if best else None
                    for field in (
                        "rmsd",
                        "head_rmsd",
                        "tail_rmsd",
                        "zinc_N",
                        "vacancy_error",
                        "Thr199_N_to_sulfonyl_O",
                    )
                },
            )


def paired_counts(cells, metadata):
    """Count comparable, ≥3-fold pairs inside each receptor/seed only. Abstention stays missing."""
    key = lambda r: (r["panel"], r["receptor"], r["seed"], r["method"])
    for (panel, rec, seed, method), group in itertools.groupby(
        sorted(cells, key=key), key
    ):
        if panel == "extra_cognate":
            continue
        for a, b in itertools.combinations(group, 2):
            a, b = sorted((a, b), key=lambda r: metadata[r["ligand"]]["Ki_nM"])
            ratio = metadata[b["ligand"]]["Ki_nM"] / metadata[a["ligand"]]["Ki_nM"]
            if ratio < 3:
                continue
            x, y = a["selected_score"], b["selected_score"]
            yield {
                "panel": panel,
                "receptor": rec,
                "seed": seed,
                "method": method,
                "stronger": a["ligand"],
                "weaker": b["ligand"],
                "ki_ratio": ratio,
                "correct": (0.5 if x == y else int(x < y))
                if x is not None and y is not None
                else None,
            }


def write_csv(path, rows):
    rows = list(rows)
    if not rows:
        raise ValueError(f"No rows: {path}")
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(
            {k: round(v, 8) if isinstance(v, float) else v for k, v in row.items()}
            for row in rows
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", required=True, type=Path)
    parser.add_argument("--late-pools", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--symmetry-poses", type=Path, help="directory from project_symmetries.py"
    )
    parser.add_argument(
        "--reflected-poses",
        type=Path,
        help="fit of the opposite ammonium state with --reflected",
    )
    parser.add_argument("--binary", required=True, type=Path)
    args = parser.parse_args()
    work = args.work.resolve()
    prepared = work / "prepared"
    metadata = json.loads((prepared / "metadata.json").read_text())
    pools = [
        ("native_search", work / "runs"),
        ("extra_cognate", work / "cognate"),
        ("native_late", args.late_pools),
    ]
    expected = {
        "native_search": set(itertools.product(("6WQ4", "6WQ8"), metadata, (42, 1729))),
        "extra_cognate": {
            (lig, lig, seed) for lig in ("6WQ5", "6WQ7", "6WQ9") for seed in (42, 1729)
        },
        "native_late": set(
            itertools.product(("6WQ4", "6WQ8"), metadata, (271828, 314159))
        ),
    }
    rows = list(pool_rows(prepared, pools, metadata))
    for panel, required in expected.items():
        actual = {
            (r["receptor"], r["ligand"], r["seed"]) for r in rows if r["panel"] == panel
        }
        if actual != required:
            raise ValueError(
                f"{panel}: missing or unexpected receptor/ligand/seed cells: {actual ^ required}"
            )
    cells = list(selected_cells(rows))
    args.output.mkdir(parents=True, exist_ok=False)
    write_csv(args.output / "candidates.csv", rows)
    write_csv(args.output / "selected.csv", cells)
    write_csv(args.output / "pairs.csv", paired_counts(cells, metadata))
    write_csv(
        args.output / "spiked.csv",
        spiked_pools(
            args.binary.resolve(),
            work,
            metadata,
            args.symmetry_poses,
            args.reflected_poses,
        ),
    )


if __name__ == "__main__":
    main()

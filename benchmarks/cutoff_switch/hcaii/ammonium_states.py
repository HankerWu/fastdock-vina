"""Matched trajectory-budget test of both assigned tertiary-ammonium configurations."""

import argparse
import hashlib
import itertools
import json
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from coordination_probe import (
    error_components,
    ligand_mapping,
    load_anchor,
    parse_atoms,
    read_models,
    write_csv,
)
from rdkit import Chem

LIGANDS = ("6WQ4", "6WQ5", "6WQ7", "6WQ8", "6WQ9")
RECEPTORS = ("6WQ4", "6WQ8")
SEEDS = (314159, 271828)


def mirror_state(prepared, ligand):
    """Reflect only entire molecules with one stereoelement, an unspecified tertiary NH+."""
    graph, serials, _ = ligand_mapping(prepared, ligand)
    elements = Chem.FindPotentialStereo(graph)
    if len(elements) != 1 or elements[0].type != Chem.StereoType.Atom_Tetrahedral:
        raise ValueError("Reflection requires exactly one potential stereoelement")
    index = elements[0].centeredOn
    atom = graph.GetAtomWithIdx(index)
    if (
        atom.GetAtomicNum() != 7
        or atom.GetFormalCharge() != 1
        or atom.GetDegree() != 3
        or atom.IsInRing()
        or atom.GetTotalNumHs() != 1
    ):
        raise ValueError("Only singly protonated acyclic tertiary amines are supported")
    source = (prepared / f"{ligand}_ligand.pdbqt").read_text()
    (smiles,) = [
        line.split(maxsplit=2)[2]
        for line in source.splitlines()
        if line.startswith("REMARK SMILES ")
        and not line.startswith("REMARK SMILES IDX")
    ]
    stated = Chem.MolFromSmiles(smiles)
    raw_stereo = Chem.FindPotentialStereo(stated)
    if (
        Chem.MolToSmiles(stated, isomericSmiles=False)
        != Chem.MolToSmiles(graph, isomericSmiles=False)
        or len(raw_stereo) != 1
        or raw_stereo[0].specified != Chem.StereoSpecified.Unspecified
    ):
        raise ValueError(
            "Source SMILES must match the graph and leave its sole stereoelement unspecified"
        )
    lines = []
    for line in source.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            line = line[:30] + f"{-float(line[30:38]):8.3f}" + line[38:]
        lines.append(line)
    signature = (
        serials[index],
        *(
            serials[n.GetIdx()]
            for n in sorted(atom.GetNeighbors(), key=lambda n: n.GetIdx())
        ),
    )
    return "\n".join(lines) + "\n", signature


def signed_volume(points, signature):
    center, *neighbours = signature
    return float(
        np.linalg.det([np.asarray(points[s]) - points[center] for s in neighbours])
    )


def xyzs(text):
    return {s: a.xyz for s, a in parse_atoms(text.splitlines()).items()}


def prepare(work, out, scaffold_anchors):
    prepared = work / "prepared"
    inputs = out / "inputs"
    inputs.mkdir()
    for lig in LIGANDS:
        mirror, _ = mirror_state(prepared, lig)
        (inputs / f"{lig}_reflected.pdbqt").write_text(mirror)
        for rec in RECEPTORS:
            destination = inputs / f"{rec}_{lig}.anchor"
            if scaffold_anchors is not None:
                files = [
                    scaffold_anchors / f"{rec}_{lig}_{seed}_scaffold.anchor"
                    for seed in SEEDS
                ]
                content = files[0].read_bytes()
                if any(path.read_bytes() != content for path in files[1:]):
                    raise ValueError("Anchor coordinates differ between trial seeds")
                destination.write_bytes(content)
            else:
                donor, coords = load_anchor(prepared, rec, lig)
                destination.write_text(
                    f"{donor} " + " ".join(f"{v:.8f}" for v in coords) + "\n"
                )


def input_path(prepared, out, lig, state):
    return (
        (out / "inputs" / f"{lig}_reflected.pdbqt")
        if state == "reflected"
        else (prepared / f"{lig}_ligand.pdbqt")
    )


def dock(item, binary, work, out):
    rec, lig, seed, slot = item
    prepared = work / "prepared"
    state = "reflected" if slot == "second_reflected" else "original"
    search_seed = seed if slot == "common" else seed + 1
    prefix = out / f"{rec}_{lig}_{seed}_{slot}"
    cmd = [
        str(binary),
        "--receptor",
        str(prepared / f"{rec}_receptor.pdbqt"),
        "--config",
        str(prepared / f"{rec}_box.txt"),
        "--ligand",
        str(input_path(prepared, out, lig, state)),
        "--ligand_bonds",
        str(work / "crystal_start" / f"{lig}.bonds"),
        "--seed",
        str(search_seed),
        "--cpu",
        "1",
        "--exhaustiveness",
        "2",
        "--max_evals",
        "50000",
        "--num_modes",
        "9",
        "--energy_range",
        "100",
        "--pose_selection",
        "late",
        "--search_anchor",
        str(out / "inputs" / f"{rec}_{lig}.anchor"),
        "--search_counters",
        "--out",
        str(prefix.with_suffix(".pdbqt")),
    ]
    with prefix.with_suffix(".log").open("w") as log:
        subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)


def fixed_terms(binary, prepared, bonds, receptor, ligand, pose):
    cmd = [
        str(binary),
        "--cpu",
        "1",
        "--seed",
        "42",
        "--receptor",
        str(prepared / f"{receptor}_receptor.pdbqt"),
        "--config",
        str(prepared / f"{receptor}_box.txt"),
        "--ligand",
        str(ligand),
        "--ligand_bonds",
        str(bonds),
        "--pose_coordinates",
        str(pose),
        "--score_only",
        "--net_polar",
        "native",
        "--verbosity",
        "0",
    ]
    output = subprocess.check_output(cmd, text=True)
    inter, intra, _ = map(
        float, re.search(r"^Net polar native: (.+)$", output, re.MULTILINE)[1].split()
    )
    factor = float(
        re.search(r"^Net polar normalization: (.+)$", output, re.MULTILINE)[1]
    )
    return inter, intra, factor


def evaluate(item, binary, work, out, metadata):
    rec, lig, seed, slot = item
    prepared = work / "prepared"
    state = "reflected" if slot == "second_reflected" else "original"
    path = out / f"{rec}_{lig}_{seed}_{slot}.pdbqt"
    ligand = input_path(prepared, out, lig, state)
    _, signature = mirror_state(prepared, lig)
    _, serials, _ = ligand_mapping(prepared, lig)
    reference = np.asarray(metadata[lig]["references"][rec]["coords"])
    true_volume = signed_volume(dict(zip(serials, reference, strict=True)), signature)
    input_volume = signed_volume(xyzs(ligand.read_text()), signature)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    log = path.with_suffix(".log").read_text()
    evaluations = sum(map(int, re.findall(r"\b(?:hunt_evals|second_evals)=(\d+)", log)))
    rows = []
    with TemporaryDirectory() as temp:
        pose = Path(temp) / "pose.pdbqt"
        for model, lines in read_models(path.read_text().splitlines()).items():
            text = "\n".join(lines) + "\n"
            pose.write_text(text)
            atoms = parse_atoms(lines)
            volume = signed_volume({s: a.xyz for s, a in atoms.items()}, signature)
            if input_volume * volume <= 0:
                raise ValueError(
                    f"Configuration changed during docking: {item}, {model}"
                )
            inter, intra, factor = fixed_terms(
                binary,
                prepared,
                work / "crystal_start" / f"{lig}.bonds",
                rec,
                ligand,
                pose,
            )
            rows.append(
                dict(
                    receptor=rec,
                    ligand=lig,
                    seed=seed,
                    slot=slot,
                    state=state,
                    model=model,
                    input_volume=input_volume,
                    bound_volume=true_volume,
                    pose_volume=volume,
                    same_sign_as_bound=int(volume * true_volume > 0),
                    **error_components(prepared, lig, atoms, reference),
                    inter=inter,
                    intra=intra,
                    total=inter + intra,
                    flex_factor=factor,
                    selected_score=inter * factor,
                    grid_evals=evaluations,
                    pool_sha256=digest,
                )
            )
    if not rows:
        raise ValueError(f"Empty docking output: {path}")
    return rows


def evaluate_all(rows, metadata):
    cells = []
    for rec, lig, seed in itertools.product(RECEPTORS, LIGANDS, SEEDS):
        for arm, second in [
            ("original_only", "second_original"),
            ("both_states", "second_reflected"),
        ]:
            pool = [
                r
                for r in rows
                if (r["receptor"], r["ligand"], r["seed"]) == (rec, lig, seed)
                and r["slot"] in ("common", second)
            ]
            selected = min(
                pool, key=lambda r: (r["total"], r["slot"] != "common", r["model"])
            )
            cells.append(
                {
                    "receptor": rec,
                    "ligand": lig,
                    "seed": seed,
                    "arm": arm,
                    "score": selected["selected_score"],
                    "top_rmsd": selected["rmsd"],
                    "best_rmsd": min(r["rmsd"] for r in pool),
                    "best_matching_sign_rmsd": min(
                        (r["rmsd"] for r in pool if r["same_sign_as_bound"]),
                        default=None,
                    ),
                    "head_rmsd": selected["head_rmsd"],
                    "tail_rmsd": selected["tail_rmsd"],
                    "selected_slot": selected["slot"],
                    "selected_model": selected["model"],
                    "chosen_sign_matches_bound": selected["same_sign_as_bound"],
                    "n_modes": len(pool),
                    "grid_evals": sum(
                        next(r["grid_evals"] for r in pool if r["slot"] == s)
                        for s in ("common", second)
                    ),
                }
            )
    pairs = []
    for rec, seed, arm in itertools.product(
        RECEPTORS, SEEDS, ("original_only", "both_states")
    ):
        block = {
            r["ligand"]: r
            for r in cells
            if (r["receptor"], r["seed"], r["arm"]) == (rec, seed, arm)
        }
        for a, b in itertools.combinations(
            sorted(LIGANDS, key=lambda lig: metadata[lig]["Ki_nM"]), 2
        ):
            if metadata[b]["Ki_nM"] / metadata[a]["Ki_nM"] < 3:
                continue
            x, y = block[a]["score"], block[b]["score"]
            pairs.append(
                {
                    "receptor": rec,
                    "seed": seed,
                    "arm": arm,
                    "stronger": a,
                    "weaker": b,
                    "correct": 0.5 if x == y else int(x < y),
                }
            )
    return cells, pairs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("work", "output", "binary"):
        parser.add_argument(f"--{key}", type=Path, required=True)
    parser.add_argument(
        "--scaffold-anchors",
        type=Path,
        help="frozen independent-donor landmark files from the scaffold experiment",
    )
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--evaluate-only", action="store_true")
    args = parser.parse_args()
    work, out, binary = (
        args.work.resolve(),
        args.output.resolve(),
        args.binary.resolve(),
    )
    tasks = list(
        itertools.product(
            RECEPTORS, LIGANDS, SEEDS, ("common", "second_original", "second_reflected")
        )
    )
    if not args.evaluate_only:
        out.mkdir(parents=True, exist_ok=False)
        prepare(work, out, args.scaffold_anchors)
        with ThreadPoolExecutor(args.jobs) as pool:
            jobs = [pool.submit(dock, item, binary, work, out) for item in tasks]
            for job in jobs:
                job.result()
    metadata = json.loads((work / "prepared" / "metadata.json").read_text())
    with ThreadPoolExecutor(args.jobs) as pool:
        jobs = [
            pool.submit(evaluate, item, binary, work, out, metadata) for item in tasks
        ]
        rows = [r for j in jobs for r in j.result()]
    cells, pairs = evaluate_all(rows, metadata)
    for name, data in [("candidates", rows), ("cells", cells), ("pairs", pairs)]:
        write_csv(out / f"{name}.csv", data)


if __name__ == "__main__":
    main()

"""Paired Late docking on fresh seeds with one common native final potential."""

import argparse
import csv
import hashlib
import itertools
import json
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from time import monotonic

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from check_cutoff_1iep import BASIC, crystal_reference


SEEDS = (314159, 271828)


def panel(series, work):
    """Yield input files and a receptor-frame RMSD evaluator for each query."""
    prep = work / "prepared"
    if series == "abl":
        refs, serials = crystal_reference(work / "raw/1IEP.pdb")
        yield (
            "1IEP",
            "1IEP",
            "",
            BASIC / "1iep_receptor.pdbqt",
            BASIC / "1iep_ligand.pdbqt",
            BASIC / "1iep_receptor.box.txt",
            serials,
            reference_rmsd(refs),
        )
        return
    if series == "t4l":
        metadata = json.loads((prep / "manifest.json").read_text())
        receptors = ("1L90", "186L")
    else:
        metadata = json.loads((prep / "metadata.json").read_text())
        receptors = ("6WQ4", "6WQ8")

    for receptor, (ligand, data) in itertools.product(
        receptors, sorted(metadata.items())
    ):
        if series == "t4l":
            measure = data["kd_uM"]
            serials = data["heavy_serials"]
            evaluate = reference_rmsd(
                data["alignments"][receptor]["reference_by_match"]
            )
        else:
            measure = data["Ki_nM"]
            serials = data["qt_serial"]
            molecule = Chem.SDMolSupplier(
                str(prep / f"{ligand}_graph.sdf"), removeHs=True
            )[0]
            for i, xyz in enumerate(data["references"][receptor]["coords"]):
                molecule.GetConformer().SetAtomPosition(i, xyz)
            evaluate = graph_rmsd(molecule)
        yield (
            receptor,
            ligand,
            measure,
            prep / f"{receptor}_receptor.pdbqt",
            prep / f"{ligand}_ligand.pdbqt",
            prep / f"{receptor}_box.txt",
            serials,
            evaluate,
        )


def reference_rmsd(references):
    reference = np.asarray(references)

    def evaluate(position):
        return float(
            np.sqrt(((reference - position) ** 2).sum(axis=2).mean(axis=1)).min()
        )

    return evaluate


def graph_rmsd(reference):
    def evaluate(position):
        probe = Chem.Mol(reference)
        for i, xyz in enumerate(position):
            probe.GetConformer().SetAtomPosition(i, xyz)
        # CalcRMS keeps the receptor frame and handles equivalent terminal groups.
        return rdMolAlign.CalcRMS(probe, reference, maxMatches=1000000)

    return evaluate


def run(binary, output, query, seed, arm):
    rec, lig, measure, receptor, ligand, box, serials, evaluate = query
    stem = output / f"{rec}_{lig}_{seed}_{arm}"
    pose = stem.with_suffix(".pdbqt")
    command = [
        str(binary),
        "--ligand",
        str(ligand),
        "--receptor",
        str(receptor),
        "--config",
        str(box),
        "--cpu",
        "1",
        "--exhaustiveness",
        "4",
        "--seed",
        str(seed),
        "--max_evals",
        "50000",
        "--pose_selection",
        "late",
        "--energy_range",
        "100",
        "--num_modes",
        "9",
        "--search_counters",
        "--out",
        str(pose),
    ]
    if arm == "mixed":
        command.append("--mixed_cutoff_search")
    start = monotonic()
    with stem.with_suffix(".log").open("w") as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    wall = monotonic() - start
    log = stem.with_suffix(".log").read_text()
    energy, rmsd = [], []
    for block in re.findall(
        r"^MODEL\s+\d+\s*\n(.*?)^ENDMDL", pose.read_text(), re.M | re.S
    ):
        atoms = {
            int(line[6:11]): [float(line[i : i + 8]) for i in (30, 38, 46)]
            for line in block.splitlines()
            if line.startswith(("ATOM  ", "HETATM"))
        }
        rmsd.append(evaluate(np.asarray([atoms[a] for a in serials])))
        energy.append(float(re.search(r"^REMARK VINA RESULT:\s+(\S+)", block, re.M)[1]))
    gradients = re.findall(
        r"^Search trajectory: .*?hunt_evals=(\d+) .*?second_evals=(\d+)", log, re.M
    )
    if not energy or len(gradients) != 4:
        raise ValueError(f"Incomplete docking result in {stem}")
    explicit = re.search(
        r"^Search counters: .*?explicit_evals=(\d+) late_prescore=(\d+)", log, re.M
    )
    return dict(
        receptor=rec,
        ligand=lig,
        measurement=measure,
        seed=seed,
        arm=arm,
        score=energy[0],
        top1_rmsd=rmsd[0],
        best_rmsd=min(rmsd),
        pool_rmsds=",".join(f"{v:.6f}" for v in rmsd),
        grid_evals=sum(int(a) + int(b) for a, b in gradients),
        explicit_evals=int(explicit[1]),
        prescore_inputs=int(explicit[2]),
        wall_seconds=round(wall, 4),
        sha256=hashlib.sha256(pose.read_bytes()).hexdigest(),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("series", choices=["t4l", "hcaii", "abl"])
    parser.add_argument("work", type=Path)
    parser.add_argument(
        "output", type=Path, help="New directory for poses, logs and rows.csv"
    )
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    binary = args.binary.resolve()
    jobs = itertools.product(
        list(panel(args.series, args.work)), SEEDS, ("native", "mixed")
    )
    with ThreadPoolExecutor(args.workers) as executor:
        rows = list(executor.map(lambda task: run(binary, args.output, *task), jobs))
    with (args.output / "rows.csv").open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()

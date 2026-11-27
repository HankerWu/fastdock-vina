"""Run the fixed DGJ transfer and verify output coordinates against the deposition."""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import gaa
from study import run_cases

# Historical preparation modules add their own benchmark directories.
sys.path.insert(0, str(HERE.parent))

CASES = {"3GXT": ("A", "NOJ"), "3S5Y": ("A", "DGJ")}
QUERY = "3S5Y"


def save(path, data):
    path.write_text(
        json.dumps(data, indent=2, default=lambda value: value.item()) + "\n"
    )


def audit(work, meta):
    graph = Chem.MolFromSmiles(meta[QUERY]["canonical_smiles"])
    states = []
    for entry, (chain, comp) in CASES.items():
        # prepare() requires a full coordinate-derived chiral match to each CCD.
        reference = Chem.MolFromSmiles(meta[entry]["canonical_smiles"])
        assert graph.GetNumAtoms() == reference.GetNumAtoms() == 11
        match = reference.HasSubstructMatch(graph, useChirality=True)
        (coords,) = gaa.cross.dep_records(
            (work / "raw" / (entry + ".pdb")).read_text().splitlines(), comp
        )
        c4 = np.array(coords["C4"])
        neighbors = np.array([coords[name] for name in ("C3", "C5", "O4")])
        states.append(
            dict(
                entry=entry,
                chain=chain,
                ccd=comp,
                matches_DGJ=match,
                signed_triple_c4_A3=float(np.linalg.det(neighbors - c4)),
            )
        )
    return states


def dock(work, binary, out):
    from summarize import summarize

    meta = json.loads((work / "metadata.json").read_text())
    inputs = audit(work, meta)
    if [row["entry"] for row in inputs if row["matches_DGJ"]] != [QUERY]:
        raise ValueError("Archived receptor and reference chemistry changed")
    jobs, alignments = [], []
    for receiver in CASES:
        rotation, shift, rmsd, count = gaa.cross.alignment(meta[QUERY], meta[receiver])
        refs = np.array(meta[QUERY]["mapped"]) @ rotation + shift
        alignments.append(dict(receiver=receiver, pocket_ca=count, rmsd_A=rmsd))
        for seed in gaa.cross.SEEDS:
            args = [
                "--ligand",
                QUERY + "/protonated.pdbqt",
                "--receptor",
                receiver + "/receptor.pdbqt",
                "--cpu",
                "1",
                "--seed",
                str(seed),
                "--exhaustiveness",
                "4",
                "--max_evals",
                "24000",
                "--num_modes",
                "9",
                "--energy_range",
                "100",
                "--pose_selection",
                "late",
                "--refine_candidates",
                "40",
                "--search_counters",
                "--out",
                "pose.pdbqt",
            ]
            for axis, center in zip("xyz", meta[receiver]["center"]):
                args.extend(["--center_" + axis, str(center), "--size_" + axis, "24"])
            jobs.append(
                (
                    ("GLA", QUERY, receiver, seed, 24000),
                    work,
                    args,
                    None,
                    meta[QUERY]["protonated_serial"],
                    refs,
                )
            )
    run_cases(jobs, binary, out)
    summarize(out)
    save(out / "geometry.json", dict(stereochemistry=inputs, alignments=alignments))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=lambda name: Path(name).resolve())
    sub = parser.add_subparsers(dest="action", required=True)
    sub.add_parser("prepare")
    for action in ("dock", "verify"):
        runner = sub.add_parser(action)
        runner.add_argument("binary", type=lambda name: Path(name).resolve())
        runner.add_argument("output", type=lambda name: Path(name).resolve())
    args = parser.parse_args()
    if args.action == "prepare":
        gaa.prepare(args.work, CASES)
    elif args.action == "dock":
        dock(args.work, args.binary, args.output)
    else:
        from verify import reconstruct

        result = reconstruct(args.output, {"GLA": args.work}, args.binary)
        save(args.output / "independent.json", result)

"""Assess ligand-excluded substructure-guided searches on the shared T4L/BRD4 receptors."""

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

import numpy as np
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
from score_physical_pool import mapped_ligand_bonds, serial_index_pairs

RECEPTORS = {"brd4": ("3MXF", "4OGI", "4OGJ"), "t4l": ("1L90", "186L")}
SEEDS = (314159, 271828)


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def ligand_graph(prepared, ligand):
    text = (prepared / f"{ligand}_ligand.pdbqt").read_text()
    mapped_ligand_bonds(text)
    lines = text.splitlines()
    (smiles,) = [
        line.split()[2]
        for line in lines
        if line.startswith("REMARK SMILES ")
        and not line.startswith("REMARK SMILES IDX ")
    ]
    graph = Chem.MolFromSmiles(smiles)
    serials = {
        i - 1: serial for i, serial in serial_index_pairs(lines, "REMARK SMILES IDX ")
    }
    return graph, serials


def build_anchors(prepared, work, panel):
    from transfer_mapping import choose_reference

    metadata = json.loads((prepared / "manifest.json").read_text())
    mapped = {lig: ligand_graph(prepared, lig) for lig in sorted(metadata)}
    molecules = {lig: value[0] for lig, value in mapped.items()}
    mappings = []
    for query in sorted(metadata):
        match = choose_reference(query, molecules)
        for rec in RECEPTORS[panel]:
            path = work / f"{rec}_{query}.anchor"
            if match is None:
                mappings.append(
                    {
                        "query": query,
                        "donor": "",
                        "receptor": rec,
                        "n": 0,
                        "query_serials": "",
                        "donor_serials": "",
                        "sha256": "",
                    }
                )
                continue
            donor, query_atoms, donor_atoms = match
            if donor == query:
                raise ValueError("Query ligand cannot supply its own anchor")
            query_serials = [mapped[query][1][i] for i in query_atoms]
            donor_serials = [mapped[donor][1][i] for i in donor_atoms]
            donor_reference = metadata[donor]["alignments"][rec]["reference_by_match"][
                0
            ]
            points = dict(
                zip(metadata[donor]["heavy_serials"], donor_reference, strict=True)
            )
            positions = [points[s] for s in donor_serials]
            path.write_text(
                "".join(
                    f"{s} " + " ".join(f"{x:.8f}" for x in xyz) + "\n"
                    for s, xyz in zip(query_serials, positions, strict=True)
                )
            )
            mappings.append(
                {
                    "query": query,
                    "donor": donor,
                    "receptor": rec,
                    "n": len(query_atoms),
                    "query_serials": " ".join(map(str, query_serials)),
                    "donor_serials": " ".join(map(str, donor_serials)),
                    "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                }
            )
    write_csv(work / "mappings.csv", mappings)
    return {(r["receptor"], r["query"]): r["n"] for r in mappings}


def dock(binary, prepared, work, cell, count):
    rec, ligand, seed, arm = cell
    stem = work / f"{rec}_{ligand}_{seed}_{arm}"
    cmd = [
        str(binary.resolve()),
        "--receptor",
        str(prepared / f"{rec}_receptor.pdbqt"),
        "--ligand",
        str(prepared / f"{ligand}_ligand.pdbqt"),
        "--config",
        str(prepared / f"{rec}_box.txt"),
        "--seed",
        str(seed),
        "--cpu",
        "1",
        "--exhaustiveness",
        "4",
        "--max_evals",
        "50000",
        "--num_modes",
        "9",
        "--energy_range",
        "100",
        "--pose_selection",
        "late",
        "--search_counters",
        "--out",
        str(stem.with_suffix(".pdbqt")),
    ]
    if arm == "transfer" and count[(rec, ligand)]:
        cmd += ["--search_anchor", str(work / f"{rec}_{ligand}.anchor")]
    with stem.with_suffix(".log").open("w") as out:
        subprocess.run(cmd, stdout=out, stderr=subprocess.STDOUT, check=True)


def mapping_audit(prepared, work):
    """After docking only, compare chosen donor anchors with the withheld query crystal."""
    metadata = json.loads((prepared / "manifest.json").read_text())
    with (work / "mappings.csv").open(newline="") as stream:
        mappings = list(csv.DictReader(stream))
    diagnostics = []
    for row in mappings:
        if not int(row["n"]):
            continue
        query, donor, rec = row["query"], row["donor"], row["receptor"]
        rows = np.loadtxt(work / f"{rec}_{query}.anchor", ndmin=2)
        query_serials = metadata[query]["heavy_serials"]
        indices = [query_serials.index(int(serial)) for serial in rows[:, 0]]
        crystal = np.array(metadata[query]["alignments"][rec]["reference_by_match"])[
            :, indices, :
        ]
        donor_points = rows[:, 1:]
        rmsds = np.sqrt(np.mean(np.sum((crystal - donor_points) ** 2, axis=2), axis=1))
        match = rmsds.argmin()
        centroid = np.linalg.norm(
            donor_points.mean(axis=0) - crystal[match].mean(axis=0)
        )
        diagnostics.append(
            {
                "query": query,
                "donor": donor,
                "receptor": rec,
                "n": int(row["n"]),
                "query_coverage": len(indices) / len(query_serials),
                "donor_coverage": len(indices) / len(metadata[donor]["heavy_serials"]),
                "anchor_rmsd": float(rmsds[match]),
                "centroid_distance": float(centroid),
            }
        )
    write_csv(work / "mapping_audit.csv", diagnostics)


def evaluate(prepared, work, cells):
    metadata = json.loads((prepared / "manifest.json").read_text())
    modes, selected = [], []
    for rec, ligand, seed, arm in cells:
        stem = work / f"{rec}_{ligand}_{seed}_{arm}"
        ref = np.array(metadata[ligand]["alignments"][rec]["reference_by_match"])
        pool = []
        path = stem.with_suffix(".pdbqt")
        for text in re.findall(
            r"^MODEL .*?^ENDMDL$", path.read_text(), re.MULTILINE | re.DOTALL
        ):
            coords = {
                int(line[6:11]): [float(line[i : i + 8]) for i in (30, 38, 46)]
                for line in text.splitlines()
                if line.startswith(("ATOM  ", "HETATM"))
            }
            xyz = np.array([coords[s] for s in metadata[ligand]["heavy_serials"]])
            rmsd = float(np.sqrt(((ref - xyz) ** 2).sum(axis=2).mean(axis=1)).min())
            score = float(
                re.search(r"^REMARK VINA RESULT:\s+(\S+)", text, re.MULTILINE)[1]
            )
            pool.append(
                {
                    "receptor": rec,
                    "ligand": ligand,
                    "seed": seed,
                    "arm": arm,
                    "model": len(pool) + 1,
                    "score": score,
                    "rmsd": round(rmsd, 6),
                }
            )
        log = stem.with_suffix(".log").read_text()
        counts = re.findall(
            r"^Search trajectory: .*?hunt_evals=(\d+) .*?second_evals=(\d+)",
            log,
            re.MULTILINE,
        )
        if len(counts) != 4 or not pool:
            raise ValueError(f"Incomplete docking {stem}")
        modes.extend(pool)
        selected.append(
            dict(
                **pool[0],
                best_rmsd=min(r["rmsd"] for r in pool),
                n_modes=len(pool),
                grid_evals=sum(int(a) + int(b) for a, b in counts),
                kd_uM=metadata[ligand]["kd_uM"],
                sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            )
        )
    pairs = []
    for key in sorted({(r["receptor"], r["seed"], r["arm"]) for r in selected}):
        group = [r for r in selected if (r["receptor"], r["seed"], r["arm"]) == key]
        for a, b in itertools.combinations(sorted(group, key=lambda x: x["kd_uM"]), 2):
            if b["kd_uM"] < 3 * a["kd_uM"]:
                continue
            gap = b["score"] - a["score"]
            pairs.append(
                {
                    "receptor": key[0],
                    "seed": key[1],
                    "arm": key[2],
                    "strong": a["ligand"],
                    "weak": b["ligand"],
                    "kd_ratio": b["kd_uM"] / a["kd_uM"],
                    "score_gap": gap,
                    "credit": float(gap > 0) + 0.5 * (gap == 0),
                }
            )
    for label, rows in (("modes", modes), ("cells", selected), ("pairs", pairs)):
        write_csv(work / (label + ".csv"), rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage", choices=("prepare", "native", "transfer", "evaluate"))
    p.add_argument("panel", choices=RECEPTORS)
    p.add_argument("prepared", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--binary", type=Path)
    p.add_argument("--workers", type=int, default=2)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.stage == "prepare":
        build_anchors(args.prepared, args.output, args.panel)
        return
    queries = sorted(json.loads((args.prepared / "manifest.json").read_text()))
    arms = ("native", "transfer") if args.stage == "evaluate" else (args.stage,)
    cells = list(itertools.product(RECEPTORS[args.panel], queries, SEEDS, arms))
    if args.stage == "evaluate":
        evaluate(args.prepared, args.output, cells)
        mapping_audit(args.prepared, args.output)
        return
    if not args.binary:
        p.error("native and transfer stages require --binary")
    count = {}
    if args.stage == "transfer":
        count = build_anchors(args.prepared, args.output, args.panel)
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        list(
            executor.map(
                lambda cell: dock(args.binary, args.prepared, args.output, cell, count),
                cells,
            )
        )


if __name__ == "__main__":
    main()

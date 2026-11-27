"""Paired cocrystal runs, same-coordinate rescoring and independent chemical graph audit."""

import argparse
import hashlib
import json
import re
import subprocess
import sys
from itertools import product
from pathlib import Path
from time import perf_counter

from rdkit import Chem
from prepare_inputs import CASES, HERE, panel

sys.path.insert(0, str(HERE.parents[1] / "tools"))
from score_physical_pool import mapped_ligand_bonds


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compare_graph(pdb, ligand, meta, args):
    molecule = Chem.MolFromSmiles(meta["canonical_smiles"])
    serial = meta["serial"]
    chemical = {
        tuple(sorted((serial[b.GetBeginAtomIdx()], serial[b.GetEndAtomIdx()])))
        for b in molecule.GetBonds()
    }
    distances = Chem.GetDistanceMatrix(molecule)
    ix = {s: i for i, s in enumerate(serial)}
    rows = []
    for arm, binary in (
        ("parent", args.parent_graph),
        ("triangle", args.triangle_graph),
    ):
        heavy, pairs = set(), set()
        for line in subprocess.check_output([binary, ligand], text=True).splitlines():
            kind, a, b = line.split()
            (heavy if kind == "bond" else pairs).add((int(a), int(b)))

        def labels(edges):
            return " ".join(f"{a}-{b}" for a, b in sorted(edges))

        short = [p for p in pairs if distances[ix[p[0]], ix[p[1]]] <= 3]
        rows.append(
            dict(
                pdb=pdb,
                arm=arm,
                heavy_bonds=len(heavy),
                missing=labels(chemical - heavy),
                extra=labels(heavy - chemical),
                scored_pairs=len(pairs),
                wrongly_scored=labels(short),
            )
        )
    return rows


def main(args):
    work = args.work.resolve()
    out = work / "runs"
    out.mkdir(exist_ok=True)
    graph, cells, poses, fixed = [], [], [], []
    for pdb in CASES:
        path = work / pdb
        ligand, rec = path / "ligand.pdbqt", path / "receptor.pdbqt"
        meta = json.loads((path / "preparation.json").read_text())
        bonds = path / "bonds.txt"
        bonds.write_text(
            "".join(f"{a} {b}\n" for a, b in mapped_ligand_bonds(ligand.read_text()))
        )
        graph.extend(compare_graph(pdb, ligand, meta, args))
        for seed, cap in product([817139, 191009], [3000, 12000]):
            common = [
                "--ligand",
                ligand,
                "--receptor",
                rec,
                "--cpu",
                "1",
                "--exhaustiveness",
                "1",
                "--max_evals",
                str(cap),
                "--num_modes",
                "9",
                "--energy_range",
                "100",
                "--seed",
                str(seed),
            ]
            for axis, xyz in zip("xyz", meta["center"]):
                common.extend(["--center_" + axis, str(xyz), "--size_" + axis, "24"])
            for arm, binary in (
                ("upstream", args.upstream),
                ("parent", args.parent),
                ("triangle", args.triangle),
                ("explicit", args.parent),
            ):
                key = f"{pdb}_{seed}_{cap}_{arm}"
                output = out / (key + ".pdbqt")
                cmd = [binary, *common, "--out", output]
                if arm != "upstream":
                    cmd += ["--search_counters"]
                if arm == "explicit":
                    cmd += ["--ligand_bonds", bonds]
                t = perf_counter()
                log = subprocess.check_output(cmd, text=True, stderr=subprocess.STDOUT)
                seconds = perf_counter() - t
                (out / (key + ".log")).write_text(log)
                row = dict(
                    key=key,
                    pdb=pdb,
                    seed=seed,
                    cap=cap,
                    arm=arm,
                    seconds=seconds,
                    sha256=sha(output),
                    **panel.paired.counts(log),
                )
                if pdb == "1XMU":
                    panel.paired.RADII.update(Mg=1.73, Zn=1.39)
                poses.extend(panel.paired.inspect(out, meta, meta, row))
                cells.append(row)
            pool = out / (f"{pdb}_{seed}_{cap}_parent.pdbqt")
            for rank, block in enumerate(pool.read_text().split("ENDMDL")[:-1], 1):
                pose = work / "pose.pdbqt"
                pose.write_text(block + "ENDMDL\n")
                for arm, binary in (
                    ("parent", args.parent),
                    ("triangle", args.triangle),
                    ("explicit", args.parent),
                ):
                    cmd = [
                        binary,
                        "--cpu",
                        "1",
                        "--seed",
                        str(seed),
                        "--score_only",
                        "--receptor",
                        rec,
                        "--ligand",
                        ligand,
                        "--pose_coordinates",
                        pose,
                    ]
                    if arm == "explicit":
                        cmd += ["--ligand_bonds", bonds]
                    log = subprocess.check_output(cmd, text=True)
                    vals = list(
                        map(
                            float,
                            re.findall(
                                r":\s+(-?\d+\.\d+)\s+\(kcal/mol\)",
                                log.split("Estimated Free Energy of Binding")[-1],
                            ),
                        )
                    )
                    if len(vals) != 10:
                        raise ValueError(log)
                    fixed.append(
                        dict(
                            pdb=pdb,
                            seed=seed,
                            cap=cap,
                            rank=rank,
                            arm=arm,
                            self_score=vals[0],
                            inter=vals[1],
                            intra=vals[4],
                            zero_reference_numerator=vals[1] + vals[4],
                        )
                    )
    for pdb in ("6WTO", "4LMN"):
        lig = work / "raw" / (pdb + ".pdbqt")
        text = lig.read_text()
        smiles = next(
            l.split(maxsplit=2)[2]
            for l in text.splitlines()
            if l.startswith("REMARK SMILES ") and " IDX " not in l
        )
        idx = [
            int(x)
            for l in text.splitlines()
            if l.startswith("REMARK SMILES IDX ")
            for x in l.split()[3:]
        ]
        mapping = dict(zip(idx[::2], idx[1::2]))
        graph.extend(
            compare_graph(
                pdb,
                lig,
                dict(
                    canonical_smiles=smiles,
                    serial=[mapping[i] for i in sorted(mapping)],
                ),
                args,
            )
        )
    for filename, rows in (
        ("graph", graph),
        ("runs", cells),
        ("modes", poses),
        ("fixed", fixed),
    ):
        panel.paired.write(work / (filename + ".csv"), rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for arg in (
        "work",
        "upstream",
        "parent",
        "triangle",
        "parent_graph",
        "triangle_graph",
    ):
        p.add_argument(arg, type=lambda s: Path(s).resolve())
    main(p.parse_args())

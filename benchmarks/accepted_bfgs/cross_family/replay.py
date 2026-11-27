"""Compare two explicit default local refinements on the historical eight-target panel."""

import argparse
import csv
import hashlib
import json
import os
import re
import resource
import subprocess
import tarfile
import time
from pathlib import Path

import numpy as np
from rdkit import Chem

TARGETS = ("5S8I", "7BLA", "7F5D", "7JGW", "7LT0", "7LZQ", "7MRH", "7MS7")
STARTS = ("projected", "63419", "63421")
INPUTS = {
    "pose/inputs.tar.xz": "fc58a62358a1917990396b96705784374aa43f6734194318aa348bb12afcf92b",
    "basins/starts_and_relaxations.tar.xz": "fafa7d335a89b75d64db78c164398974ecb67ec07f29122263fb5b60eeacb642",
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coords(path):
    atoms = {}
    for line in path.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            serial = int(line[6:11])
            if serial in atoms:
                raise ValueError(f"Duplicate serial {serial} in {path}")
            atoms[serial] = [float(line[k : k + 8]) for k in (30, 38, 46)]
    return atoms


def chemical_graph(ligand, output):
    """Decode the heavy-atom graph and all retained hydrogen bonds from Meeko."""
    lines = ligand.read_text().splitlines()
    smiles = next(
        line[14:]
        for line in lines
        if line.startswith("REMARK SMILES ")
        and not line.startswith("REMARK SMILES IDX")
    )

    def pairs(prefix):
        values = [
            int(x)
            for line in lines
            if line.startswith(prefix)
            for x in line[len(prefix) :].split()
        ]
        return list(zip(values[::2], values[1::2], strict=True))

    indexed = dict(pairs("REMARK SMILES IDX "))
    mol = Chem.MolFromSmiles(smiles)
    if set(indexed) != set(range(1, mol.GetNumAtoms() + 1)):
        raise ValueError(f"Missing SMILES-to-PDBQT atom mapping in {ligand}")
    edges = [
        (indexed[b.GetBeginAtomIdx() + 1], indexed[b.GetEndAtomIdx() + 1])
        for b in mol.GetBonds()
    ]
    edges.extend(
        (indexed[parent], hydrogen) for parent, hydrogen in pairs("REMARK H PARENT ")
    )
    bonds = sorted(tuple(sorted(e)) for e in edges)
    if (
        len(bonds) != len(set(bonds))
        or {v for e in edges for v in e} != coords(ligand).keys()
    ):
        raise ValueError(f"Incomplete graph in {ligand}")
    output.write_text("".join(f"{a} {b}\n" for a, b in bonds))


def rmsd(path, metadata):
    atoms = coords(path)
    xyz = np.asarray([atoms[i] for i in metadata["serial"]])
    return float(
        np.min(
            np.sqrt(
                np.mean(
                    np.sum((np.asarray(metadata["mapped"]) - xyz) ** 2, axis=2), axis=1
                )
            )
        )
    )


def run(exe, args, output, work):
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    start = time.perf_counter()
    with output.open("w") as log:
        subprocess.run(
            [exe, *args], cwd=work, stdout=log, stderr=subprocess.STDOUT, check=True
        )
    elapsed = time.perf_counter() - start
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu = (after.ru_utime - before.ru_utime) + (after.ru_stime - before.ru_stime)
    return elapsed, cpu


def measurements(exe, common, prepared, pose, prefix, work):
    score_log = work / f"{prefix}.score.log"
    run(
        exe,
        [
            *common,
            "--ligand",
            prepared,
            "--pose_coordinates",
            pose,
            "--score_only",
            "--unbound_energy",
            "0",
        ],
        score_log,
        work,
    )
    report = score_log.read_text()
    if not re.search(r"Unbound System.s Energy\s*:\s*0\.000", report):
        raise ValueError(f"Missing explicit zero reference: {score_log}")
    score = float(
        re.search(r"Estimated Free Energy of Binding\s*:\s*([-+.\d]+)", report)[1]
    )
    static = Path(f"{prefix}.one.pdbqt")
    obj_log = work / f"{prefix}.one.log"
    run(
        exe,
        [
            *common,
            "--ligand",
            pose,
            "--local_only",
            "--refinement_optimizer",
            "armijo",
            "--max_evals",
            "1",
            "--out",
            static,
        ],
        obj_log,
        work,
    )
    match = re.search(r"gradient_evals=1 objective_kcal_mol=(\S+)", obj_log.read_text())
    if not match or coords(work / pose) != coords(work / static):
        raise ValueError(f"Objective probe moved or failed: {obj_log}")
    return score, float(match[1])


def write_rows(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def replay(repo, old, new, parser, work):
    work.mkdir(parents=True, exist_ok=False)
    raw = work / "raw"
    raw.mkdir()
    for (name, sha), dest in zip(INPUTS.items(), ("fixtures", "starts"), strict=True):
        archive = repo / "benchmarks/geometry_comparison" / name
        if digest(archive) != sha:
            raise ValueError(f"Expected frozen {archive}")
        with tarfile.open(archive) as tar:
            tar.extractall(work / dest, filter="data")
    if {p.name for p in (work / "starts").iterdir()} != set(TARGETS):
        raise ValueError("Historical panel mismatch")
    protocol = {
        "inputs_sha256": INPUTS,
        "binary_sha256": {
            k: digest(p)
            for k, p in (("parent", old), ("fixed", new), ("topology_parser", parser))
        },
        "pdb": TARGETS,
        "starts": STARTS,
        "cpu_affinity": sorted(os.sched_getaffinity(0)),
    }
    (work / "inputs.json").write_text(json.dumps(protocol, indent=2) + "\n")
    rows = []
    for idx, pdb in enumerate(TARGETS):
        fixture = Path("fixtures") / pdb
        prepared = fixture / "ligand.pdbqt"
        metadata = json.loads((work / fixture / "metadata.json").read_text())
        bonds = Path("raw") / f"{pdb}.bonds"
        chemical_graph(work / prepared, work / bonds)
        signature = subprocess.check_output([parser, prepared, bonds], cwd=work)
        (raw / f"{pdb}.topology.txt").write_bytes(signature)
        common = [
            "--receptor",
            fixture / "receptor.pdbqt",
            "--cpu",
            "1",
            "--seed",
            "63419",
            "--ligand_bonds",
            bonds,
        ]
        for axis, center in zip("xyz", metadata["center"], strict=True):
            common.extend([f"--center_{axis}", str(center), f"--size_{axis}", "24"])
        for pos, label in enumerate(STARTS):
            start = Path("starts") / pdb / f"{label}.pdbqt"
            if subprocess.check_output([parser, start, bonds], cwd=work) != signature:
                raise ValueError(f"Inconsistent prepared topology in {start}")
            arms = [("parent", old), ("fixed", new)]
            if (idx * len(STARTS) + pos) % 2:
                arms.reverse()
            for arm, exe in arms:
                prefix = Path("raw") / f"{pdb}-{label}-{arm}"
                end = prefix.with_suffix(".pdbqt")
                wall, cpu = run(
                    exe,
                    [*common, "--ligand", start, "--local_only", "--out", end],
                    work / prefix.with_suffix(".local.log"),
                    work,
                )
                if subprocess.check_output([parser, end, bonds], cwd=work) != signature:
                    raise ValueError(f"Inconsistent output topology in {end}")
                e0, obj0 = measurements(
                    exe, common, prepared, start, f"{prefix}-start", work
                )
                e1, obj1 = measurements(
                    exe, common, prepared, end, f"{prefix}-end", work
                )
                rows.append(
                    {
                        "pdb": pdb,
                        "target": metadata["target"],
                        "start": label,
                        "arm": arm,
                        "start_sha256": digest(work / start),
                        "end_sha256": digest(work / end),
                        "start_rmsd_A": rmsd(work / start, metadata),
                        "end_rmsd_A": rmsd(work / end, metadata),
                        "start_zero_score": e0,
                        "end_zero_score": e1,
                        "start_objective": obj0,
                        "end_objective": obj1,
                        "wall_s": wall,
                        "cpu_s": cpu,
                    }
                )
    write_rows(work / "rows.csv", rows)
    with tarfile.open(work / "raw.tar.xz", "w:xz") as tar:
        tar.add(raw, arcname="raw")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    for name in ("repo", "old", "new", "parser", "work"):
        ap.add_argument(name, type=lambda value: Path(value).resolve())
    args = ap.parse_args()
    replay(args.repo, args.old, args.new, args.parser, args.work)

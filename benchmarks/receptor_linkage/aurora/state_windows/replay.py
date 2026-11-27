"""Evaluate identical fixed ligand centers in two rigid experimental receptor structures."""
import argparse
import csv
import hashlib
import io
import json
import resource
import subprocess
import sys
import tarfile
import time
from pathlib import Path

import numpy as np
from rdkit import Chem

AURORA = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(AURORA))
import geometry as geom
from pose_probe import mapped_rmsd, store_pose

KINASE = AURORA.parents[1] / "torsion_kinase_affinity/evidence.tar.xz"
SOURCES = {
    AURORA / "probe_raw.tar.xz": "7066127b96ade03756571adfa9ef711294b2e9ac80d4f0ed735c460120fd9b9b",
    AURORA / "docking.tar.xz": "783f8d17ae31ed32d9fbe167e2f50e455ea54255f033c4506651577fb0a8af7d",
    AURORA / "references.tar.xz": "a653a77c1fff3201668cf87ad3dfe859c8d93307ec8e397ef2a383ad9673d9fc",
    AURORA / "solver_replay/outputs.tar.xz": "1924587e8f10a3551eabc941d9f169ac55c859b227ebba263fb74729c42c9353",
    KINASE: "72c45e52eda2b79d5241b00b4efe79b2e2c19dc100f01a3ad52356d9c50cf56e",
}
STATES = ("2J4Z", "2X6E")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_inputs(work):
    for archive, expect in SOURCES.items():
        if digest(archive) != expect:
            raise ValueError(f"Input archive differs: {archive}")
        if archive != KINASE:
            with tarfile.open(archive) as tf:
                tf.extractall(work, filter="data")
    with tarfile.open(KINASE) as tf:
        whole = {
            r: geom.ca_atoms(tf.extractfile(f"raw/{r}.pdb").read().decode().splitlines(), "A")
            for r in STATES
        }
    return {"pocket": geom.receiver_frames(KINASE), "whole": whole}


def run(args):
    p = subprocess.run([str(x) for x in args], capture_output=True, text=True)
    if p.returncode:
        raise ValueError(f"{args}: {p.stdout} {p.stderr}")
    return p.stdout


def first_model(path):
    lines = path.read_text().split("ENDMDL")[0].splitlines()
    return [s for s in lines if s.strip() and not s.startswith("MODEL")]


def localize(binary, work, origin, drug, bonds, out):
    prefix = "sdf" if drug == "at9283" else "native"
    start = work / "raw" / f"{origin}_{drug}_{prefix}_legal_placed.pdbqt"
    box = json.loads((work / "prepared" / f"{origin}.json").read_text())["center"]
    box_args = [arg for axis, v in zip("xyz", box, strict=True)
                for arg in (f"--center_{axis}", v, f"--size_{axis}", 26)]
    message = run([binary, "--receptor", work / "prepared" / f"{origin}.pdbqt", "--ligand", start,
                   "--ligand_bonds", bonds, "--cpu", 1, "--seed", 60913, *box_args,
                   "--local_only", "--out", out])
    out.with_suffix(".log").write_text(message)


def make_pool(poses, rotation, shift, pool):
    parts = []
    for i, path in enumerate(poses, start=1):
        lines = first_model(path)
        pts = np.array([geom.xyz(s) for s in lines if s.startswith(("ATOM  ", "HETATM"))])
        single = pool.with_name(f"{pool.stem}_{i}.pdbqt")
        store_pose(single, lines, pts @ rotation + shift)
        parts.append(f"MODEL {i}\n{single.read_text()}ENDMDL\n")
    pool.write_text("".join(parts))


def experiment(binary, work, out, counts):
    work.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    wall, cpu = time.perf_counter(), time.process_time()
    frames = load_inputs(work)
    prep = work / "prepared"
    outputs = work / "coordinates"
    outputs.mkdir(exist_ok=True)
    rows, alignments, centers = [], {}, []
    for drug, entries in geom.REFERENCES.items():
        reference, component, _ = entries[0]
        crystal = (work / f"{reference}.pdb").read_text().splitlines()
        ligand = prep / f"{drug}.pdbqt"
        mol = Chem.RemoveHs(Chem.MolFromMolFile(str(ligand.with_suffix(".sdf")), removeHs=False))
        targets, _ = geom.reference_mappings(mol, geom.ligand_atoms(crystal, component, "A"),
                                            geom.ccd_data(work, component), crystal)
        serials, _ = geom.serial_map(ligand, mol)
        graph = "chemical" if drug == "at9283" else "native"
        bonds = work / "raw" / f"{drug}_{graph}.bonds"
        for origin in STATES:
            key = f"{origin}_{drug}"
            witness = outputs / f"{key}_witness.pdbqt"
            localize(binary, work, origin, drug, bonds, witness)
            poses = [witness] + [work / "pools" / f"pro_{key}_{s}.pdbqt" for s in geom.SEEDS]
            labels = ["witness", *map(str, geom.SEEDS)]
            r, t, _ = geom.transform(frames["pocket"][origin], geom.ca_atoms(crystal, "A"))
            for label, path in zip(labels, poses, strict=True):
                centers.append(dict(origin=origin, ligand=drug, pose=label,
                                    reference=reference, source_sha256=digest(path),
                                    rmsd_A=mapped_rmsd(geom.mode_coords(path, serials)[0], targets @ r + t)))
            for state in STATES:
                for frame, ca in frames.items():
                    if origin == state and frame == "whole":
                        continue
                    r, t, fit = geom.transform(ca[state], ca[origin])
                    alignments[f"{frame}_{origin}_{state}"] = fit
                    pool = outputs / f"{key}_{state}_{frame}.pdbqt"
                    make_pool(poses, r, t, pool)
                    for n in counts:
                        stdout = run([binary, "--receptor", prep / f"{state}.pdbqt", "--ligand", ligand,
                                      "--ligand_bonds", bonds, "--cpu", 1, "--score_only",
                                      "--rigid_window_pool", pool, "--rigid_window_samples", n, "--verbosity", 0])
                        records = csv.DictReader(io.StringIO(stdout[stdout.index("model,inter,intra,"):]))
                        for i, (label, row) in enumerate(zip(labels, records, strict=True), start=1):
                            if int(row["model"]) != i:
                                raise ValueError("Unexpected model order")
                            score = (float(row["inter"]) + float(row["intra"])) * float(row["slope"])
                            rows.append(dict(origin=origin, ligand=drug, pose=label, state=state,
                                             registration=frame, nodes_per_scale=n, center=score,
                                             window=score + float(row["window_delta"]), **row))
    geom.write_csv(out, "state_scores.csv", rows)
    geom.write_csv(out, "centers.csv", centers)
    (out / "registration.json").write_text(json.dumps(alignments, indent=2) + "\n")
    sub = resource.getrusage(resource.RUSAGE_CHILDREN)
    (out / "cost.json").write_text(json.dumps(dict(
        binary_sha256=digest(binary), nodes_per_scale=counts,
        quadrature_node_evaluations=len(rows) // len(counts) * 3 * sum(counts),
        local_calls=4, wall_s=time.perf_counter()-wall,
        cpu_s=time.process_time()-cpu+sub.ru_utime+sub.ru_stime,
        child_maxrss_KiB=sub.ru_maxrss,
        input_archives={str(p.relative_to(AURORA.parents[1])):h for p,h in SOURCES.items()},
    ), indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("binary", type=Path, help="binary built from 5015481")
    p.add_argument("work", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("--nodes", nargs="+", type=int, default=[2048, 8192, 32768, 65536])
    a = p.parse_args()
    experiment(a.binary.resolve(), a.work.resolve(), a.out.resolve(), a.nodes)

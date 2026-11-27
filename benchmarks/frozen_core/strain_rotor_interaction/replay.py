"""Replay the frozen first-pose strain × contact-rotor comparison without rerunning search."""
import argparse
import csv
import json
from decimal import Decimal
import lzma
import sys
import tarfile
from pathlib import Path

from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
CORE = HERE.parent
BENCH = CORE.parent
sys.path.insert(0, str(CORE / "affinity"))
import test_entropy as entropy


def load(path):
    opener = lzma.open if path.suffix == ".xz" else open
    with opener(path, "rt", newline="") as stream:
        return list(csv.DictReader(stream))


def save(path, rows):
    opener = lzma.open if path.suffix == ".xz" else open
    with opener(path, "wt", newline="") as stream:
        writer = csv.DictWriter(stream, rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def index(rows, keys):
    result = {tuple(r[k] for k in keys): r for r in rows}
    if len(result) != len(rows):
        raise ValueError(f"Duplicate row for {keys}")
    return result


def observation(panel, lab, pose, scale, unit, native, contact, raw):
    """Remove native torsion scaling before combining the two substitutions."""
    native, unit, contact, raw, scale = map(float, (native, unit, contact, raw, scale))
    if abs(raw * scale - native) > .00151:
        raise ValueError(f"Incompatible normalization in {panel}: {lab}, {pose}")
    return dict(panel=panel, **lab, receptor=pose["receptor"], seed=pose["seed"],
                native=native, strain=unit, rotor=contact,
                joint=contact + (unit - native) / scale, strain_unscaled_no_rotor=raw + (unit - native) / scale,
                raw=raw, scale=scale)


def kinase_bace():
    labels = index(load(BENCH / "geometry_comparison/affinity/panel.csv"), ["target", "assay", "compound"])
    rotor_meta = index(load(CORE / "affinity/native_rotors.csv"), ["compound"])
    keys = ["receptor", "compound", "seed"]
    pooled = index([r for r in load(CORE / "affinity/poses.csv.xz")
                    if r["rank"] == "1" and r["optimizer"] == "armijo"], ["budget", *keys])
    fresh = index(load(CORE / "native/cells.csv"), keys)
    expected = {"native_12k": 140, "armijo_12000": 140, "armijo_60000": 140}
    rows = []
    for file in ("fresh_mmff_modes.csv.xz", "affinity_modes.csv.xz"):
        for r in load(CORE / "strain" / file):
            if r["rank"] != "1":
                continue
            mol = r["compound"]
            key = tuple(r[k] for k in keys)
            n = rotor_meta[mol,]
            m = labels[r["target"], r["assay"], mol]
            lab = dict(target=r["target"], assay=r["assay"], endpoint=m["endpoint"], compound=mol,
                       measured_nM=m["value_nM"], heavy=n["heavy_atoms"], charge=m["formal_charge"])
            if "cap" in r:
                arm, unit, panel = pooled[r["cap"], *key], r["full_unit"], "armijo_" + r["cap"]
                scale = float(r["flexibility_scale"])
            else:
                arm, unit, panel = fresh[key], r["unit_8"], "native_12k"
                scale = 1 / (1 + .05846 * float(n["num_tors"]))
            if float(arm["native"]) != float(r["native"]):
                raise ValueError(f"Native pose mismatch: {panel} {key}")
            rows.append(observation(panel, lab, r, scale, unit, r["native"], arm["contact"], arm["raw"]))
    for panel, count in expected.items():
        if sum(r["panel"] == panel for r in rows) != count:
            raise ValueError(f"Missing frozen first poses for {panel}")
    return rows


def fxa():
    measured = index([r for r in load(BENCH / "affinity_partition/fxa_series.csv") if r["included"] == "1"], ["ligand"])
    keys = ["compound", "receptor", "seed"]
    contacts = index([r for r in load(CORE / "affinity/fxa_poses.csv.xz") if r["rank"] == "1"], keys)
    repair = index([dict(**r, compound="CHEMBL369480") for r in
                    load(CORE / "strain_transfer_audit/redock_cells.csv") if r["arm"] == "random_coords"], keys)
    probe = index(load(BENCH / "terminal_rotor_score/results/distorted_input.csv"), ["input"])
    repaired_tors = float(probe["random_coords.pdbqt",]["old"])
    repaired_scale = 1 / (1 + .05846 * repaired_tors)
    # The repaired coordinate package changes exactly one ligand and keeps the same case keys.
    with tarfile.open(CORE / "strain_transfer_audit/redock_raw.tar.xz") as poses, tarfile.open(
            CORE / "affinity/fxa_frozen.tar.xz") as receivers:
        corrected_topology = entropy.ligand_graph(poses.extractfile("random_coords.pdbqt").read().decode())
        trees = {}
        for rec in ("1FJS", "2P16"):
            lines = receivers.extractfile(f"{rec}/{rec}.pdbqt").read().decode()
            atoms = map(entropy.atom, entropy.atom_lines(lines))
            trees[rec] = cKDTree([a[2] for a in atoms if entropy.heavy(a)])
        rows = []
        for r in load(CORE / "strain/fxa_modes.csv.xz"):
            if r["rank"] != "1":
                continue
            key = tuple(r[k] for k in keys)
            compound, rec, seed = key
            m, arm = measured[compound,], contacts[key]
            lab = dict(target="FXa", assay="CHEMBL661804", endpoint="Ki", compound=compound,
                       measured_nM=m["ki_nm"], heavy=m["heavy"], charge=m["charge"])
            if r["native"] != arm["native"]:
                raise ValueError(f"Different FXa first poses at {key}")
            saved = observation("fxa_original", lab, r, r["flexibility_scale"], r["full_unit"],
                                r["native"], arm["contact"], arm["raw"])
            rows.append(saved)
            if key in repair:
                fixed = repair[key]
                raw = poses.extractfile(f"random_coords_{rec}_{seed}.pdbqt").read().decode().split("ENDMDL")[0]
                contact = entropy.pose(raw, trees[rec], corrected_topology, repaired_tors)
                if float(fixed["native"]) != contact["native"]:
                    raise ValueError("Corrected native pose mismatch")
                # For a native first mode I=U; check the actual block before using the energy ledger.
                intra = next(float(l.split()[-1]) for l in raw.splitlines() if l.startswith("REMARK INTRA:"))
                unbound = next(float(l.split()[-1]) for l in raw.splitlines() if l.startswith("REMARK UNBOUND:"))
                native = float(fixed["native"])
                delta = float(fixed["bound"]) - float(fixed["free"]) - (intra - unbound)
                saved = observation("fxa_repaired", lab, r, repaired_scale, native + repaired_scale * delta,
                                    fixed["native"], contact["contact"], contact["raw"])
            rows.append({**saved, "panel": "fxa_repaired"})
    if len(rows) != 224 or len(repair) != 4:
        raise ValueError("Incomplete FXa comparison")
    return rows


def parp():
    folder = CORE / "transfer/parp_affinity"
    keys = ["compound", "receptor", "seed"]
    mmff = index(load(folder / "mmff.csv"), keys)
    rot = index(load(folder / "native_rotors.csv"), ["compound"])
    with tarfile.open(folder / "inputs.tar.xz") as prep:
        chemistry = {r["compound"]: r for r in json.load(prep.extractfile("chemistry.json"))}
    rows = []
    for r in load(folder / "selected.csv"):
        k = tuple(r[key] for key in keys)
        unit = mmff[k]
        mol = r["compound"]
        if r["pool_sha256"] != unit["pool_sha256"]:
            raise ValueError("Different PARP1 pose pools")
        lab = dict(target="PARP1", assay="Ryan_2021_SPR", endpoint="Kd", compound=mol,
                   measured_nM=r["kd_nM"], heavy=rot[mol,]["heavy_atoms"], charge=chemistry[mol]["charge"])
        rows.append(observation("parp_spr", lab, r, unit["scale"], unit["unit_score"],
                                r["native"], r["contact"], r["raw"]))
    if len(rows) != 16:
        raise ValueError("Incomplete PARP1 comparison")
    return rows


def transfer(path):
    rows = []
    for r in load(path):
        if r["rotor_status"] != "ok":
            continue
        gal = r["target"] == "Gal3"
        lab = dict(target=r["target"], assay="Sorme_2005_FP" if gal else "Davis-S6K1", endpoint="Kd",
                   compound=r["compound"], measured_nM=Decimal(r["measurement"]) * (1000 if gal else 1),
                   heavy=r["heavy"], charge=r["charge"])
        old = observation("transfer_archived", lab, r, r["scale"], r["unit"], r["native"], r["contact"], r["raw"])
        rows.append(old)
        f = 1 / (1 + .05846 * float(r["active_heavy_tors"]))
        raw = float(r["raw"])
        strain_delta = (old["strain"] - old["native"]) / old["scale"]
        rows.append(observation("transfer_active_torsions", lab, r, f, f * (raw + strain_delta), f * raw,
                                r["contact"], raw))
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("out", type=Path)
    p.add_argument("--transfer", type=Path, default=HERE / "transfer.csv")
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    rows = kinase_bace() + fxa() + parp() + transfer(args.transfer)
    index(rows, ["panel", "target", "assay", "receptor", "compound", "seed"])
    save(args.out / "first_poses.csv.xz", rows)


if __name__ == "__main__":
    main()

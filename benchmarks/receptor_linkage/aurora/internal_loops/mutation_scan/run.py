"""Frozen search-order comparison. Only the evaluator receives crystal coordinates."""

import argparse
import importlib.util
import itertools
import json
import os
import re
import resource
import subprocess
import sys
import time
from pathlib import Path

import gemmi
import numpy as np
from rdkit import Chem

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[3]
sys.path.insert(0, str(HERE.parent))
import replay
from pose_probe import get_tree, mapped_rmsd


def load_path(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def source_graph(ligand, mol):
    _, ids, indices, _, _ = get_tree(ligand)
    serial = dict(zip(indices, ids, strict=True))
    return sorted(
        tuple(sorted((serial[a], serial[b])))
        for bond in mol.GetBonds()
        for a, b in [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())]
        if a in serial and b in serial
    )


def config(aurora, brd9, out):
    cells = {}
    for drug in ("at9283", "danusertib"):
        ligand, bonds, ids, refs = replay.parameters(aurora, drug)
        mol = Chem.MolFromMolFile(str(ligand.with_suffix(".sdf")), removeHs=False)
        if sorted(
            tuple(sorted(map(int, line.split())))
            for line in bonds.read_text().splitlines()
        ) != source_graph(ligand, mol):
            raise ValueError(f"Explicit graph differs from SDF for {drug}")
        cells[drug] = dict(
            ligand=ligand,
            bonds=bonds,
            receptor=aurora / replay.PROTEINS["1OL7_released"],
            center=json.loads((aurora / "prepared/2X6E.json").read_text())["center"],
            size=26,
            serials=ids,
            refs=refs,
            reference="2W1G:A" if drug == "at9283" else "2J50:A",
        )
    cross = load_path("crossdocking", BENCH / "paired_crossdocking/run.py")
    torsions = load_path("torsions_paired", BENCH / "torsion_proposals/paired.py")
    meta = json.loads((brd9 / "metadata.json").read_text())
    q, r = meta["6YQW"], meta["6Y7L"]
    ligand = brd9 / "6YQW/ligand.pdbqt"
    ccd = brd9 / "raw/82I.cif"
    if replay.digest(ccd) != q["ccd_sha256"]:
        raise ValueError("BRD9 CCD identity")
    descriptors = (
        gemmi.cif.read_file(str(ccd))
        .sole_block()
        .get_mmcif_category("_pdbx_chem_comp_descriptor.")
    )
    smiles = next(
        s
        for typ, s in zip(descriptors["type"], descriptors["descriptor"])
        if typ == "SMILES_CANONICAL"
    )
    canonical = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
    if canonical != q["canonical_smiles"]:
        raise ValueError("Unexpected query chemistry")
    edges = source_graph(ligand, Chem.AddHs(Chem.MolFromSmiles(canonical)))
    bonds = out / "82I.bonds"
    bonds.write_text("".join(f"{a} {b}\n" for a, b in edges))
    rot, t, _, _ = cross.alignment(*torsions.aligned_receivers(q, r))
    cells["82I"] = dict(
        ligand=ligand,
        bonds=bonds,
        receptor=brd9 / "6Y7L/receptor.pdbqt",
        center=r["center"],
        size=24,
        serials=q["serial"],
        refs=np.asarray(q["mapped"]) @ rot + t,
        reference="6YQW:A",
    )
    return cells


def command(binary, cell, seed, output, official=False, cap=25000):
    cmd = [binary, "--receptor", cell["receptor"], "--ligand", cell["ligand"]]
    for dim, center in zip("xyz", cell["center"], strict=True):
        cmd += [f"--center_{dim}", center, f"--size_{dim}", cell["size"]]
    cmd += [
        "--seed",
        seed,
        "--exhaustiveness",
        2,
        "--cpu",
        1,
        "--num_modes",
        9,
        "--energy_range",
        100,
        "--max_evals",
        cap,
        "--out",
        output,
    ]
    if not official:
        cmd += [
            "--ligand_bonds",
            cell["bonds"],
            "--pose_selection",
            "late",
            "--search_counters",
        ]
    return list(map(str, cmd))


def metrics(cell, pose):
    errors = [
        float(mapped_rmsd(xyz, cell["refs"]))
        for xyz in replay.geom.mode_coords(pose, cell["serials"])
    ]
    scores = [
        float(line.split()[3])
        for line in pose.read_text().splitlines()
        if line.startswith("REMARK VINA RESULT:")
    ]
    return list(zip(errors, scores, strict=True))


def study(args):
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    inputs = config(args.aurora.resolve(), args.brd9.resolve(), out)
    raw = out / "raw"
    raw.mkdir()
    manifests = {
        k: dict(
            **{key: replay.digest(v[key]) for key in ["ligand", "bonds", "receptor"]},
            center=v["center"],
            size=v["size"],
            reference=v["reference"],
        )
        for k, v in inputs.items()
    }
    manifests["binaries"] = {
        k: replay.digest(getattr(args, k)) for k in ("native", "scan")
    }
    (out / "inputs.json").write_text(json.dumps(manifests, indent=2) + "\n")
    cells, modes, draws = [], [], []
    for drug, cell in inputs.items():
        seeds = [61231, 71261, 71263] if drug == "at9283" else [71261, 71263]
        for iteration, seed in enumerate(seeds):
            ordering = ["native", "iid", "scan"]
            for arm in ordering[iteration % 3 :] + ordering[: iteration % 3]:
                binary = args.native if arm == "native" else args.scan
                pose = raw / f"{drug}_{seed}_{arm}.pdbqt"
                log = pose.with_suffix(".log")
                env = {
                    k: v
                    for k, v in os.environ.items()
                    if k not in ("FD_SCAN_IID", "FD_SCAN_TRACE")
                }
                if arm in ("iid", "scan"):
                    env["FD_SCAN_TRACE"] = "1"
                if arm == "iid":
                    env["FD_SCAN_IID"] = "1"
                start = time.monotonic()
                before = resource.getrusage(resource.RUSAGE_CHILDREN)
                with log.open("w") as f:
                    subprocess.run(
                        command(binary.resolve(), cell, seed, pose),
                        env=env,
                        stdout=f,
                        stderr=subprocess.STDOUT,
                        check=True,
                    )
                after = resource.getrusage(resource.RUSAGE_CHILDREN)
                wall = time.monotonic() - start
                cpu = (
                    after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime
                )
                label = dict(drug=drug, seed=seed, arm=arm)
                ranked = metrics(cell, pose)
                for rank, (rmsd, score) in enumerate(ranked, 1):
                    modes.append(dict(**label, rank=rank, rmsd_A=rmsd, score=score))
                errors = [r for r, _ in ranked]
                text = log.read_text()
                evals = sum(
                    map(
                        int, re.findall(r"\b(?:hunt|second|explicit)_evals=(\d+)", text)
                    )
                )
                cells.append(
                    dict(
                        **label,
                        n=len(errors),
                        top1=errors[0],
                        best5=min(errors[:5]),
                        oracle=min(errors),
                        regret=errors[0] - min(errors),
                        evals=evals,
                        wall_s=wall,
                        cpu_s=cpu,
                        pool_sha256=replay.digest(pose),
                    )
                )
                trails = re.findall(r"Mutation draws n=(\d+):([ 0-9]+)", text)
                if arm != "native" and len(trails) != 2:
                    raise ValueError(f"Missing draws {label}")
                hunts = [
                    int(n)
                    for n in re.findall(
                        r"Search trajectory: index=\d+ steps=(\d+)", text
                    )
                ]
                for idx, (n, raw_seq) in enumerate(trails):
                    n = int(n)
                    seq = np.fromstring(raw_seq, sep=" ", dtype=int)
                    if len(seq) != hunts[idx] or seq.min() < 0 or seq.max() >= n:
                        raise ValueError(f"Wrong number/range of draws {label}")
                    blocks = seq[: len(seq) // n * n].reshape(-1, n)
                    missing = sum(n - len(set(block)) for block in blocks)
                    if arm == "scan" and missing:
                        raise ValueError("Permutation invariant failed")
                    for entity in range(n):
                        hits = np.flatnonzero(seq == entity)
                        longest = np.diff(np.r_[-1, hits, len(seq)]).max() - 1
                        draws.append(
                            dict(
                                **label,
                                trajectory=idx + 1,
                                entity=entity,
                                entities=n,
                                steps=len(seq),
                                selected=len(hits),
                                longest_omission=int(longest),
                                complete_blocks=len(blocks),
                                missing_opportunities=int(missing),
                            )
                        )
            group = {
                r["arm"]: r for r in cells if r["drug"] == drug and r["seed"] == seed
            }
            for key in ("evals", "pool_sha256"):
                if group["native"][key] != group["iid"][key]:
                    raise ValueError(f"Instrumented IID diverged {drug} {seed} {key}")
    for name, rows in (
        ("cells.csv", cells),
        ("modes.csv", modes),
        ("draws.csv", draws),
    ):
        replay.geom.write_csv(out, name, rows)
    distributions = []
    for arm, drug in itertools.product(("native", "scan"), ("all", *inputs.keys())):
        selected = [
            r
            for r in cells
            if r["arm"] == arm and r["seed"] != 61231 and drug in (r["drug"], "all")
        ]
        top = np.array([r["top1"] for r in selected])
        q = np.percentile(top, [25, 50, 75, 90, 95])
        distributions.append(
            dict(
                arm=arm,
                drug=drug,
                cells=len(top),
                mean=top.mean(),
                **dict(zip(["q25", "median", "q75", "p90", "p95"], q)),
                **{
                    f"mean_{k}": np.mean([r[k] for r in selected])
                    for k in ("best5", "oracle", "regret")
                },
                severe_ge_8=int((top >= 8).sum()),
                empty=sum(r["n"] == 0 for r in selected),
            )
        )
    replay.geom.write_csv(out, "distributions.csv", distributions)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("aurora", "brd9", "native", "scan", "out"):
        p.add_argument(arg, type=Path)
    study(p.parse_args())

"""Evaluate fixed, published pools only. Read the HUP crystal after docking."""
import argparse
import csv
import json
import lzma
import tarfile
import tempfile
from pathlib import Path

import numpy as np

from dock import ARCHIVE, ARMS, RECEIVERS, SEEDS, TRANSFER, chemistry, pocket, save_json, sha, table, tick


def references(work):
    source = lzma.decompress((TRANSFER / "4EY5.pdb.xz").read_bytes())
    past = json.loads((TRANSFER / "manifest.json").read_text())
    if sha(source) != past["hashes"]["4EY5.pdb"]:
        raise ValueError("Changed crystallographic reference")
    rows = source.decode().splitlines()
    source_pocket = pocket(rows, "HUP")
    mapped = chemistry.ligand(TRANSFER / "HUP.cif", chemistry.dep_records(rows, "HUP"), work / "mapped.pdbqt")
    if (work / "mapped.pdbqt").read_bytes() != (TRANSFER / "HUP.pdbqt").read_bytes():
        raise ValueError("Reference-to-query chemistry mapping disagrees")
    frames, fits = {}, {}
    with tarfile.open(ARCHIVE) as archive:
        for rec, code in RECEIVERS.items():
            lines = archive.extractfile(f"novel/raw/{rec}.pdb").read().decode().splitlines()
            rot, trans, residual, count = chemistry.alignment(source_pocket, pocket(lines, code))
            frames[rec] = np.array(mapped["mapped"]) @ rot + trans
            fits[rec] = dict(rmsd_A=residual, shared_ca=count, mappings=len(frames[rec]),
                             rotation=rot.tolist(), translation_A=trans.tolist())
    return frames, mapped["serial"], fits


def distribution(values):
    if not len(values):
        return {k: None for k in ("mean", "q1", "median", "q3", "p90", "p95")}
    q = np.quantile(values, [.25, .5, .75, .90, .95], method="linear")
    return dict(mean=np.mean(values), **dict(zip(("q1", "median", "q3", "p90", "p95"), q, strict=True)))


def report(result):
    start = tick()
    with tempfile.TemporaryDirectory() as tmp:
        ref, serial, fits = references(Path(tmp))
    cells, modes, all_xyz = [], [], {}
    for run in csv.DictReader((result / "runs.csv").open()):
        rec, arm, seed = run["receiver"], run["arm"], int(run["seed"])
        label = dict(receiver=rec, arm=arm, seed=seed)
        directory = result / "raw" / rec / f"{arm}_{seed}"
        model_errors = []
        all_xyz[rec, arm, seed] = []
        if (directory / "dry.pdbqt").is_file():
            data = (directory / "dry.pdbqt").read_bytes()
            if sha(data) != run["output_sha256"]:
                raise ValueError(f"Pose digest disagrees for {directory}")
            for rank, block in enumerate(data.decode().split("ENDMDL")[:-1], 1):
                lines = block.splitlines()
                score, = [float(l.split()[3]) for l in lines if l.startswith("REMARK VINA RESULT:")]
                sites = {int(l[6:11]): chemistry.coords(l) for l in lines if l.startswith(("ATOM", "HETATM"))}
                xyz = np.array([sites[i] for i in serial])
                d = float(np.sqrt(np.sum((xyz-ref[rec])**2, axis=2).mean(axis=1)).min())
                all_xyz[rec, arm, seed].append(xyz)
                model_errors.append(d)
                modes.append(dict(**label, mode=rank, score=score, rmsd_A=d))
        d = model_errors
        cells.append(dict(**label, modes=len(d), status=run["status"], top1_A=d[0] if d else None,
                          best5_A=min(d[:5]) if d else None, oracle_A=min(d) if d else None,
                          oracle_mode=int(np.argmin(d)+1) if d else None, regret_A=d[0]-min(d) if d else None,
                          serious=d[0]>5 if d else None))
    table(result / "modes.csv", modes)
    table(result / "cells.csv", cells)
    groups = []
    for (rec, arm) in dict.fromkeys((c["receiver"], c["arm"]) for c in cells):
        sample = [c for c in cells if (c["receiver"], c["arm"]) == (rec, arm)]
        present = [c for c in sample if c["modes"]]
        for field in ("top1_A", "best5_A", "oracle_A", "regret_A"):
            groups.append(dict(receiver=rec, arm=arm, measure=field, planned=len(sample), empty=len(sample)-len(present),
                               top1_gt5=sum(c["serious"] for c in present), **distribution([c[field] for c in present])))
    table(result / "groups.csv", groups)
    index = {(c["receiver"], c["arm"], c["seed"]): c for c in cells}
    pairs, summaries = [], []
    decisions = [("4EY7", arm, "4EY6", arm) for arm in ARMS]
    decisions += [(r, a, r, "pro24k") for r in RECEIVERS for a in ("pro48k","official24k")]
    decisions += [("4EY7","pro24k_mapped","4EY7","pro24k"), ("4EY7","pro24k_mapped","4EY6","pro24k")]
    for r1, a1, r0, a0 in decisions:
        contrasts = []
        for seed in SEEDS:
            c1, c0 = index[r1, a1, seed], index[r0, a0, seed]
            if not c1["modes"] or not c0["modes"]:
                continue
            values = {k: c1[k]-c0[k] for k in ("top1_A","oracle_A","best5_A","regret_A")}
            label = dict(receiver1=r1, arm1=a1, receiver0=r0, arm0=a0)
            pairs.append(dict(**label, seed=seed, **values))
            contrasts.append(values)
        if contrasts:
            for k in contrasts[0]:
                summaries.append(dict(**label, measure=k, complete=len(contrasts), **distribution([c[k] for c in contrasts])))
    table(result / "paired.csv", pairs)
    table(result / "paired_summary.csv", summaries)
    # Distinct pose identities within a receiver, irrespective of the scores used to publish them.
    # Compare all 18 heavy-atom coordinates in fixed identity order, allowing only export rounding.
    repeats = []
    for r in RECEIVERS:
        for seed in SEEDS:
            low = all_xyz[r, "pro24k", seed]
            if low:
                for arm in ("pro48k","official24k", *(("pro24k_mapped",) if r == "4EY7" else ())):
                    other = all_xyz[r, arm, seed]
                    if other:
                        d = np.sqrt(np.sum((np.array(other)-low[0])**2, axis=2).mean(axis=1))
                        repeats.append(dict(receiver=r, seed=seed, source="pro24k_top1", compared=arm,
                                            closest_rank=int(np.argmin(d)+1), rmsd_A=d.min(), identical=(d.min()==0)))
    table(result / "coordinate_repeats.csv", repeats)
    save_json(result / "evaluation_receipt.json", dict(direct_protein_fits=fits, total_wall_cpu_s=(tick()-start).tolist()))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("result", type=Path)
    report(p.parse_args().result)

"""Reproduce the map search / final-score comparison on frozen macrocycle inputs."""
import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path
import re
import subprocess
import sys

import numpy as np

BENCH = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(BENCH.parent))
from benchmarks.concerted_ring.run import command, timed
from external_score import poses, reference
from external_run import write


def closure_energy(distances, continued, table_limit):
    # Production minimization linearly interpolates sampled values in r^2.
    # Input here is rounded PDBQT, so this independent reconstruction is approximate.
    total = 0.
    for r in distances:
        x = 32 * r**2
        if not continued and r >= 20:
            continue
        if continued and r >= table_limit:
            value = 50 * r
            total += 1000 * value / (1000 + value)
            continue
        i = int(x)
        a, b = (50 * np.sqrt(k / 32) if continued or k < 32 * 20**2 else 0. for k in (i, i+1))
        value = a + (x - i) * (b - a)
        total += 1000 * value / (1000 + value)
    return total


def run_case(binary, inputs, output, pdb, score, stage, seed, trajectories, cap, modes):
    casename = f"{pdb}-{score}-{stage}-{seed}-{trajectories}-{cap}"
    directory = inputs / pdb
    box = output / (pdb + ".box")
    center = json.loads((directory / "preparation.json").read_text())["center"]
    box.write_text("".join(f"center_{d} = {v}\nsize_{d} = 26\n" for d, v in zip("xyz", center)))
    pose = output / "raw" / (casename + ".pdbqt")
    log = output / "raw" / (casename + ".log")
    args = command((pdb, "", directory / "receptor.pdbqt", directory / "ligand.pdbqt", box),
                   seed, cap, binary, pose)
    args += ["--scoring", score]
    if stage == "map":
        args.append("--no_refine")
    for option, value in (("--exhaustiveness", trajectories), ("--num_modes", modes),
                          ("--verbosity", 2), ("--energy_range", 1e7)):
        args[args.index(option) + 1] = str(value)
    elapsed, cpu, _ = timed(args, log)
    return dict(cell=casename, pdb=pdb, scoring=score, stage=stage, seed=seed,
                trajectories=trajectories, cap_per_trajectory=cap, requested_modes=modes,
                cpu_s=round(cpu, 6), wall_s=round(elapsed / 1e9, 6))


def analyze(case, raw, ref, continued):
    key = case["cell"]
    path = raw / (key + ".pdbqt")
    listing = path.read_text()
    log = (raw / (key + ".log")).read_text()
    coords = list(poses(path))
    reported = list(map(float, re.findall(r"^REMARK VINA RESULT:\s+(\S+)", listing, re.M)))
    unscaled = list(map(float, re.findall(r"^REMARK INTER \+ INTRA:\s+(\S+)", listing, re.M)))
    grid_energy = list(map(float, re.findall(r"^ENERGY FROM SEARCH:\s+(\S+)", log, re.M)))
    use_grid = case["stage"] == "map"
    if not (len(coords) == len(reported) == len(unscaled)):
        raise ValueError(f"{key}: missing coordinates or scores")
    if use_grid and len(coords) != len(grid_energy):
        raise ValueError(f"{key}: cannot match all search energies to saved poses")
    counts = re.findall(r"^Search trajectory: .*?hunt_evals=(\d+) .*?second_evals=(\d+)", log, re.M)
    if len(counts) != int(case["trajectories"]):
        raise ValueError(f"{key}: trajectory count mismatch")
    explicit = int(re.search(r"^Search counters: .*?explicit_evals=(\d+)", log, re.M)[1])
    atom_ids, crystal, sdf_order, sdf, ends, partners, _, crystal_bond = ref
    score_order = sorted(range(len(coords)), key=reported.__getitem__)
    records = []
    for i, pose in enumerate(coords):
        xyz = np.array([pose[a][1] for a in atom_ids])
        rmsd = np.sqrt(((xyz - crystal)**2).sum(axis=-1).mean(axis=-1)).min()
        sdf_xyz = np.array([pose[a][1] for a in sdf_order])
        sdf_rmsd = np.sqrt(((sdf_xyz - sdf)**2).sum(axis=-1).mean(axis=-1)).min()
        if abs(rmsd - sdf_rmsd) > 0.0015:
            raise ValueError(f"{key}: crystal / independent SDF mapping discrepancy")
        bond = np.linalg.norm(pose[ends[0]][1] - pose[ends[1]][1])
        distances = [np.linalg.norm(pose[a][1] - pose[b][1]) for a, b in partners.items()]
        glue = closure_energy(distances, continued, 20.)
        delta = grid_energy[i] - unscaled[i] if use_grid else None
        records.append(dict(cell=key, saved_rank=i+1, reported_score_rank=score_order.index(i)+1,
                            reported_score=reported[i], final_inter_plus_intra=unscaled[i],
                            symmetry_RMSD_A=round(rmsd, 6), sdf_RMSD_A=round(sdf_rmsd, 6),
                            closure_bond_A=round(bond, 6), crystal_bond_A=round(crystal_bond, 6),
                            max_marker_error_A=round(max(distances), 6), severe_bond=int(bond > 3),
                            rounded_pose_glue=round(glue, 5), search_energy=grid_energy[i] if use_grid else "",
                            search_minus_final=round(delta, 5) if use_grid else "",
                            residual_after_glue=round(delta-glue, 5) if use_grid else ""))
    inversions = [(a,b) for a,b in itertools.combinations(range(len(coords)),2) if reported[a] > reported[b]]
    oracle = min((r["symmetry_RMSD_A"] for r in records), default="")
    row = dict(**case, modes=len(records), inversions=len(inversions),
               inversions_with_severe=sum(bool(records[a]["severe_bond"] or records[b]["severe_bond"])
                                           for a,b in inversions),
               oracle_A=oracle, severe_all=sum(r["severe_bond"] for r in records),
               hunt_evals=sum(int(a) for a,b in counts), save_evals=sum(int(b) for a,b in counts),
               explicit_evals=explicit)
    for method, indices in (("saved", list(range(len(records)))), ("score", score_order)):
        top5 = [records[i] for i in indices[:5]]
        row.update({f"{method}_top1_A":top5[0]["symmetry_RMSD_A"] if top5 else "",
                    f"{method}_regret_A":round(top5[0]["symmetry_RMSD_A"]-oracle,6) if top5 else "",
                    f"{method}_best5_A":min((r["symmetry_RMSD_A"] for r in top5), default=""),
                    f"{method}_severe_top5":sum(r["severe_bond"] for r in top5)})
    return row, records


def acyclic_controls(output, binary=None):
    basic = BENCH.parent / "examples" / "basic"
    rows = []
    for score, seed in itertools.product(("vina", "vinardo"), (191009, 458921)):
        key = f"1IEP-{score}-map-{seed}"
        path = output / "raw" / (key + ".pdbqt")
        log = output / "raw" / (key + ".log")
        if binary is not None:
            case = ("1IEP", "", basic / "1iep_receptor.pdbqt", basic / "1iep_ligand.pdbqt",
                    basic / "1iep_receptor.box.txt")
            cmd = command(case, seed, 12000, binary, path)
            cmd.extend(["--scoring", score, "--no_refine"])
            for option, value in (("--exhaustiveness", 4), ("--num_modes", 20),
                                  ("--verbosity", 2), ("--energy_range", 1e7)):
                cmd[cmd.index(option) + 1] = str(value)
            with log.open("w") as handle:
                subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT, check=True)
        values = list(map(float, re.findall(r"^REMARK VINA RESULT:\s+(\S+)", path.read_text(), re.M)))
        if not values:
            raise ValueError(key + ": missing poses")
        inversions = sum(a > b for a, b in itertools.combinations(values, 2))
        rows.append(dict(cell=key, modes=len(values), inversions=inversions))
    write(output / "acyclic.csv", rows)


def distributions(summaries):
    rows = []
    for stage, trajectories in (("map", 1), ("map", 4), ("explicit", 4)):
        cells = [r for r in summaries if r["stage"] == stage and int(r["trajectories"]) == trajectories]
        for metric in ("saved_top1_A", "score_top1_A", "saved_best5_A", "score_best5_A",
                       "oracle_A", "saved_regret_A", "score_regret_A"):
            values = [r[metric] for r in cells if r[metric] != ""]
            q = np.percentile(values, [25, 50, 75, 90, 95])
            rows.append(dict(stage=stage, trajectories=trajectories, metric=metric,
                             valid_cells=len(values), missing_cells=len(cells)-len(values),
                             mean=round(np.mean(values), 6),
                             **{f"p{k}":round(v, 6) for k, v in zip((25, 50, 75, 90, 95), q)}))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--binary", type=Path, help="Run docking first; omission reanalyzes archived output")
    parser.add_argument("--continued", action="store_true", help="Use the untruncated closure rule for reconstruction")
    args = parser.parse_args()
    runfile = args.output / "runs.csv"
    if args.binary:
        (args.output / "raw").mkdir(parents=True, exist_ok=True)
        binary = args.binary.resolve()
        cases = []
        for pdb, seed in itertools.product(("3VHA", "6QBH", "8H6P"), (191009, 458921)):
            for cap in (3000, 12000):
                cases.append(run_case(binary, args.inputs, args.output, pdb, "vina", "map", seed, 1, cap, 5))
            for score, stage in itertools.product(("vina", "vinardo"), ("map", "explicit")):
                cases.append(run_case(binary, args.inputs, args.output, pdb, score, stage, seed, 4, 12000, 20))
        write(runfile, cases)
        (args.output / "binary.sha256").write_text(hashlib.sha256(binary.read_bytes()).hexdigest() + "\n")
    acyclic_controls(args.output, args.binary.resolve() if args.binary else None)
    with runfile.open() as stream:
        cases = list(csv.DictReader(stream))
    refs = {p: reference(args.inputs, p) for p in sorted({c["pdb"] for c in cases})}
    summaries, poses_table = [], []
    for case in cases:
        summary, records = analyze(case, args.output / "raw", refs[case["pdb"]], args.continued)
        summaries.append(summary)
        poses_table.extend(records)
    write(args.output / "cells.csv", summaries)
    write(args.output / "poses.csv", poses_table)
    write(args.output / "distributions.csv", distributions(summaries))


if __name__ == "__main__":
    main()

"""Replay same-map AD4 docking for all four archived side-chain partitions (Linux)."""
import argparse
from collections import Counter
import itertools
import json
import os
from pathlib import Path
import re
import subprocess
import sys
from time import monotonic

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "cubic_grid" / "fkbp_ad4"))
sys.path.insert(0, str(HERE.parent / "cubic_grid" / "flex_replay"))
import run_official
from run_flex import ligand_only, measure, reference, write_csv

GROUPS = {"FABP4": ("2NNQ", "3FR2"), "RBP4": ("1RBP", "3FMZ")}
SEEDS = (62483, 87911, 120997)


def atoms(text):
    return [line for line in text.splitlines() if line.startswith(("ATOM  ", "HETATM"))]


def coordinates(line):
    return [float(line[i:i + 8]) for i in (30, 38, 46)]


def prepare_maps(inputs, flex, out, autogrid):
    rows = []
    for group in GROUPS.values():
        for pdb in group:
            folder = out / "partitions" / pdb
            folder.mkdir(parents=True, exist_ok=True)
            rigid = flex / f"{pdb}_rigid.pdbqt"
            side = flex / f"{pdb}_flex.pdbqt"
            root = inputs / pdb
            def records(path):
                return Counter(line[12:26] + line[30:54] + line[70:] for line in atoms(path.read_text()))
            if records(rigid) & records(side) or records(rigid) + records(side) != records(root / "receptor.pdbqt"):
                raise ValueError(f"Invalid receptor partition for {pdb}")
            for source, name in ((rigid, "receptor.pdbqt"), (root / "preparation.json", "preparation.json")):
                (folder / name).write_bytes(source.read_bytes())
            types = {line[77:].strip() for path in [side] + [inputs / q / "ligand.pdbqt" for q in group]
                     for line in atoms(path.read_text())}
            run_official.TYPES = sorted(types)
            started = monotonic()
            hashes = run_official.prep_grid(out / "partitions", out, autogrid, pdb)
            rows.append(dict(pdb=pdb, seconds=monotonic()-started, maps=hashes, types=sorted(types),
                             side_sha256=run_official.sha(side), rigid_sha256=run_official.sha(rigid)))
    (out / "maps.json").write_text(json.dumps(rows, indent=2) + "\n")


def run_timed(command, logfile):
    with logfile.open("w") as log:
        start = monotonic()
        proc = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
        _, status, stats = os.wait4(proc.pid, 0)
        proc.returncode = os.waitstatus_to_exitcode(status)
        return dict(status=proc.returncode, wall_s=monotonic()-start,
                    cpu_s=stats.ru_utime+stats.ru_stime, peak_kib=stats.ru_maxrss)


def dock_ad4(binary, ligand, maps, flex, seed, pose, counters):
    command = [binary.resolve(), "--ligand", ligand, "--maps", maps, "--flex", flex, "--scoring", "ad4",
               "--cpu", "1", "--exhaustiveness", "1", "--max_evals", "12000", "--num_modes", "9",
               "--energy_range", "100", "--seed", str(seed), "--out", pose]
    if counters:
        command.append("--search_counters")
    return run_timed(command, pose.with_suffix(".log"))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("inputs", "flex", "out", "autogrid", "old", "new"):
        p.add_argument(arg, type=Path)
    p.add_argument("--official", type=Path)
    p.add_argument("--maps-only", action="store_true")
    p.add_argument("--existing-maps", action="store_true")
    a = p.parse_args()
    inputs, flex, out = [x.resolve() for x in (a.inputs, a.flex, a.out)]
    pools = out / "pools"
    pools.mkdir(parents=True, exist_ok=True)
    if not a.existing_maps:
        prepare_maps(inputs, flex, out, a.autogrid.resolve())
    if a.maps_only:
        return
    arms = {"old": a.old, "new": a.new}
    if a.official:
        arms["official"] = a.official
    (out / "binaries.json").write_text(json.dumps({k:run_official.sha(v) for k,v in arms.items()}, indent=2)+"\n")
    cases, modes = [], []
    for target, group in GROUPS.items():
        for query, receiver, seed in itertools.product(group, group, SEEDS):
            get = lambda pdb: json.loads((inputs / pdb / "preparation.json").read_text())
            q, r = get(query), get(receiver)
            expected = reference(q, r)
            side = flex / f"{receiver}_flex.pdbqt"
            fixed = [line for line in atoms(side.read_text()) if line[12:16].strip() in ("CA", "CB")]
            anchors = np.asarray([coordinates(line) for line in fixed])
            assert len(anchors) == 2
            for arm in (list(arms) if seed != SEEDS[1] else list(reversed(arms))):
                stem = pools / f"{query}-{receiver}-{seed}-{arm}"
                pose, log = stem.with_suffix(".pdbqt"), stem.with_suffix(".log")
                case = dict(target=target, query=query, receiver=receiver, seed=seed, arm=arm)
                cost = dock_ad4(arms[arm], inputs/query/"ligand.pdbqt", out/"maps"/receiver/"receptor",
                                side, seed, pose, arm != "official")
                trace = log.read_text()
                text = pose.read_text() if pose.exists() and cost["status"] == 0 else ""
                if cost["status"] and "WARNING: Zero poses" not in trace:
                    raise RuntimeError(trace)
                rmsd = []
                if "ENDMDL" in text:
                    ligand = ligand_only(text)
                    only = stem.with_suffix(".ligand")
                    only.write_text(ligand)
                    rmsd = measure(only, q["serial"], expected)
                    only.unlink()
                    for rank,(err, block) in enumerate(zip(rmsd, ligand.split("ENDMDL")), 1):
                        xyz = np.array([coordinates(line) for line in atoms(block) if int(line[6:11]) in q["serial"]])
                        nearest = np.linalg.norm(xyz[:,None,:]-anchors[None,:,:], axis=2).min()
                        energy = float(re.search(r"REMARK VINA RESULT:\s*([-+.\d]+)", block)[1])
                        modes.append(dict(**case, rank=rank, rmsd=err, min_anchor_A=nearest, reported_score=energy))
                counts = {key: (int(m[1]) if (m:=re.search(r"\b"+key+r"=(\d+)",trace)) else "")
                          for key in ("hunt_evals", "second_evals", "explicit_evals")}
                cases.append(dict(**case, **cost, **counts, reported=len(rmsd), top1=rmsd[0] if rmsd else "",
                                  best5=min(rmsd[:5]) if rmsd else "", best_reported=min(rmsd) if rmsd else "",
                                  regret=rmsd[0]-min(rmsd) if rmsd else ""))
    write_csv(out/"cases.csv", cases)
    write_csv(out/"modes.csv", modes)


if __name__ == "__main__":
    main()

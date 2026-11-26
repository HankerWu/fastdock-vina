"""Equal-search-budget dry, mixed and hydrated pools with common-coordinate scoring."""
import argparse
import csv
import hashlib
import json
import re
import subprocess
import tarfile
import tempfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
SEEDS = (95701, 95713, 95717)
ARMS = ("dry", "polar")
PLANS = {"dry+dry": ("dry", "dry"), "dry+polar": ("dry", "polar"), "polar+polar": ("polar", "polar")}


def write(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def unpack(source, dest):
    dest.mkdir(parents=True, exist_ok=True)
    with tarfile.open(source) as archive:
        archive.extractall(dest, filter="data")


def sources(work):
    bptf, brdt = work / "BPTF", work / "BRDT"
    unpack(BASE / "bptf_water/raw.tar.xz", bptf)
    unpack(BASE / "pose/inputs.tar.xz", bptf / "base")
    unpack(BASE / "water_transfer/inputs.tar.xz", brdt)
    unpack(BASE / "water_transfer/outputs.tar.xz", brdt)
    cells = read(BASE / "bptf_water/confirm_cells.csv") + read(BASE / "water_transfer/cells.csv")
    keep = [r for r in cells if r["arm"] in ARMS and int(r["seed"]) in SEEDS and int(r["max_evals"]) == 120000]
    modes = read(BASE / "bptf_water/confirm_modes.csv") + read(BASE / "water_transfer/modes.csv")
    by_source = defaultdict(list)
    for row in modes:
        by_source[row["target"], row["query"], row["receiver"], row["arm"], row["seed"]].append(row)
    return keep, by_source


def score_cell(cell, binary, work, by_source):
    target, q, r, arm, seed = (cell[k] for k in ("target", "query", "receiver", "arm", "seed"))
    directory = work / target
    if target == "BPTF":
        folder = "confirm"
        ligdir = directory / ("base" if q in ("7F5D", "7RWO") else "external") / q
        recdir = directory / "base" / r
    else:
        folder = "poses"
        ligdir = directory / "inputs" / q
        recdir = directory / "inputs" / r
    pdbqt = directory / arm / folder / f"{q}-{r}-{seed}-{arm}.pdbqt"
    assert hashlib.sha256(pdbqt.read_bytes()).hexdigest() == cell["pdbqt_sha256"]
    blocks = re.findall(r"^MODEL[^\n]*\n(.*?)^ENDMDL$", pdbqt.read_text(), re.M | re.S)
    source_rows = by_source[target, q, r, arm, seed]
    assert len(blocks) == len(source_rows) == int(cell["modes"]) <= 9
    center = json.loads((recdir / "metadata.json").read_text())["center"]
    records = []
    with tempfile.TemporaryDirectory() as tmp:
        pose = Path(tmp) / "pose.pdbqt"
        for i, (block, source) in enumerate(zip(blocks, source_rows), 1):
            assert int(source["rank"]) == i
            pose.write_text(block)
            row = dict(target=target, query=q, receiver=r, kind=cell["kind"], seed=seed, source_arm=arm, rank=i,
                       rmsd_A=float(source["rmsd_A"]), reported_score=float(source["score"]))
            for scoring in ARMS:
                cmd = [binary, "--score_only", "--receptor", directory / scoring / r / "receptor.pdbqt",
                       "--ligand", ligdir / "ligand.pdbqt", "--pose_coordinates", pose, "--unbound_energy", "0", "--cpu", "1"]
                for axis, x in zip("xyz", center):
                    cmd += ["--center_" + axis, str(x), "--size_" + axis, "24"]
                output = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
                def value(key):
                    m = re.search(key + r"\s*:\s*([-+0-9.eE]+) \(kcal/mol\)", output)
                    if not m:
                        raise ValueError((q, r, scoring, key, output))
                    return float(m[1])
                assert value("Unbound System's Energy") == 0
                row[scoring + "_common_score"] = value("Estimated Free Energy of Binding")
                row[scoring + "_inter"] = value("Final Intermolecular Energy")
                row[scoring + "_intra"] = value("Final Total Internal Energy")
            assert row["dry_intra"] == row["polar_intra"]
            records.append(row)
    return records


def evaluate(cells, energies):
    source = {(c["target"], c["query"], c["receiver"], int(c["seed"]), c["arm"]): c for c in cells}
    poses = defaultdict(list)
    groups = {}
    for p in energies:
        key = p["target"], p["query"], p["receiver"]
        groups[key] = p["kind"]
        poses[*key, int(p["seed"]), p["source_arm"]].append(p)
    folds = []
    for (target, q, r), kind in sorted(groups.items()):
        for i, first in enumerate(SEEDS):
            second = SEEDS[(i + 1) % len(SEEDS)]
            for plan, (a, b) in PLANS.items():
                keys = [(target, q, r, first, a), (target, q, r, second, b)]
                merged = poses[keys[0]] + poses[keys[1]]
                assert len(merged) == sum(int(source[k]["modes"]) for k in keys)
                for scoring in ARMS:
                    best = min(p[scoring + "_common_score"] for p in merged)
                    winners = [p for p in merged if p[scoring + "_common_score"] == best]
                    top = winners[0]  # fixed first source, then native rank for exact reported ties
                    folds.append(dict(target=target, kind=kind, query=q, receiver=r, first_seed=first, second_seed=second,
                                      plan=plan, scoring=scoring, candidates=len(merged),
                                      top_source=top["source_arm"], top_source_seed=top["seed"], top_source_rank=top["rank"],
                                      top1_A=top["rmsd_A"], oracle_A=min(p["rmsd_A"] for p in merged), score=best,
                                      tied_top=len(winners), ambiguous_hit=len({p["rmsd_A"] <= 2 for p in winners}) > 1,
                                      search_evals=sum(int(source[k]["search_evals"]) for k in keys),
                                      explicit_evals=sum(int(source[k]["explicit_evals"]) for k in keys),
                                      source_wall_s=sum(float(source[k]["wall_s"]) for k in keys)))
    grouped = defaultdict(list)
    for f in folds:
        grouped[f["target"], f["kind"], f["plan"], f["scoring"]].append(f)
    summary = []
    for (target, kind, plan, scoring), c in sorted(grouped.items()):
        summary.append(dict(target=target, kind=kind, plan=plan, scoring=scoring, folds=len(c),
                            top_hits=sum(r["top1_A"] <= 2 for r in c), oracle_hits=sum(r["oracle_A"] <= 2 for r in c),
                            mean_candidates=sum(r["candidates"] for r in c) / len(c),
                            mean_search_evals=sum(r["search_evals"] for r in c) / len(c),
                            mean_explicit_evals=sum(r["explicit_evals"] for r in c) / len(c),
                            mean_source_wall_s=sum(r["source_wall_s"] for r in c) / len(c),
                            ambiguous_ties=sum(r["ambiguous_hit"] for r in c)))
    return folds, summary


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("binary", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--jobs", type=int, default=4)
    args = p.parse_args()
    binary = args.binary.resolve()
    args.output.mkdir(exist_ok=True, parents=True)
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        cells, modes = sources(work)
        with ThreadPoolExecutor(args.jobs) as pool:
            energies = [record for rows in pool.map(lambda c: score_cell(c, binary, work, modes), cells) for record in rows]
    folds, summary = evaluate(cells, energies)
    for name, rows in (("energies", energies), ("folds", folds), ("summary", summary)):
        write(args.output / (name + ".csv"), rows)

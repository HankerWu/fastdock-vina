"""External assay transfer of the prespecified full and historical strain arms."""

import argparse
import itertools
import statistics
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from pose import HERE, METHODS, ROOT, physical, read, score_pose, unpack, write


def run(job):
    binary, work, comp = job
    receptor0 = work / "2P16"
    pool = physical.read_poses(receptor0 / (comp + "_271828.pdbqt"))
    mol, identity = physical.decode_ligand(pool[0])
    free, n = physical.gas_reference(mol)
    scale = physical.evaluate_coordinates(
        binary, receptor0 / (comp + ".pdbqt"), receptor0 / "2P16.pdbqt", None, 90
    )["scale"]
    rows = []
    for rec, seed in itertools.product(("2P16", "1FJS"), (271828, 314159)):
        path = work / rec / (comp + f"_{seed}.pdbqt")
        poses = physical.read_poses(path)
        physical.fixed_coordinate_inputs(poses, work / rec / (comp + ".pdbqt"), None)
        for p in poses:
            rows.append(
                dict(
                    compound=comp,
                    receptor=rec,
                    seed=seed,
                    **score_pose(p, identity, free, n, scale),
                )
            )
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("binary", type=Path)
    p.add_argument("work", type=Path)
    p.add_argument("--jobs", type=int, default=1)
    a = p.parse_args()
    unpack(HERE.parent / "affinity/fxa_frozen.tar.xz", a.work)
    folder = ROOT / "benchmarks/affinity_partition"
    labels = {
        r["ligand"]: r for r in read(folder / "fxa_series.csv") if r["included"] == "1"
    }
    with ProcessPoolExecutor(a.jobs) as pool:
        rows = [
            r
            for v in pool.map(run, [(a.binary.resolve(), a.work, c) for c in labels])
            for r in v
        ]
    write(a.work / "modes.csv", rows)
    summarize(a.work, folder, labels, rows)


def summarize(work, folder, labels, rows):
    sels = {}
    for key, v in itertools.groupby(
        rows, key=lambda r: (r["compound"], r["receptor"], int(r["seed"]))
    ):
        states = list(v)
        for arm in METHODS:
            sels[*key, arm, "all"] = min(float(s[arm]) for s in states)
            sels[*key, arm, "first"] = float(
                next(s for s in states if int(s["rank"]) == 1)[arm]
            )
    for r in read(folder / "results/fxa_rotor_cells.csv"):
        if r["ligand"] in labels and sels[
            r["ligand"], r["receptor"], int(r["seed"]), "native", "all"
        ] != float(r["native"]):
            raise ValueError("Native reference diverges")
    edges = []
    for stronger, weaker in itertools.combinations(
        sorted(labels, key=lambda c: float(labels[c]["ki_nm"])), 2
    ):
        a0, b = labels[stronger], labels[weaker]
        if float(b["ki_nm"]) / float(a0["ki_nm"]) < 3:
            continue
        equal = (a0["charge"], a0["heavy"]) == (b["charge"], b["heavy"])
        for receiver in ("2P16", "1FJS", "both"):
            recs = ("2P16", "1FJS") if receiver == "both" else (receiver,)
            for method, policy in itertools.product(METHODS, ("first", "all")):
                diffs = [
                    sels[weaker, rec, s, method, policy]
                    - sels[stronger, rec, s, method, policy]
                    for rec, s in itertools.product(recs, (271828, 314159))
                ]
                gap = statistics.median(diffs)
                edges.append(
                    {
                        "receptor": receiver,
                        "stronger": stronger,
                        "weaker": weaker,
                        "same_heavy_charge": equal,
                        "method": method,
                        "policy": policy,
                        "gap": gap,
                        "credit": 0.5 if abs(gap) < 1e-10 else float(gap > 0),
                    }
                )
    write(work / "pairs.csv", edges)
    summary = []
    for rec, arm, policy, subset in itertools.product(
        ("2P16", "1FJS", "both"),
        METHODS,
        ("first", "all"),
        ("all", "equal_heavy_charge"),
    ):
        rr = [
            r
            for r in edges
            if r["method"] == arm
            and r["receptor"] == rec
            and r["policy"] == policy
            and (subset == "all" or r["same_heavy_charge"])
        ]
        summary.append(
            {
                "receptor": rec,
                "method": arm,
                "policy": policy,
                "stratum": subset,
                "pairs": len(rr),
                "credit": sum(r["credit"] for r in rr),
            }
        )
    write(work / "summary.csv", summary)


if __name__ == "__main__":
    main()

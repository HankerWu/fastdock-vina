"""Apply the frozen scorer to the separate previously assayed 28-ligand human FXa series."""

import argparse
import hashlib
import itertools
import statistics
import tarfile
from pathlib import Path
from scipy.spatial import cKDTree

import test_entropy as probe

HERE = Path(__file__).resolve().parent
LABELS = HERE.parents[1] / "affinity_partition"
RECS = ("2P16", "1FJS")
SEEDS = (271828, 314159)


def pack(work, output):
    with tarfile.open(output, "w:xz") as tar:
        for state, folder in zip(RECS, work):
            for p in sorted(folder.glob("*.pdbqt")):
                if p.name == f"{state}.pdbqt" or p.name.startswith("CHEMBL"):
                    tar.add(p, arcname=f"{state}/{p.name}")


def run(archive, dest):
    labels = {
        r["ligand"]: r
        for r in probe.read_csv(LABELS / "fxa_series.csv")
        if r["included"] == "1"
    }
    native = {
        (r["receptor"], r["ligand"], int(r["seed"])): float(r["native"])
        for r in probe.read_csv(LABELS / "results" / "fxa_rotor_cells.csv")
    }
    rotors = {
        r["compound"]: float(r["num_tors"])
        for r in probe.read_csv(HERE / "fxa_native_rotors.csv")
    }
    mode_rows, cells, edges, tables = [], [], [], []
    for state in RECS:
        inputs = dict(probe.contents(archive, state + "/"))
        receptor = cKDTree(
            [
                a[2]
                for a in map(probe.atom, probe.atom_lines(inputs[state]))
                if probe.heavy(a)
            ]
        )
        for c in labels:
            topology = probe.ligand_graph(inputs[c])
            for seed in SEEDS:
                pool = inputs[f"{c}_{seed}"].split("ENDMDL")
                models = [v for v in pool if v.strip()]
                modes = [
                    dict(rank=i + 1, **probe.pose(v, receptor, topology, rotors[c]))
                    for i, v in enumerate(models)
                ]
                assert modes[0]["native"] == native[state, c, seed]
                mode_rows.extend(
                    dict(receptor=state, compound=c, seed=seed, **r) for r in modes
                )
                for policy, score in itertools.product(("first", "all"), probe.SCORES):
                    best = min(
                        modes if policy == "all" else modes[:1],
                        key=lambda x: (x[score], x["rank"]),
                    )
                    cells.append(
                        dict(
                            receptor=state,
                            compound=c,
                            seed=seed,
                            policy=policy,
                            score=score,
                            rank=best["rank"],
                            value=best[score],
                        )
                    )
    vals = {
        (r["receptor"], r["compound"], r["seed"], r["policy"], r["score"]): r["value"]
        for r in cells
    }
    for receptor in (*RECS, "both"):
        chosen = RECS if receptor == "both" else [receptor]
        for policy, score in itertools.product(("first", "all"), probe.SCORES):
            group = []
            for a, b in itertools.combinations(
                sorted(labels, key=lambda c: float(labels[c]["ki_nm"])), 2
            ):
                fold = float(labels[b]["ki_nm"]) / float(labels[a]["ki_nm"])
                if fold < 3:
                    continue
                # Match the pre-existing FXa primary estimator: median of paired differences across receptor/seed.
                gap = round(
                    statistics.median(
                        vals[s, b, k, policy, score] - vals[s, a, k, policy, score]
                        for s in chosen
                        for k in SEEDS
                    ),
                    10,
                )
                group.append(
                    dict(
                        receptor=receptor,
                        policy=policy,
                        score=score,
                        stronger=a,
                        weaker=b,
                        fold=fold,
                        gap=gap,
                        credit=1 if gap > 0 else 0 if gap < 0 else 0.5,
                        size_charge_match=int(
                            all(
                                labels[a][x] == labels[b][x]
                                for x in ("heavy", "charge")
                            )
                        ),
                    )
                )
            edges.extend(group)
            for subset in ("all", "size_charge_match"):
                use = (
                    group
                    if subset == "all"
                    else [r for r in group if r["size_charge_match"]]
                )
                tables.append(
                    dict(
                        receptor=receptor,
                        policy=policy,
                        score=score,
                        subset=subset,
                        pairs=len(use),
                        correct=sum(r["credit"] for r in use),
                    )
                )
    old = probe.read_csv(LABELS / "results" / "size_controls.csv")
    for r in tables:
        if r["receptor"] == "both" and r["score"] == "native":
            (ref,) = [
                v
                for v in old
                if v["target"] == "fxa"
                and v["subset"]
                == ("all" if r["subset"] == "all" else "equal_size_charge")
            ]
            assert float(ref["native"]) == r["correct"]
    for name, rows in [
        ("fxa_poses.csv.xz", mode_rows),
        ("fxa_selected.csv.xz", cells),
        ("fxa_pairs.csv.xz", edges),
        ("fxa_summary.csv", tables),
    ]:
        probe.write(dest / name, rows)
    print(
        "FXa archive",
        hashlib.sha256(archive.read_bytes()).hexdigest(),
        len(mode_rows),
        "poses",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pack", type=Path, nargs=2, metavar=("2P16_work", "1FJS_work")
    )
    parser.add_argument("--archive", type=Path, default=HERE / "fxa_frozen.tar.xz")
    parser.add_argument("--out", type=Path, default=HERE)
    args = parser.parse_args()
    if args.pack:
        pack(args.pack, args.archive)
    run(args.archive, args.out)

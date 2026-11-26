"""Score the deposited affinity poses. Coordinates, assay rows and model order stay frozen."""

import argparse
import csv
import hashlib
import itertools
import lzma
import math
import tarfile
from collections import defaultdict
from pathlib import Path
from statistics import median

import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parents[1] / "geometry_comparison" / "affinity"
COST = 0.00198720425864 * 298.15 * math.log(3)
SCORES = ("native", "raw", "uniform", "contact")
CAPS = (12000, 60000)


def read_csv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def write(path, rows):
    opener = lzma.open if path.suffix == ".xz" else open
    with opener(path, "wt", newline="") as stream:
        out = csv.DictWriter(stream, fieldnames=rows[0])
        out.writeheader()
        out.writerows(rows)


def contents(archive, prefix):
    with tarfile.open(archive) as tar:
        for m in tar:
            if m.isfile() and m.name.startswith(prefix) and m.name.endswith(".pdbqt"):
                yield Path(m.name).stem, tar.extractfile(m).read().decode()


def atom(line):
    return (
        int(line[6:11]),
        line[77:].strip(),
        [float(line[i : i + 8]) for i in (30, 38, 46)],
    )


def atom_lines(block):
    return (
        line for line in block.splitlines() if line.startswith(("ATOM  ", "HETATM"))
    )


def heavy(at):
    return at[1] not in {"H", "HD", "HS"}


def ligand_graph(block):
    """Associate each BRANCH with the full downstream heavy subtree (cut both sides)."""
    serials, branches, stack = [], {}, []
    for line in block.splitlines():
        fields = line.split()
        if line.startswith("BRANCH"):
            bond = tuple(map(int, fields[1:]))
            assert bond not in branches
            stack.append(bond)
            branches[bond] = set()
        elif line.startswith("ENDBRANCH"):
            assert tuple(map(int, fields[1:])) == stack.pop()
        elif line.startswith(("ATOM  ", "HETATM")):
            serial, typ, xyz = atom(line)
            if heavy((serial, typ, xyz)):
                serials.append(serial)
                for bond in stack:
                    branches[bond].add(serial)
    assert not stack and len(serials) == len(set(serials))
    index = {s: i for i, s in enumerate(serials)}
    cuts = []
    for (a, b), side in branches.items():
        assert b in side and a not in side
        # A terminal donor heteroatom may have only hydrogens beyond the axis.
        if len(side) > 1 and len(serials) - len(side) > 1:
            cuts.append((index[a], index[b], np.array([s in side for s in serials])))
    return serials, cuts


def pose(block, receptor, topology, native_tors):
    fields = dict(
        line[7:].split(":", 1)
        for line in block.splitlines()
        if line.startswith("REMARK ")
        and line.split(":", 1)[0]
        in {"REMARK INTER + INTRA", "REMARK VINA RESULT", "REMARK UNBOUND"}
    )
    native = float(fields["VINA RESULT"].split()[0])
    raw = round(float(fields["INTER + INTRA"]) - float(fields["UNBOUND"]), 3)
    reconstructed = raw / (1 + 0.05846 * native_tors)
    # Each of the two printed raw terms and reported affinity is rounded to 0.001.
    assert abs(reconstructed - native) <= 0.001501, (raw, native, native_tors)
    atoms = [a for a in map(atom, atom_lines(block)) if heavy(a)]
    serials, cuts = topology
    assert [a[0] for a in atoms] == serials
    xyz = np.array([a[2] for a in atoms])
    close = receptor.query(xyz)[0] <= 4.0
    ncontact = 0
    for a, b, side in cuts:
        axis = xyz[b] - xyz[a]
        dist_to_line = np.linalg.norm(
            np.cross(xyz - xyz[a], axis), axis=1
        ) / np.linalg.norm(axis)
        eligible = close & (dist_to_line >= 1.0)
        ncontact += int(np.any(eligible & side) and np.any(eligible & ~side))
    return dict(
        native=native,
        raw=raw,
        uniform=raw + COST * len(cuts),
        contact=raw + COST * ncontact,
        num_tors=native_tors,
        eligible_rotors=len(cuts),
        constrained=ncontact,
        native_rebuilt=round(reconstructed, 7),
    )


def load(source, rotor_file):
    # Verify the published pool/input byte checksums before deriving results.
    for line in (source / "SHA256SUMS").read_text().splitlines():
        digest, rel = line.split()
        assert hashlib.sha256((source / rel).read_bytes()).hexdigest() == digest
    panel = read_csv(source / "panel.csv")
    rotors = {r["compound"]: float(r["num_tors"]) for r in read_csv(rotor_file)}
    graphs, receptors = {}, {}
    for name, block in contents(source / "inputs.tar.xz", "prepared/"):
        if name in rotors:
            graphs[name] = ligand_graph(block)
        elif len(name) == 4:
            rec = [a[2] for a in map(atom, atom_lines(block)) if heavy(a)]
            receptors[name] = cKDTree(rec)
    poses = []
    for cap in CAPS:
        for name, pool in contents(source / f"cap{cap}" / "outputs.tar.xz", "pools/"):
            state, ligand, seed, optimizer = name.split("_")
            models = [b.strip() for b in pool.split("ENDMDL") if b.strip()]
            for rank, block in enumerate(models, 1):
                assert block.splitlines()[0] == f"MODEL {rank}"
                result = pose(block, receptors[state], graphs[ligand], rotors[ligand])
                poses.append(
                    dict(
                        budget=cap,
                        optimizer=optimizer,
                        receptor=state,
                        compound=ligand,
                        seed=int(seed),
                        rank=rank,
                        **result,
                    )
                )
    return panel, poses


def summarize(source, panel, poses):
    pools = defaultdict(list)
    for row in poses:
        pools[
            tuple(
                row[k] for k in ("budget", "optimizer", "receptor", "compound", "seed")
            )
        ].append(row)
    assert len(pools) == 560
    calls = []
    for (cap, optimizer, receptor, compound, seed), rows in sorted(pools.items()):
        for policy in ("first", "all"):
            for score in SCORES:
                selected = min(
                    rows[:1] if policy == "first" else rows,
                    key=lambda x: (x[score], x["rank"]),
                )
                calls.append(
                    dict(
                        budget=cap,
                        optimizer=optimizer,
                        receptor=receptor,
                        compound=compound,
                        seed=seed,
                        policy=policy,
                        score=score,
                        selected=selected["rank"],
                        value=selected[score],
                    )
                )
    groups = defaultdict(list)
    for row in calls:
        groups[
            tuple(
                row[k]
                for k in (
                    "budget",
                    "optimizer",
                    "policy",
                    "score",
                    "receptor",
                    "compound",
                )
            )
        ].append(row["value"])
    preds = {k: median(v) for k, v in groups.items()}
    assert all(len(v) == 2 for v in groups.values())
    by_assay = defaultdict(dict)
    for r in panel:
        by_assay[r["target"], r["assay"]][r["compound"]] = float(r["value_nM"])
    edges, summary = [], []
    for cap in CAPS:
        native_check = read_csv(source / f"cap{cap}" / "summary.csv")
        for (target, assay), labels in by_assay.items():
            recs = sorted({r["receptor"] for r in native_check if r["assay"] == assay})
            for opt, policy, score in itertools.product(
                ("geometry", "armijo"), ("first", "all"), SCORES
            ):
                for receptor in [*recs, "mean_receptors"]:
                    states = recs if receptor == "mean_receptors" else [receptor]
                    values = {
                        c: sum(preds[cap, opt, policy, score, s, c] for s in states)
                        / len(states)
                        for c in labels
                    }
                    key = dict(
                        budget=cap,
                        optimizer=opt,
                        policy=policy,
                        score=score,
                        target=target,
                        assay=assay,
                        receptor=receptor,
                    )
                    credits = []
                    for a, b in itertools.combinations(
                        sorted(labels, key=labels.get), 2
                    ):
                        if labels[b] / labels[a] < 3:
                            continue
                        delta = round(values[b] - values[a], 10)
                        credit = 1 if delta > 0 else 0 if delta < 0 else 0.5
                        credits.append(credit)
                        edges.append(
                            dict(
                                **key,
                                stronger=a,
                                weaker=b,
                                fold=labels[b] / labels[a],
                                gap=delta,
                                credit=credit,
                            )
                        )
                    rank_corr = spearmanr(
                        [labels[c] for c in values], list(values.values())
                    ).statistic
                    summary.append(
                        dict(
                            **key,
                            pairs=len(credits),
                            correct=sum(credits),
                            spearman=round(rank_corr, 6),
                        )
                    )
                    if receptor != "mean_receptors" and score == "native":
                        (prior,) = [
                            r
                            for r in native_check
                            if r["assay"] == assay
                            and r["receptor"] == receptor
                            and r["arm"] == opt
                            and r["fold_min"] == "3"
                        ]
                        assert (len(credits), sum(credits)) == (
                            int(prior["pairs"]),
                            float(prior["correct"]),
                        )
    return calls, edges, summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, default=SOURCE)
    p.add_argument("--out", type=Path, default=HERE)
    args = p.parse_args()
    args.out.mkdir(exist_ok=True, parents=True)
    panel, poses = load(args.source, HERE / "native_rotors.csv")
    calls, edges, summary = summarize(args.source, panel, poses)
    for name, rows in [
        ("poses.csv.xz", poses),
        ("selected.csv.xz", calls),
        ("pairs.csv.xz", edges),
        ("summary.csv", summary),
    ]:
        write(args.out / name, rows)
    print(
        f"{len(poses)} poses, {len(calls)} selected rows, {len(edges)} evaluated pair rows"
    )


if __name__ == "__main__":
    main()

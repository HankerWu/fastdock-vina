"""Matched HA vs carboxylate cross-docks; crystals are read only by the separate evaluator."""

import argparse
import csv
import json
from pathlib import Path
import sys
import tarfile

import numpy as np
from rdkit import Chem

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent))
from pose_runs import digest, invoke, source_record, table
from generate import profiles

BENCH = ROOT / "benchmarks/cubic_grid"
CASES = {"2NNQ": "3FR2", "3FR2": "2NNQ", "3FMZ": "1RBP"}
SEED = 91813


def chemical_graph(path, out):
    text = path.read_text().splitlines()
    mol = Chem.MolFromSmiles(
        next(l.split()[2] for l in text if l.startswith("REMARK SMILES "))
    )
    order = {}
    for l in text:
        if l.startswith("REMARK SMILES IDX "):
            values = list(map(int, l.split()[3:]))
            order.update(zip(values[::2], values[1::2], strict=True))
    edges = [
        (order[b.GetBeginAtomIdx() + 1], order[b.GetEndAtomIdx() + 1])
        for b in mol.GetBonds()
    ]
    for l in text:
        if l.startswith("REMARK H PARENT "):
            values = list(map(int, l.split()[3:]))
            edges.extend(
                (order[i], j) for i, j in zip(values[::2], values[1::2], strict=True)
            )
    records = {int(l[6:11]) for l in text if l.startswith(("ATOM  ", "HETATM"))}
    assert set(j for pair in edges for j in pair) == records
    out.write_text(
        "".join(f"{i} {j}\n" for i, j in sorted(tuple(sorted(pair)) for pair in edges))
    )


def prepare(p):
    for folder in ("raw", "sources", "inputs"):
        (p / folder).mkdir(exist_ok=True, parents=True)
    dest = p / "inputs"
    files = [
        BENCH / "acid_states/acid-pools.tar.xz",
        BENCH / "human_panel/prepared.tar.xz",
    ]
    cases = {}
    with tarfile.open(files[0]) as acid, tarfile.open(files[1]) as prep:
        for q, r in CASES.items():
            prep.extract(r + "/receptor.pdbqt", dest, filter="data")
            center = json.load(prep.extractfile(r + "/preparation.json"))["center"]
            for state in ("HA", "A"):
                tag = q + "-" + state
                lig = dest / (tag + ".pdbqt")
                graph = lig.with_suffix(".bonds")
                acid.extract(lig.name, dest, filter="data")
                chemical_graph(lig, graph)
                cases[tag] = dict(
                    ligand=lig,
                    receptor=dest / r / "receptor.pdbqt",
                    bonds=graph,
                    center=center,
                    size=24,
                )
    (p / "inputs.json").write_text(
        json.dumps(
            dict(
                archives={str(f.relative_to(ROOT)): digest(f) for f in files},
                cells={
                    k: {
                        str(j): digest(v) if isinstance(v, Path) else v
                        for j, v in c.items()
                    }
                    for k, c in cases.items()
                },
            ),
            indent=2,
        )
        + "\n"
    )
    return cases


def run(p, exe, official, counter):
    assert (
        digest(official)
        == "f31f774f723bba7bbe6e9d1c47577020eea9a8da16424284c043d22593570644"
    )
    cases = prepare(p)
    cost = []
    sources = []
    for tag, c in cases.items():
        root = p / "raw" / (tag + "_base")
        invoke(
            p,
            cost,
            root.name,
            exe,
            c,
            SEED,
            env=dict(
                FD_SOURCE_TRAJECTORY=1,
                FD_SOURCE_STEP=1,
                FD_SOURCE_XYZ=root.with_suffix(".xyz"),
                FD_SOURCE_Q=root.with_suffix(".q"),
            ),
        )
        lines = root.with_suffix(".xyz").read_text().splitlines()
        serial = list(map(int, lines[0][1:].split()))
        x = np.array(lines[1].split()[3:], float).reshape(-1, 3)
        q = np.loadtxt(root.with_suffix(".q"))[2:]
        sources.append(source_record(p, tag, c, SEED, 1, 1, x, q, serial))
    table(p / "sources.csv", sources)
    profiles(p / "sources.csv", p / "candidates")
    for row in csv.DictReader((p / "selections.csv").open()):
        name = row["source"].split("_")[0]
        arm = row["arm"]
        stem = name + "_" + arm
        base = p / "raw" / stem
        base.with_suffix(".input").write_text(
            "1 1 "
            + " ".join(row[k] for k in ("qx", "qy", "qz", "qw", "qi", "qj", "qk"))
            + "\n"
        )
        invoke(
            p,
            cost,
            stem,
            exe,
            cases[name],
            SEED,
            env=dict(
                FD_RIGID_INTERVENTION=base.with_suffix(".input"),
                FD_RIGID_LOG=base.with_suffix(".inject"),
                FD_RIGID_COORDS=base.with_suffix(".stages"),
                FD_RIGID_HISTORY=base.with_suffix(".history"),
            ),
        )
    overhead = []
    for tag, c in cases.items():
        invoke(p, cost, tag + "_high", exe, c, SEED, 50000)
        for name, budget in (("official", 25000), ("official-double", 50000)):
            stem = tag + "_" + name
            pred = invoke(p, cost, stem, official, c, SEED, budget, official=True)
            check = invoke(
                p, overhead, stem + "_count", counter, c, SEED, budget, official=True
            )
            assert pred.read_bytes() == check.read_bytes()
            for key in (
                "hunt_calls",
                "hunt_evals",
                "second_calls",
                "second_evals",
                "explicit_calls",
                "explicit_evals",
                "late_prescore",
            ):
                cost[-1][key] = overhead[-1][key]
            assert cost[-1]["hunt_evals"] > 0
    table(p / "cost.csv", cost)
    table(p / "official_counter.csv", overhead)
    (p / "binaries.json").write_text(
        json.dumps(
            {
                k: digest(f)
                for k, f in (
                    ("observer", exe),
                    ("official", official),
                    ("counter", counter),
                )
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    a = argparse.ArgumentParser(description=__doc__)
    for k in ("work", "exe", "official", "counter"):
        a.add_argument(k, type=lambda s: Path(s).resolve())
    args = a.parse_args()
    run(args.work, args.exe, args.official, args.counter)

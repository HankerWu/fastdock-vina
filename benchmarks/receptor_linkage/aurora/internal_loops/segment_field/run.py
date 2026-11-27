"""Replay receptor-field and ligand-graph interventions with the original producer."""

import argparse
import json
import sys
from pathlib import Path

from rdkit import Chem

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import replay
from pose_probe import get_tree
from summarize import read

FRESH = (45101, 45102, 45103, 45104)
SEGMENTS = {"long": range(281, 292), "short": range(303, 306)}
MOVE = {
    "far_both": {"short", "long"},
    "long_only": {"short"},
    "short_only": {"long"},
    "both": set(),
}


def prepare_arms(work, out):
    original = (work / replay.PROTEINS["dry"]).read_bytes()
    full = (work / replay.PROTEINS["1OL7_released"]).read_bytes()
    if not full.startswith(original):
        raise ValueError("Expected unchanged original receptor as a byte prefix")
    records = full[len(original) :].decode("ascii").splitlines(keepends=True)
    proteins = {"dry": work / replay.PROTEINS["dry"]}
    for arm, move in MOVE.items():
        suffix = []
        for line in records:
            residue = int(line[22:26])
            (segment,) = (
                seg for seg, numbers in SEGMENTS.items() if residue in numbers
            )
            if segment in move:
                line = line[:30] + f"{float(line[30:38]) + 100:8.3f}" + line[38:]
            suffix.append(line)
        path = out / (arm + ".pdbqt")
        path.write_bytes(original + "".join(suffix).encode("ascii"))
        proteins[arm] = path
    if proteins["both"].read_bytes() != full:
        raise ValueError("Both-at-original-position arm changed the source")
    proteins.update({arm + "_graph": proteins[arm] for arm in ("far_both", "both")})
    return proteins


def validate_graph(ligand, graph):
    _, serials, sdf_indices, _, _ = get_tree(ligand)
    mapping = dict(zip(sdf_indices, serials, strict=True))
    mol = Chem.MolFromMolFile(str(ligand.with_suffix(".sdf")), removeHs=False)
    expected = set()
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in mapping and b in mapping:
            expected.add(tuple(sorted((mapping[a], mapping[b]))))
    explicit = {
        tuple(sorted(map(int, line.split()))) for line in graph.read_text().splitlines()
    }
    if expected != explicit:
        raise ValueError(f"Explicit edge file differs from SDF for {ligand.name}")
    return {
        "n_edges": len(explicit),
        "pdbqt_sha256": replay.digest(ligand),
        "sdf_sha256": replay.digest(ligand.with_suffix(".sdf")),
        "edges_sha256": replay.digest(graph),
    }


def run(work, out, pro, official):
    out.mkdir(parents=True, exist_ok=False)
    proteins = prepare_arms(work, out)
    graphs = {}
    for drug in replay.geom.REFERENCES:
        ligand, bonds, _, _ = replay.parameters(work, drug)
        graphs[drug] = validate_graph(ligand, bonds)
    (out / "inputs.json").write_text(
        json.dumps(
            {
                "receptors_sha256": {k: replay.digest(v) for k, v in proteins.items()},
                "graphs": graphs,
                "pro_sha256": replay.digest(pro),
                "official_sha256": replay.digest(official),
            },
            indent=2,
        )
        + "\n"
    )
    locked = ("far_both_graph", "both_graph")
    for group, seeds, arms in [
        ("historical", replay.geom.SEEDS, (*MOVE, *locked)),
        ("fresh", FRESH, (*MOVE, *locked, "official", "double")),
    ]:
        dest = out / group
        dest.mkdir()
        replay.blind(
            work,
            pro,
            official,
            dest,
            arms,
            proteins=proteins,
            seeds=seeds,
            bond_arms=locked,
        )
        if group == "historical":
            archived = {
                (r["drug"], r["seed"], r["arm"]): r
                for r in read(HERE.parent / "results/cells.csv")
            }
            control = {"far_both": "dry", "both": "1OL7_released"}
            for row in read(dest / "cells.csv"):
                if row["arm"] in control:
                    prev = archived[row["drug"], row["seed"], control[row["arm"]]]
                    if any(row[k] != prev[k] for k in ("pool_sha256", "evals")):
                        raise ValueError(
                            "Historical production control failed; outputs retained in "
                            + str(dest)
                        )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "out", "pro", "official"):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    run(*(v.resolve() for v in (a.work, a.out, a.pro, a.official)))

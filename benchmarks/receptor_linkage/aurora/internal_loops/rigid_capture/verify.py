"""Independent crystal-frame RDKit/Gemmi and heavy-contact checks from PDB, PDBQT and raw coordinates."""

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from rdkit.Chem import rdMolAlign
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import check_predictions as audit


def rows(p):
    with p.open() as f:
        return list(csv.DictReader(f))


def check(work, result, raw):
    crystal, anchors = audit.refmol(work, "at9283")
    query, order = audit.frame.read_query(work, "at9283")
    protein = []
    for line in (
        (work / "template1OL7/trans_model/released/receptor.pdbqt")
        .read_text()
        .splitlines()
    ):
        if line.startswith(("HETATM", "ATOM  ")) and not line[77:].strip().startswith(
            "H"
        ):
            protein.append(tuple(float(line[i : i + 8]) for i in (30, 38, 46)))
    tree = cKDTree(protein)
    qraw = np.load(result / "quench.npz")
    serials = qraw["serials"].tolist()
    remap = [serials.index(s) for s in order]

    def rms(x):
        for i, y in enumerate(x):
            query.GetConformer().SetAtomPosition(i, tuple(y))
        return rdMolAlign.CalcRMS(
            query, crystal, symmetrizeConjugatedTerminalGroups=False
        )

    deltas = []
    costs = rows(result / "costs.csv")
    assert all(int(c["main_evals"]) > 0 for c in costs)
    inputs = rows(result / "quench.csv")
    for i, item in enumerate(inputs):
        for key in ("start", "return"):
            x = qraw[key + "_xyz"][i]
            deltas += [
                abs(rms(x[remap]) - float(item[key + "_A"])),
                abs(
                    float(tree.query(x)[0].min())
                    - float(item[key + "_nearest_heavy_A"])
                ),
            ]
            assert sum(map(len, tree.query_ball_point(x, 2.0))) == int(
                item[key + "_heavy_pairs_lt_2A"]
            )
    # Check every isolated rigid intervention starts from the same shape and
    # keeps the original torsions byte-for-byte, including the historical witness.
    rigid_discrepancies = []
    for source in ("complete", "inferred"):
        original = np.load(
            HERE.parent / "proposal_crossover/results" / (source + "_on_complete.npz")
        )
        lookup = {
            (int(t), int(k)): i
            for i, (t, k) in enumerate(
                zip(original["trajectory"], original["step"], strict=True)
            )
        }
        for i, item in enumerate(inputs):
            if item["source"] != source:
                continue
            j = lookup[int(item["trajectory"]), int(item["step"])]
            np.testing.assert_array_equal(
                qraw["generalized"][i, 7:], original["generalized"][j, 7:]
            )
            x, y = original["start_xyz"][j], qraw["start_xyz"][i]
            distance = lambda v: np.linalg.norm(v[:, None] - v[None, :], axis=2)
            rigid_discrepancies.append(np.max(np.abs(distance(x) - distance(y))))
    pub = {
        (int(r["trajectory"]), int(r["step"]), int(r["rank"])): r
        for r in rows(result / "modes.csv")
    }
    stage = {
        (
            int(r["trajectory"]),
            int(r["step"]),
            r["stage"],
            int(r["scope"]),
            int(r["index"]),
        ): float(r["rmsd_A"])
        for r in rows(result / "stages.csv")
    }
    branch_rows = rows(result / "branches.csv")
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in query.GetBonds()]
    template = {}
    for line in (work / "prepared/at9283.pdbqt").read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            template[int(line[6:11])] = [float(line[k : k + 8]) for k in (30, 38, 46)]

    def lengths(x):
        return np.array([np.linalg.norm(x[i] - x[j]) for i, j in bonds])

    base_edges = lengths(np.array([template[s] for s in order]))
    bond_errors = []
    for a in branch_rows:
        tr, k = int(a["trajectory"]), int(a["step"])
        path = raw / f"fit_t{tr}_s{k}.pdbqt"
        assert hashlib.sha256(path.read_bytes()).hexdigest() == a["hash"]
        values = []
        for rank, mode in enumerate(audit.frame.positions(path), 1):
            x = np.array([mode[s] for s in order])
            bond_errors.append(np.max(np.abs(lengths(x) - base_edges)))
            values.append(rms(x))
            deltas.append(abs(values[-1] - float(pub.pop((tr, k, rank))["rmsd_A"])))
        for field, value in (
            ("top1_A", values[0]),
            ("best5_A", min(values[:5])),
            ("oracle_A", min(values)),
            ("regret_A", values[0] - min(values)),
        ):
            deltas.append(abs(float(a[field]) - value))
        observations = {}
        with path.with_suffix(".xyz").open() as f:
            assert list(map(int, next(f).split()[1:])) == serials
            for line in f:
                name, owner, index, *coords = line.split()
                value = rms(np.array(coords, float).reshape(-1, 3)[remap])
                deltas.append(
                    abs(value - stage.pop((tr, k, name, int(owner), int(index))))
                )
                if name != "end_trajectory" or int(owner) == tr:
                    observations[name] = min(observations.get(name, np.inf), value)
        for name, value in observations.items():
            deltas.append(abs(value - float(a[name + "_min_A"])))
    assert (
        not stage and not pub and max(deltas) < 1e-9 and max(rigid_discrepancies) < 1e-9
    )
    return dict(
        receptor_CA_anchors=anchors,
        receptor_heavy_atoms=len(protein),
        graph_edges=len(bonds),
        modes_checked=len(bond_errors),
        branches=len(branch_rows),
        stages_checked=len(rows(result / "stages.csv")),
        quench_checks=2 * len(inputs),
        common_shapes=len(inputs),
        max_shape_pair_distance_drift_A=float(max(rigid_discrepancies)),
        max_abs_metric_error_A=float(max(deltas)),
        max_published_bond_length_change_from_prepared_A=float(max(bond_errors)),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    for k in ("work", "result", "raw"):
        cli.add_argument(k, type=Path)
    print(json.dumps(check(**vars(cli.parse_args())), indent=2))

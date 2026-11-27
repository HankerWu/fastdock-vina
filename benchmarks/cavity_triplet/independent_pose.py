"""Rebuild cross-crystal frames and audit search-stage coordinates without the placement code."""

import argparse
import csv
import gzip
import json
from pathlib import Path
import sys
import tempfile

import numpy as np
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(
    0, str(ROOT / "benchmarks/receptor_linkage/aurora/internal_loops/mutation_controls")
)
import run as independently


def records(path):
    with gzip.open(path, "rt") if path.suffix == ".gz" else path.open() as f:
        return list(csv.DictReader(f))


def validate(work, output):
    with tempfile.TemporaryDirectory() as tmp:
        data = independently.inputs(work / "aurora", work / "brd", Path(tmp))
    serials, references, mapping_counts = {}, {}, {}
    for drug, (_, lig, _, _, _, mol, serial, ref) in data.items():
        smi = next(
            l.split()[2]
            for l in lig.read_text().splitlines()
            if l.startswith("REMARK SMILES ")
        )
        typed = Chem.MolFromSmiles(smi)
        maps = ref.GetSubstructMatches(mol, uniquify=False)
        autos = typed.GetSubstructMatches(typed, uniquify=False, useChirality=True)
        # For these exact molecules, resolving all bond orders and stereo centres
        # excludes none of the element/connectivity mappings onto the deposition.
        assert {tuple(maps[0][i] for i in p) for p in autos} == set(maps)
        serials[drug] = serial
        references[drug] = ref.GetConformer().GetPositions()[np.asarray(maps)]
        mapping_counts[drug] = len(maps)

    def rms(x, refs):
        return np.sqrt(
            np.min(
                np.mean(np.sum((x[..., None, :, :] - refs) ** 2, axis=-1), axis=-1),
                axis=-1,
            )
        )

    modes = records(work / "modes.csv")
    maximum = {
        "published_A": 0.0,
        "stage_A": 0.0,
        "candidates_A": 0.0,
        "pool_endpoints_A": 0.0,
    }
    costs = {}
    nstage = 0
    for row in records(work / "cells.csv"):
        stem = row["stem"]
        pool = work / "raw" / (stem + ".pdbqt")
        drug = row["drug"]
        _, lig, _, _, _, mol, serial, ref = data[drug]
        measured = independently.measure(pool, lig, mol, serial, ref)
        expected = [m for m in modes if m["stem"] == stem]
        assert len(measured) == int(row["modes"]) == len(expected)
        errors = []
        for (err, score, drift), r in zip(measured, expected, strict=True):
            maximum["published_A"] = max(
                maximum["published_A"], abs(err - float(r["rmsd_A"]))
            )
            assert score == float(r["score"]) and drift < 0.004
            errors.append(err)
        if errors:
            stats = (errors[0], min(errors[:5]), min(errors), errors[0] - min(errors))
            maximum["pool_endpoints_A"] = max(
                maximum["pool_endpoints_A"],
                *[
                    abs(v - float(row[k]))
                    for v, k in zip(
                        stats,
                        ("top1_A", "best5_A", "oracle_A", "regret_A"),
                        strict=True,
                    )
                ],
            )
        costs[stem] = drug
    stage_keys = {
        (r["stem"], r["stage"], int(r["scope"]), int(r["index"])): float(r["rmsd_A"])
        for r in records(work / "stages.csv")
    }
    for path in sorted((work / "raw").glob("*.stages")):
        stem = path.stem
        if stem not in costs:
            continue
        drug = costs[stem]
        lines = path.read_text().splitlines()
        order = [list(map(int, lines[0][1:].split())).index(s) for s in serials[drug]]
        for l in lines[1:]:
            state, scope, index, *xyz = l.split()
            lab = stage_keys.pop((stem, state, int(scope), int(index)), None)
            if lab is None:
                continue
            x = np.array(xyz, float).reshape(-1, 3)[order]
            maximum["stage_A"] = max(
                maximum["stage_A"], abs(rms(x, references[drug]) - lab)
            )
            nstage += 1
    assert not stage_keys
    crows = records(work / "candidate_analysis/candidate_metrics.csv.gz")
    from collections import defaultdict

    grouped = defaultdict(list)
    for row in crows:
        grouped[row["source"], row["arm"]].append(row)
    for (name, arm), items in grouped.items():
        drug = name.split("_")[0]
        with (
            np.load(work / "sources" / (name + ".npz")) as s,
            np.load(work / "candidates" / (name + ".npz")) as p,
        ):
            xyz = s["xyz"]
            index = [s["serials"].tolist().index(sid) for sid in serials[drug]]
            coords = (xyz - xyz.mean(0)) @ p[arm + "_rotation"].transpose(0, 2, 1) + p[
                "center"
            ][:, None, :]
            values = rms(coords[:, index], references[drug])
            differences = [
                abs(values[int(row["index"])] - float(row["rmsd_A"])) for row in items
            ]
            maximum["candidates_A"] = max(maximum["candidates_A"], max(differences))
    assert max(maximum.values()) < 1e-8
    output.write_text(
        json.dumps(
            dict(
                pools=len(costs),
                published=len(modes),
                stage_states=nstage,
                candidate_states=len(crows),
                chemically_equivalent_mappings=mapping_counts,
                maximum_absolute_difference=maximum,
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("output", type=Path)
    a = p.parse_args()
    validate(a.work.resolve(), a.output)

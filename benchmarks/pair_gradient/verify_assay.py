"""Independently reconstruct assay decisions from archived PDBQT text and published Kd."""

import csv
import hashlib
import re
import sys
import tarfile
import tempfile
from decimal import Decimal
from itertools import combinations
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "galectin_affinity"))
from sugar import crystal_core, ligand_core, xyz

ARMS = ("upstream", "parent", "linear", "cubic", "shift")
RECEIVERS = {"Gal1": ("1GZW", "4Q26"), "Gal3": ("1KJL", "3T1L")}


def read_csv(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def main():
    labels = {
        int(r["compound"]): r
        for r in read_csv(HERE.parent / "galectin_affinity/table1.csv")
    }
    cells = read_csv(HERE / "galectins/cells.csv")
    summaries = read_csv(HERE / "galectins/summary.csv")
    indexed = {
        (r["receptor"], int(r["compound"]), int(r["seed"]), r["arm"]): r for r in cells
    }
    assert len(indexed) == len(cells) == 520
    observed = {}
    hashes = {}
    # A post-hoc partial-core diagnostic. The six experimental ring coordinates
    # never enter the predictive search or the affinity score.
    with tempfile.TemporaryDirectory() as tmp:
        fixtures = Path(tmp)
        with tarfile.open(HERE.parent / "galectin_affinity/evidence.tar.xz") as src:
            files = [
                m
                for m in src
                if (m.name.startswith("raw/") and m.name.endswith(".pdb"))
                or (m.name.startswith("prepared/") and m.name.endswith(".pdbqt"))
            ]
            src.extractall(fixtures, members=files, filter="data")
        crystal = {
            code: crystal_core(fixtures, code)
            for states in RECEIVERS.values()
            for code in states
        }
        rings = {i: ligand_core(fixtures, i) for i in labels}
    geometry = []
    with tarfile.open(HERE / "galectins/raw.tar.xz") as archive:
        for member in archive:
            if not member.name.endswith(".pdbqt"):
                continue
            receptor, compound, seed, arm = Path(member.name).stem.split("-")
            key = receptor, int(compound), int(seed), arm
            data = archive.extractfile(member).read()
            first = Decimal(
                re.search(rb"REMARK VINA RESULT:\s+(-?[0-9.]+)", data)[1].decode()
            )
            assert first == Decimal(indexed[key]["top_score"])
            assert int(indexed[key]["modes"]) == data.count(b"ENDMDL")
            # Work in exact integer millienergy units to preserve ties when pooling.
            observed[key] = int(first * 1000)
            hashes[key] = hashlib.sha256(data).hexdigest()
            distances = []
            for block in data.decode().split("ENDMDL")[:-1]:
                atoms = {
                    int(l[6:11]): xyz(l)
                    for l in block.splitlines()
                    if l.startswith(("ATOM  ", "HETATM"))
                }
                coords = np.array([atoms[i] for i in rings[int(compound)]])
                distances.append(
                    np.sqrt(((coords - crystal[receptor]) ** 2).sum(axis=1).mean())
                )
            geometry.append(
                {
                    "receptor": receptor,
                    "compound": compound,
                    "seed": seed,
                    "arm": arm,
                    "top_core_A": distances[0],
                    "best_core_A": min(distances),
                }
            )
    with (HERE / "galectins/core_geometry.csv").open("w", newline="") as file:
        w = csv.DictWriter(file, list(geometry[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(geometry)
    assert set(observed) == set(indexed)
    output = []
    for target, states in RECEIVERS.items():
        kd = {i: Decimal(row[target + "_Kd_uM"]) for i, row in labels.items()}
        pairs = [
            (strong, weak)
            for strong, weak in combinations(sorted(kd, key=kd.get), 2)
            if kd[weak] >= 3 * kd[strong]
        ]
        credit = {}
        for arm in ARMS:
            total = {
                i: sum(observed[s, i, k, arm] for k in (75467, 98251) for s in states)
                for i in labels
            }
            credit[arm] = np.array(
                [
                    int(total[a] < total[b]) + 0.5 * int(total[a] == total[b])
                    for a, b in pairs
                ]
            )
            ledger = next(
                r
                for r in summaries
                if r["target"] == target
                and r["arm"] == arm
                and r["subset"] == "all"
                and r["state"] == "consensus"
            )
            assert len(pairs) == int(ledger["pairs"])
            assert sum(credit[arm]) == float(ledger["correct"])
        rng = np.random.default_rng(20261127)
        counts = rng.multinomial(
            len(labels), [1 / len(labels)] * len(labels), size=10000
        )
        weights = np.column_stack(
            [counts[:, a - 5] * counts[:, b - 5] for a, b in pairs]
        )
        denom = weights.sum(axis=1)
        for arm in ARMS:
            delta = credit[arm] - credit["parent"]
            reps = (weights[denom > 0] @ delta) / denom[denom > 0]
            bounds = np.quantile(reps, [0.025, 0.975])
            output.append(
                {
                    "target": target,
                    "arm": arm,
                    "unique_ligands": len(labels),
                    "eligible_pairs": len(pairs),
                    "correct": float(sum(credit[arm])),
                    "delta_pp": round(100 * delta.mean(), 3),
                    "compound_boot_low_pp": round(100 * bounds[0], 3),
                    "compound_boot_high_pp": round(100 * bounds[1], 3),
                    "usable_draws": len(reps),
                }
            )
    with (HERE / "galectins/compound_sensitivity.csv").open("w", newline="") as file:
        w = csv.DictWriter(file, list(output[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(output)
    print("520 output pools verified, all ten consensus pair totals recovered.")
    keys = [k[:3] for k in observed if k[-1] == "parent"]
    print(
        "Identical upstream/parent serialized pools:",
        sum(hashes[*k, "parent"] == hashes[*k, "upstream"] for k in keys),
        "/",
        len(keys),
    )
    for k in keys:
        for field in ("hunt_evals", "second_evals", "hunt_calls", "second_calls"):
            assert indexed[*k, "shift"][field] == indexed[*k, "linear"][field]


if __name__ == "__main__":
    main()

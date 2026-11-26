"""Rebuild handed assignments from donor crystal atom names, and check selected poses and Ki order."""

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "benchmarks/cutoff_switch/hcaii"), str(ROOT / "tools")]
from check_ammonium_states import ccd_vectors, handed, read_points
from check_ammonium_states import verify as verify_rmsd
from select_zinc_pose import read_models


def table(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def audit(data, graphs, folder, relax):
    meta = json.loads((data / "prepared/metadata.json").read_text())
    catalog = json.loads((graphs / "graphs.json").read_text())
    for row in table(folder / "filter.csv"):
        rec, lig, donor, state = (
            row[k] for k in ("receptor", "ligand", "donor", "state")
        )
        sig, _ = ccd_vectors(data, lig, meta[lig])
        path = (
            data / "prepared" / f"{lig}_ligand.pdbqt"
            if state == "original"
            else folder / "inputs" / f"{lig}_reflected.pdbqt"
        )
        points = read_points(path.read_text())
        qsign = handed([points[s] for s in sig])
        qi = [meta[lig]["qt_serial"].index(s) for s in sig]
        d_sig, _ = ccd_vectors(data, donor, meta[donor])
        di = [meta[donor]["qt_serial"].index(s) for s in d_sig]
        deposited = {
            x[12:16].strip(): np.array([float(x[k : k + 8]) for k in (30, 38, 46)])
            for x in (data / "raw" / f"{donor}.pdb").read_text().splitlines()
            if x.startswith("HETATM")
            and x[21] == "A"
            and int(x[22:26]) == 303
            and x[17:20] == meta[donor]["ccd"]
            and x[16] in (" ", "A")
        }
        source = (graphs / f"{rec}_{lig}.multi").read_text().strip().split("\nOR\n")
        decisions = []
        wanted = []
        for items, block in zip(catalog[lig]["maps"], source, strict=True):
            m = dict(items)
            if not all(j in m for j in qi):
                status = "partial"
            else:
                dest = [m[j] for j in qi]
                if dest[0] != di[0] or set(dest) != set(di):
                    raise ValueError("Map crossed distinct donor chemical environments")
                points = [deposited[meta[donor]["names"][d]] for d in dest]
                status = "compatible" if qsign == handed(points) else "opposite"
            decisions.append(status)
            if status != "opposite":
                wanted.append(block)
        if decisions != row["decisions"].split():
            raise ValueError("Sign audit differs from generation")
        if not wanted:
            wanted = source
            if relax:
                wanted = list(
                    dict.fromkeys(
                        "\n".join(
                            l for l in b.splitlines() if int(l.split()[0]) not in sig
                        )
                        for b in source
                    )
                )
        actual = (
            (folder / "inputs" / f"{rec}_{lig}_{state}.anchor")
            .read_text()
            .strip()
            .split("\nOR\n")
        )
        if actual != wanted:
            raise ValueError("State anchor contains an incorrect point")
    # Independent rank calculations on every emitted candidate. Work in log-space.
    cand = {
        (r["receptor"], r["ligand"], r["seed"], r["slot"], int(r["model"])): r
        for r in table(folder / "candidates.csv")
    }
    by_cell = defaultdict(list)
    for output in sorted(folder.glob("*.pdbqt")):
        rec, lig, seed, slot = output.stem.split("_", 3)
        state = "reflected" if slot == "second_reflected" else "original"
        blocks = (
            (folder / "inputs" / f"{rec}_{lig}_{state}.anchor")
            .read_text()
            .strip()
            .split("\nOR\n")
        )
        arrays = [
            np.array([[float(z) for z in l.split()] for l in block.splitlines()])
            for block in blocks
        ]
        for n, lines in read_models(output.read_text().splitlines()).items():
            x = read_points("\n".join(lines))
            energies = []
            for arr in arrays:
                offsets = np.stack([x[int(s)] for s in arr[:, 0]]) - arr[:, 1:]
                energies.append(
                    np.mean(4 * (np.sqrt(1 + np.sum(offsets * offsets, axis=1)) - 1))
                )
            prior = -0.6 * (
                np.logaddexp.reduce(-np.array(energies) / 0.6) - np.log(len(energies))
            )
            row = cand[rec, lig, seed, slot, n]
            for arm in (
                ("original_only", "both_states")
                if slot == "common"
                else ("original_only",)
                if state == "original"
                else ("both_states",)
            ):
                by_cell[rec, lig, seed, arm].append(
                    (float(row["total"]), prior, slot != "common", n, row)
                )
    scores = {}
    for row in table(folder / "rank_cells.csv"):
        key = tuple(row[k] for k in ("receptor", "ligand", "seed", "arm"))
        pool = by_cell[key]
        method = row["method"]

        def score(a, method=method):
            return (
                a[0]
                if method == "release"
                else a[1]
                if method == "geometry_only"
                else a[0] + a[1],
                a[2],
                a[3],
            )

        top = min(pool, key=score)
        raw = min(pool, key=lambda a: (a[0], a[2], a[3]))
        expected = (top[0] - float(raw[4]["intra"])) * float(top[4]["flex_factor"])
        if (top[4]["slot"], top[3]) != (row["slot"], int(row["model"])) or abs(
            expected - float(row["native_affinity"])
        ) > 2e-7:
            raise ValueError("Independent rank or native affinity differed")
        scores[(*key, method)] = expected
    pairs = table(folder / "rank_pairs.csv")
    pair_count = 0
    for r in pairs:
        strong, weak = r["stronger"], r["weaker"]
        if meta[weak]["Ki_nM"] < 3 * meta[strong]["Ki_nM"]:
            raise ValueError("Wrong Ki pair")
        a, b = (
            scores[r["receptor"], lig, r["seed"], r["arm"], r["method"]]
            for lig in (strong, weak)
        )
        if (a < b) + 0.5 * (a == b) != float(r["correct"]):
            raise ValueError("Ki credit disagreement")
        pair_count += 1
    if pair_count != 168:
        raise ValueError(
            "Expected seven distinct pairs in each of 24 arm/condition/method combinations"
        )
    verify_rmsd(data, folder)
    print(
        f"20 donor/state assignments, {len(cand)} candidates, {len(scores)} pose selections and {pair_count} Ki pair rows independently checked"
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("data", "graphs", "folder"):
        p.add_argument(name, type=Path)
    p.add_argument("--release-conflicting-neighborhood", action="store_true")
    a = p.parse_args()
    audit(a.data, a.graphs, a.folder, a.release_conflicting_neighborhood)


if __name__ == "__main__":
    main()

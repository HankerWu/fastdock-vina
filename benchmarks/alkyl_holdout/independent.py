"""Independently verify raw deposited structures against every emitted heavy-atom pose."""

import argparse
import csv
import json
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.spatial.transform import Rotation


def xyz(line):
    return np.fromstring(line[30:54], sep=" ")


def read_frame(lines):
    return {
        int(l[22:26]): xyz(l)
        for l in lines
        if l.startswith("ATOM  ")
        and l[12:16] == " CA "
        and l[16] in (" ", "A")
        and l[21] == "A"
        and 80 <= int(l[22:26]) <= 160
        and not 107 <= int(l[22:26]) <= 115
    }


def carbon_graph(points):
    # These six ligands are unbranched aryl alkanes containing only carbon heavy atoms.
    d = squareform(pdist(points))
    g = (1.15 < d) & (d < 1.75)
    if g.sum() != 2 * len(points):
        raise ValueError("Distance graph is not a monocyclic hydrocarbon")
    return g


def isomorphisms(a, b):
    if len(a) != len(b):
        raise ValueError("Different atom counts")
    degree_a, degree_b = a.sum(axis=1), b.sum(axis=1)
    candidates = [np.where(degree_b == v)[0].tolist() for v in degree_a]
    order = sorted(range(len(a)), key=lambda j: len(candidates[j]))
    mapping = {}

    def search(depth):
        if depth == len(a):
            yield [mapping[j] for j in range(len(a))]
            return
        j = order[depth]
        for k in candidates[j]:
            if k not in mapping.values() and all(
                a[j, l] == b[k, v] for l, v in mapping.items()
            ):
                mapping[j] = k
                yield from search(depth + 1)
                del mapping[j]

    return list(search(0))


def raw_reference(work, ligand):
    raw = work / "raw"
    donor = read_frame((raw / "4W57.pdb").read_text().splitlines())
    lines = (raw / f"{ligand}.pdb").read_text().splitlines()
    frame = read_frame(lines)
    if len(donor) != 72 or set(frame) != set(donor):
        raise ValueError("Bad alignment frame")
    fixed = np.stack([donor[k] for k in sorted(frame)])
    moving = np.stack([frame[k] for k in sorted(frame)])
    rotation, _ = Rotation.align_vectors(
        fixed - fixed.mean(axis=0), moving - moving.mean(axis=0)
    )
    lig = [
        l
        for l in lines
        if l.startswith("HETATM") and l[21:26] == "A 200" and l[76:78].strip() == "C"
    ]
    alts = sorted({l[16] for l in lig} - {" "}) or [" "]
    original = {
        a.strip() or ".": np.stack([xyz(l) for l in lig if l[16] in (" ", a)])
        for a in alts
    }
    return {
        a: rotation.apply(coords - moving.mean(axis=0)) + fixed.mean(axis=0)
        for a, coords in original.items()
    }


def check(work, out):
    records = {
        (r["receptor_state"], r["ligand"], r["seed"], r["arm"], r["rank"]): r
        for r in csv.DictReader((out / "poses.csv").open())
    }
    refs = {}
    counts = 0
    worst = 0
    max_score_error = 0
    raw_cells = {}
    for posefile in sorted((work / "runs").glob("*.pdbqt")):
        state, lig, seed, arm = posefile.stem.split("_")
        if lig not in refs:
            refs[lig] = raw_reference(work, lig)
        reference = refs[lig]
        for model in posefile.read_text().strip().split("ENDMDL")[:-1]:
            lines = model.strip().splitlines()
            rank = lines[0].split()[1]
            atoms = [l for l in lines if l.startswith("ATOM")]
            coords = np.stack([xyz(l) for l in atoms])
            graph = carbon_graph(coords)
            row = records.pop((state, lig, seed, arm, rank))
            by_alt = json.loads(row["all_alt_rmsd"])
            score = next(
                float(l.partition(":")[2].split()[0])
                for l in lines
                if l.startswith("REMARK VINA RESULT")
            )
            max_score_error = max(max_score_error, abs(score - float(row["score"])))
            if rank == "1":
                raw_cells[state, seed, arm, lig] = score
            counts += 1
            for alt, crystal in reference.items():
                maps = isomorphisms(graph, carbon_graph(crystal))
                if len(maps) != (12 if lig == "4W52" else 2):
                    raise ValueError(f"Bad graph symmetry {posefile}")
                d = min(
                    np.linalg.norm(coords - crystal[mapping]) / np.sqrt(len(crystal))
                    for mapping in maps
                )
                worst = max(worst, abs(d - by_alt[alt]))
    if records or not counts or worst > 1e-6 or max_score_error != 0:
        raise ValueError((counts, worst, max_score_error))
    # Direct Kd ratios from Table 1; independently reconstruct the qualifying pairs.
    kd = {"4W52": 175, "4W53": 102, "4W54": 68, "4W55": 18}
    predicted = {}
    for r in csv.DictReader((out / "pairs.csv").open()):
        st, seed, arm = r["receptor_state"], r["seed"], r["arm"]
        a, b = r["stronger"], r["weaker"]
        margin = raw_cells[st, seed, arm, b] - raw_cells[st, seed, arm, a]
        assert kd[b] / kd[a] >= 3 and abs(margin - float(r["score_margin"])) < 1e-7
        predicted[st, seed, arm, a, b] = margin
    eligible = [
        (a, b) for a, b in combinations(sorted(kd, key=kd.get), 2) if kd[b] / kd[a] >= 3
    ]
    if len(predicted) != 3 * 3 * 2 * len(eligible):
        raise ValueError("Incomplete matched-pair comparison")
    report = {
        "poses": counts,
        "score_error": max_score_error,
        "maximum_per_altloc_rmsd_discrepancy_angstrom": worst,
        "matched_pair_rows": len(predicted),
        "unique_pairs": len(eligible),
    }
    (out / "verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(report)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    check(a.work, a.out)

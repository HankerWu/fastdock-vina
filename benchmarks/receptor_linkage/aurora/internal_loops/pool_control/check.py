"""Independent CLI, coordinate, provenance, selection and counterfactual checks."""

import argparse
import csv
import hashlib
import json
import re
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from replay import PROTEINS


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def atoms(text):
    return [l for l in text.splitlines() if l.startswith(("ATOM", "HETATM"))]


def xyz(records, heavy=True):
    return np.array(
        [
            [float(l[i : i + 8]) for i in (30, 38, 46)]
            for l in records
            if not heavy or l.split()[-1] not in ("H", "HD")
        ]
    )


def check(work, out, pro, jobs):
    data = rows(out / "scores.csv")
    cli = out / "cli"
    cli.mkdir(exist_ok=True)
    center = json.loads((work / "prepared/2X6E.json").read_text())["center"]
    box = [
        s
        for c, v in zip("xyz", center)
        for s in ("--center_" + c, str(v), "--size_" + c, "26")
    ]
    dry_records = atoms((work / PROTEINS["dry"]).read_text())
    signatures = {l[11:26] for l in dry_records}
    trees = {}
    for a, rec in PROTEINS.items():
        added = [
            l for l in atoms((work / rec).read_text()) if l[11:26] not in signatures
        ]
        if a != "dry":
            assert len(xyz(added)) == 101
            trees[a] = cKDTree(xyz(added))
    source = {
        (r["arm"], r["drug"], r["seed"], r["rank"]): float(r["rmsd_A"])
        for f in ("base", "templates")
        for r in rows(work / f / "modes.csv")
        if r["arm"] in PROTEINS
    }
    pool_cache = {}
    for d in data:
        a, q, seed, k = d["source"], d["drug"], d["seed"], d["rank"]
        key = (a, q, seed)
        if key not in pool_cache:
            name = (
                work
                / ("base" if a == "dry" else "templates")
                / "blind"
                / f"{a}_{q}_{seed}.pdbqt"
            )
            pool_cache[key] = [
                s
                for s in re.split(r"(?m)^MODEL\s+\d+\s*$", name.read_text())
                if s.strip()
            ]
        assert abs(source[a, q, seed, k] - float(d["rmsd_A"])) < 1e-12
        lines = atoms(pool_cache[key][int(k) - 1])
        assert lines == atoms((out / "poses" / (d["label"] + ".pdbqt")).read_text())
        if d["receptor"] != "dry":
            nearest = trees[d["receptor"]].query(xyz(lines))[0].min()
            assert abs(nearest - float(d["nearest_added_A"])) < 1e-12

    def score(d):
        path = cli / (d["label"] + "_" + d["receptor"] + ".log")
        cmd = [
            pro,
            "--receptor",
            work / PROTEINS[d["receptor"]],
            "--ligand",
            work / "prepared" / (d["drug"] + ".pdbqt"),
            "--pose_coordinates",
            out / "poses" / (d["label"] + ".pdbqt"),
            "--unbound_energy",
            "0",
            "--score_only",
            "--cpu",
            "1",
            *box,
        ]
        with path.open("w") as output:
            subprocess.run(cmd, stdout=output, stderr=subprocess.STDOUT, check=True)
        text = path.read_text()
        outputs = [
            float(re.search(pattern, text, re.MULTILINE)[1])
            for pattern in (
                r"Estimated Free Energy of Binding\s*:\s*([-+.\d]+)",
                r"^    Ligand - Receptor\s*:\s*([-+.\d]+)",
                r"^    Ligand\s*:\s*([-+.\d]+)",
            )
        ]
        return max(
            abs(float(d[k]) - v)
            for k, v in zip(("zero_score", "inter", "intra"), outputs)
        )

    with ThreadPoolExecutor(max_workers=jobs) as ex:
        discrepancies = list(ex.map(score, data))
    assert max(discrepancies) <= 0.000501
    indices = defaultdict(list)
    for d in data:
        for src in (d["source"], "union"):
            indices[d["drug"], d["seed"], src, d["receptor"]].append(d)
    for r in rows(out / "pools.csv"):
        group = indices[r["drug"], r["seed"], r["source"], r["receptor"]]
        minimum = min(float(d["zero_score"]) for d in group)
        expected = {d["label"] for d in group if float(d["zero_score"]) == minimum}
        assert expected == set(r["winners"].split("|"))
        assert min(float(d["rmsd_A"]) for d in group) == float(r["oracle_A"])
    far_inputs = atoms((out / "raw_far/far.pdbqt").read_text())
    assert far_inputs[: len(dry_records)] == dry_records
    original = atoms((work / PROTEINS["1OL7_released"]).read_text())
    tail = original[len(dry_records) :]
    new = far_inputs[len(dry_records) :]
    assert [l[:30] + l[54:] for l in tail] == [l[:30] + l[54:] for l in new]
    assert np.allclose(
        xyz(new, False) - xyz(tail, False), [100, 0, 0], rtol=0, atol=1e-10
    )
    key = r"(?m)^(Search trajectory:|Hunt chain:|Pocket start:|Search counters:).*$"

    def get_lines(path):
        return [
            m.group() for m in re.finditer(key, path.with_suffix(".log").read_text())
        ]

    equal = 0
    for r in rows(out / "far.csv"):
        dry = work / "base/blind" / f"dry_{r['drug']}_{r['seed']}.pdbqt"
        displaced = out / "raw_far" / f"far_{r['drug']}_{r['seed']}.pdbqt"
        assert digest(dry) == r["dry_sha256"]
        assert digest(displaced) == r["far_sha256"]
        equal += dry.read_bytes() == displaced.read_bytes()
        assert get_lines(dry) == get_lines(displaced) and get_lines(dry)
    print(
        json.dumps(
            {
                "cli_cases": len(data),
                "max_error_to_three_decimals": max(discrepancies),
                "distance_checks": sum(d["receptor"] != "dry" for d in data),
                "original_rmsd_labels": len(source),
                "fixed_pool_winner_checks": len(indices),
                "complete_far_pools_identical": equal,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    for key in ("work", "results", "pro"):
        ap.add_argument(key, type=Path)
    ap.add_argument("--jobs", type=int, default=2)
    a = ap.parse_args()
    check(a.work.resolve(), a.results.resolve(), a.pro.resolve(), a.jobs)

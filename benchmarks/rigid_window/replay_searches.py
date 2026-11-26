"""Replay fixed-budget Late searches and match printed coordinates irrespective of rank."""

import argparse
import hashlib
import itertools
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
from scipy.optimize import linear_sum_assignment
from run import write


def models(path):
    result = []
    for block in path.read_text().split("ENDMDL")[:-1]:
        atoms = [l for l in block.splitlines() if l.startswith(("ATOM", "HETATM"))]
        xyz = np.asarray([[float(l[k : k + 8]) for k in (30, 38, 46)] for l in atoms])
        inter = float(
            next(
                x.split(":")[1]
                for x in block.splitlines()
                if x.startswith("REMARK INTER:")
            )
        )
        score = float(
            next(
                x.split(":")[1].split()[0]
                for x in block.splitlines()
                if x.startswith("REMARK VINA RESULT:")
            )
        )
        result.append((xyz, inter, score))
    return result


def jobs(base):
    for panel, data in (
        ("brd4", "fdv-brd4"),
        ("t4l", "fdv-pro-fixtures"),
        ("hcaii", "fdv-second"),
        ("alkyl", "fdv-alkyl"),
    ):
        root = (
            base / ("fdv-map-" + panel) if panel != "alkyl" else base / "fdv-alkyl/runs"
        )
        prepared = base / data / "prepared"
        for path in sorted(root.glob("*.pdbqt")):
            rec, lig, seed, arm = path.stem.split("_")
            extra = []
            if arm != "native":
                extra = [
                    "--search_anchor",
                    (
                        root / f"{rec}_{lig}.{arm}"
                        if panel != "alkyl"
                        else prepared / f"{lig}.{arm}"
                    ),
                ]
            yield (
                panel,
                path,
                [
                    "--ligand",
                    prepared / f"{lig}_ligand.pdbqt",
                    "--receptor",
                    prepared / f"{rec}_receptor.pdbqt",
                    "--config",
                    prepared / (f"{rec}_box.txt" if panel != "alkyl" else "box.txt"),
                    "--cpu",
                    "1",
                    "--seed",
                    seed,
                    "--exhaustiveness",
                    "4",
                    "--max_evals",
                    "50000",
                    "--num_modes",
                    "9",
                    "--energy_range",
                    "100",
                    "--pose_selection",
                    "late",
                    "--search_counters",
                    *extra,
                ],
            )


def run(entry, binary, out):
    panel, old, cmd = entry
    dest = out / panel / old.name
    with dest.with_suffix(".log").open("w") as log:
        subprocess.run(
            [binary, *cmd, "--out", dest],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    same = old.read_bytes() == dest.read_bytes()
    a, b = models(old), models(dest)
    distances = np.array(
        [[np.linalg.norm(x[0] - y[0], axis=1).max() for y in b] for x in a]
    )
    left, right = linear_sum_assignment(distances)
    # A sub-milliangstrom rigid repair changes each printed component by at
    # most two 0.001-A PDBQT steps. Distinct retained modes can enter/leave.
    matches = {
        i: j for i, j in zip(left, right, strict=True) if distances[i, j] < 0.004
    }
    changes = []

    def record(i, j):
        return dict(
            panel=panel,
            pool=old.stem,
            old_rank=i + 1 if i is not None else "",
            new_rank=j + 1 if j is not None else "",
            max_atom_delta=distances[i, j] if i is not None and j is not None else "",
            old_inter=a[i][1] if i is not None else "",
            new_inter=b[j][1] if j is not None else "",
            old_affinity=a[i][2] if i is not None else "",
            new_affinity=b[j][2] if j is not None else "",
        )

    if not same:
        for i in range(len(a)):
            j = matches.get(i)
            if j is None or i != j or distances[i, j] > 0 or a[i][1:] != b[j][1:]:
                changes.append(record(i, j))
        changes.extend(
            record(None, j) for j in range(len(b)) if j not in matches.values()
        )
    return dict(
        panel=panel,
        pool=old.stem,
        old_modes=len(a),
        new_modes=len(b),
        byte_identical=int(same),
        top_matched_rank=int(matches[0]) + 1 if 0 in matches else "",
        top_score_delta=b[0][2] - a[0][2],
        old_sha=hashlib.sha256(old.read_bytes()).hexdigest(),
        new_sha=hashlib.sha256(dest.read_bytes()).hexdigest(),
    ), changes


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("base", type=Path)
    p.add_argument("binary", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("--jobs", type=int, default=4)
    a = p.parse_args()
    for panel in ("brd4", "t4l", "hcaii", "alkyl"):
        (a.out / panel).mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(a.jobs) as pool:
        results = list(
            pool.map(lambda e: run(e, a.binary.resolve(), a.out), jobs(a.base))
        )
    write(a.out / "replay.csv", [r for r, _ in results])
    diffs = list(itertools.chain.from_iterable(v for _, v in results))
    if diffs:
        write(a.out / "changed_models.csv", diffs)
    print(len(results), sum(r["byte_identical"] for r, _ in results), len(diffs))

"""Redock all 282 Gal3 cells and compare released output to a printed-coordinate spacing oracle."""

import argparse
import hashlib
import math
import sys
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "accepted_retention/affinity"))
from run import COUNT, cells, dock, read, write


def blocks(data):
    return [
        s.split(b"\n", 1)[1] + b"ENDMDL\n" for s in data.split(b"ENDMDL\n") if s.strip()
    ]


def positions(model):
    return [
        (int(s[6:11]), tuple(float(s[k : k + 8]) for k in (30, 38, 46)))
        for s in model.splitlines()
        if s.startswith((b"ATOM  ", b"HETATM")) and s.split()[-1] not in {b"HD", b"H"}
    ]


def expected(data):
    candidates = blocks(data)
    keep, coords = [], []
    uncertain = False
    for i, block in enumerate(candidates):
        pos = positions(block)
        distances = []
        for old in coords:
            if [a[0] for a in pos] != [a[0] for a in old]:
                raise ValueError("Heavy atom order differs")
            distances.append(
                math.sqrt(
                    sum(math.dist(a[1], b[1]) ** 2 for a, b in zip(pos, old)) / len(pos)
                )
            )
        if any(abs(d - 1.0) <= math.sqrt(3) * 0.001 for d in distances):
            uncertain = True
        if all(d >= 1.0 for d in distances):
            keep.append(i)
            coords.append(pos)
    return [candidates[i] for i in keep], keep, uncertain


def main(args):
    dest = args.work.resolve()
    dest.mkdir(parents=True, exist_ok=True)
    rows = []
    for name, ntraj in (("single", 1), ("four", 4)):
        study = dest / name
        study.mkdir()
        with tarfile.open(HERE.parent / "affinity" / (name + ".tar.xz")) as arc:
            arc.extractall(study, filter="data")
        work = dest / (name + "-release")
        (work / "raw").mkdir(parents=True)
        configs = {(c[1], c[2]): c for c in cells(work, "Gal3")}
        for row in read(study / "runs.csv"):
            if row["arm"] != "candidate":
                continue
            cell = configs[row["receptor"], row["compound"]]
            result = dock(
                cell,
                work,
                str(args.binary.resolve()),
                "candidate",
                int(row["seed"]),
                int(row["cap"]),
                True,
                trajectories=ntraj,
            )
            key = result["key"]
            old = (study / "raw" / (key + ".pdbqt")).read_bytes()
            new = (work / "raw" / (key + ".pdbqt")).read_bytes()
            wanted, keep, ambiguous = expected(old)
            same = blocks(new) == wanted
            if (
                not same
                or blocks(new)[0] != blocks(old)[0]
                or any(int(row[k]) != int(result[k]) for k in (*COUNT, "joint"))
            ):
                raise ValueError(
                    f"Output or work not explained by final separation: {key}"
                )
            rows.append(
                {
                    "layout": name,
                    "key": key,
                    "identical": old == new,
                    "kept_original_ranks": " ".join(str(i + 1) for i in keep),
                    "before": len(blocks(old)),
                    "after": len(wanted),
                    "near_rounding_boundary": ambiguous,
                    "first_score": row["first_score"],
                    "joint": result["joint"],
                    "elapsed_ns": result["elapsed_ns"],
                    "sha256": hashlib.sha256(new).hexdigest(),
                }
            )
    write(dest / "release.csv", rows)
    print(
        len(rows),
        "checked;",
        sum(not r["identical"] for r in rows),
        "changed;",
        sum(r["near_rounding_boundary"] for r in rows),
        "rounded boundary cases",
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("binary", type=Path)
    main(p.parse_args())

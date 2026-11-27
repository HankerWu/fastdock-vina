"""Diagnose active proposals at their actual incoming states and verify first paired differences."""

import argparse
import csv
import functools
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "macrocycle_proposal" / "external"))
from external_score import poses
from external_run import write


def load(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def fitted(x, ref):
    a = x - x.mean(axis=0)
    b = ref - ref.mean(axis=0)
    u, _, v = np.linalg.svd(a.T @ b)
    u[:, -1] *= np.linalg.det(u @ v)
    return np.mean(np.sum((a @ u @ v - b) ** 2, axis=1)) ** 0.5


def main(inputs, stages, out):
    meta = {}

    @functools.lru_cache(maxsize=32)
    def coords(key, stage):
        pdb = key.split("_")[0]
        if pdb not in meta:
            meta[pdb] = json.loads((inputs / pdb / "preparation.json").read_text())
        f = stages / "raw" / (key + "." + stage + ".pdbqt")
        ids = [
            int(l.split()[1])
            for l in f.read_text().splitlines()
            if l.startswith("MODEL ")
        ]
        serial = meta[pdb]["serial"]
        return {i: np.array([p[s][1] for s in serial]) for i, p in zip(ids, poses(f))}

    def metrics(x, refs):
        d = np.sqrt(np.mean(np.sum((x - refs) ** 2, axis=2), axis=1))
        return (
            d.min(),
            min(fitted(x, r) for r in refs),
            np.linalg.norm(x.mean(axis=0) - refs[0].mean(axis=0)),
        )

    events = load(stages / "events.csv")
    trials = []
    for r in events:
        key = r["key"]
        if not key.endswith("_concerted"):
            continue
        step = int(r["step"])
        source = coords(key, "before")[step]
        refs = np.array(meta[key.split("_")[0]]["mapped"])
        row = {**r, "accepted": int(step in coords(key, "accepted"))}
        for stage in ("before", "proposal", "hunt"):
            x = coords(key, stage)[step]
            for k, v in zip(("pose_A", "internal_A", "centroid_A"), metrics(x, refs)):
                row[stage + "_" + k] = v
        trials.append(row)
    write(out / "trial_changes.csv", trials)
    first = []
    for key in sorted({t["key"] for t in trials}):
        control = key.replace("_concerted", "_null")
        new, old = coords(key, "proposal"), coords(control, "proposal")
        for step in sorted(set(new) & set(old)):
            if np.array_equal(new[step], old[step]):
                continue
            if not np.array_equal(
                coords(key, "before")[step], coords(control, "before")[step]
            ):
                raise ValueError("A preceding difference is already visible at " + key)
            er = {
                r["key"]: r
                for r in events
                if r["key"] in (key, control) and int(r["step"]) == step
            }
            if not (
                len(er) == 2
                and er[key]["draw_rad"] == er[control]["draw_rad"]
                and er[key]["axis"] == er[control]["axis"]
            ):
                raise ValueError(
                    "First difference is not a shared incoming ring draw: " + key
                )
            refs = np.array(meta[key.split("_")[0]]["mapped"])
            t = {
                "key": key,
                "step": step,
                "axis": er[key]["axis"],
                "draw_rad": er[key]["draw_rad"],
            }
            for arm, k in (("null", control), ("concerted", key)):
                for stage in ("proposal", "hunt"):
                    t[arm + "_" + stage + "_A"] = metrics(coords(k, stage)[step], refs)[
                        0
                    ]
            first.append(t)
            break
    write(out / "first_difference.csv", first)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("inputs", "stages", "out"):
        p.add_argument(k, type=Path)
    a = p.parse_args()
    main(a.inputs, a.stages, a.out)

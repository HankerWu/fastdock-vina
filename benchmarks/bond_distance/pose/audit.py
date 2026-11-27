"""Recompute phase counts and symmetry-corrected per-mode receptor-frame RMSD from raw files."""

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

import numpy as np


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def independent_rmsds(query, receiver, path):
    refs = np.asarray(query["mapped"])
    center_q, center_r = map(np.asarray, (query["center"], receiver["center"]))
    keys = sorted(
        k
        for k in query["ca"].keys() & receiver["ca"].keys()
        if np.linalg.norm(query["ca"][k] - center_q) < 14
        and np.linalg.norm(receiver["ca"][k] - center_r) < 14
    )
    p, q = [np.array([m["ca"][k] for k in keys]) for m in (query, receiver)]
    u, _, v = np.linalg.svd((p - p.mean(0)).T @ (q - q.mean(0)))
    d = np.diag([1, 1, np.sign(np.linalg.det(u @ v))])
    refs = (refs - p.mean(0)) @ (u @ d @ v) + q.mean(0)
    result = []
    for model in path.read_text().split("ENDMDL")[:-1]:
        atoms = {
            int(l[6:11]): [float(l[k : k + 8]) for k in (30, 38, 46)]
            for l in model.splitlines()
            if l.startswith(("ATOM", "HETATM"))
        }
        xyz = np.array([atoms[i] for i in query["serial"]])
        result.append(
            float(np.min(np.linalg.norm(xyz - refs, axis=(1, 2))) / np.sqrt(len(xyz)))
        )
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    w = p.parse_args().work
    meta = {
        p.parent.name: json.loads(p.read_text()) for p in w.glob("*/preparation.json")
    }
    report = {}
    for runfile, modefile in (
        ("cells", "modes"),
        ("cost", "cost_modes"),
        ("confirmation", "confirmation_modes"),
    ):
        rows = read(w / (runfile + ".csv"))
        modes = {(r["key"], int(r["mode"])): r for r in read(w / (modefile + ".csv"))}
        observed = 0
        hits = {}
        max_error = 0
        for r in rows:
            key = r["key"]
            file = w / "runs" / (key + ".pdbqt")
            if r.get("sha256"):
                assert sha(file) == r["sha256"], key
            log = (
                w / "runs" / (key + (".counts" if r["arm"] == "upstream" else ".log"))
            ).read_text()
            pairs = dict(
                re.findall(
                    r"(hunt_calls|second_calls|explicit_calls|hunt_evals|second_evals|explicit_evals)=(\d+)",
                    log,
                )
            )
            assert len(pairs) == 6 and all(r[k] == v for k, v in pairs.items()), key
            vals = independent_rmsds(meta[r["query"]], meta[r["receiver"]], file)
            assert len(vals) == int(r["modes"]), key
            for i, x in enumerate(vals, 1):
                max_error = max(max_error, abs(float(modes[key, i]["rmsd_A"]) - x))
            assert abs(vals[0] - float(r["top1_A"])) < 1e-7, key
            assert abs(min(vals[:5]) - float(r["best5_A"])) < 1e-7, key
            assert abs(min(vals) - float(r["oracle_A"])) < 1e-7, key
            observed += len(vals)
            h = hits.setdefault(
                r["arm"], {"cells": 0, "top1": 0, "top5": 0, "oracle": 0}
            )
            for k, v in (
                ("cells", 1),
                ("top1", vals[0] < 2),
                ("top5", min(vals[:5]) < 2),
                ("oracle", min(vals) < 2),
            ):
                h[k] += v
        assert max_error < 1e-7
        assert observed == len(modes)
        report[runfile] = {
            "runs": len(rows),
            "modes": observed,
            "max_absolute_RMSD_disagreement_A": max_error,
            "hits": hits,
        }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

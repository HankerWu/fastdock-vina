"""Run or independently re-read the four frozen cross-docking arms."""

import argparse
import csv
import hashlib
import itertools
import json
import re
import resource
import subprocess
import time
from pathlib import Path

import numpy as np

PDBS = ("4W9H", "5LLI")
ARMS = ("official", "parent", "fixed", "parent_double")


def write_csv(path, rows):
    with path.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def frame(query, receptor):
    keys = sorted(query["ca"].keys() & receptor["ca"].keys())
    assert len(keys) >= 8
    a, b = [np.array([x["ca"][k] for k in keys]) for x in (query, receptor)]
    x, y = a - a.mean(0), b - b.mean(0)
    u, _, v = np.linalg.svd(x.T @ y)
    rot = (u * np.array([1, 1, np.linalg.det(u @ v)])) @ v
    offset = b.mean(0) - a.mean(0) @ rot
    return (
        np.asarray(query["crystal"]) @ rot + offset,
        len(keys),
        np.sqrt(((a @ rot + offset - b) ** 2).sum(1).mean()),
    )


def run(work, bins):
    meta = {p: json.loads((work / p / "metadata.json").read_text()) for p in PDBS}
    runs = work / "runs"
    runs.mkdir(exist_ok=True)
    modes, cells, costs, binary_info = [], [], [], {}
    if bins:
        for arm, path in bins.items():
            binary_info[arm] = {
                "sha256": sha(path),
                "version": subprocess.check_output(
                    [path, "--version"], text=True
                ).strip(),
            }
        (work / "binaries.json").write_text(json.dumps(binary_info, indent=2) + "\n")
    for n, ((q, r), seed) in enumerate(
        itertools.product(itertools.permutations(PDBS), [27183, 91439])
    ):
        ref, n_ca, res = frame(meta[q], meta[r])
        for arm in ARMS[n % 4 :] + ARMS[: n % 4]:
            key = f"{q}_{r}_{seed}_{arm}"
            out = runs / (key + ".pdbqt")
            log = runs / (key + ".log")
            if bins:
                cmd = [
                    str(bins[arm]),
                    "--receptor",
                    f"{r}/receptor.pdbqt",
                    "--ligand",
                    f"{q}/ligand.pdbqt",
                    "--cpu",
                    "1",
                    "--seed",
                    str(seed),
                    "--exhaustiveness",
                    "2" if arm == "parent_double" else "1",
                    "--max_evals",
                    "12000",
                    "--num_modes",
                    "9",
                    "--energy_range",
                    "100",
                    "--out",
                    f"runs/{out.name}",
                ]
                for name, val in zip("xyz", meta[r]["center"]):
                    cmd.extend(
                        [f"--center_{name}", f"{val:.7f}", f"--size_{name}", "24"]
                    )
                if arm != "official":
                    cmd.append("--search_counters")
                before = resource.getrusage(resource.RUSAGE_CHILDREN)
                start = time.monotonic()
                with log.open("w") as trace:
                    subprocess.run(
                        cmd,
                        cwd=work,
                        stdout=trace,
                        stderr=subprocess.STDOUT,
                        check=True,
                    )
                stop = time.monotonic()
                after = resource.getrusage(resource.RUSAGE_CHILDREN)
                costs.append(
                    dict(
                        key=key,
                        wall_s=stop - start,
                        cpu_s=(after.ru_utime + after.ru_stime)
                        - (before.ru_utime + before.ru_stime),
                        arguments=json.dumps(cmd[1:]),
                    )
                )
            text = log.read_text()
            gradients = re.findall(r"hunt_evals=(\d+).*?second_evals=(\d+)", text)
            explicit = re.findall(r"explicit_evals=(\d+)", text)
            if arm != "official":
                assert gradients and explicit
            local = []
            for block in out.read_text().split("ENDMDL"):
                if not block.strip():
                    continue
                coords = {
                    int(l[6:11]): [float(l[k : k + 8]) for k in (30, 38, 46)]
                    for l in block.splitlines()
                    if l.startswith(("ATOM", "HETATM"))
                }
                xyz = np.array([coords[k] for k in meta[q]["serial"]])
                minima = np.sqrt(np.sum((ref - xyz) ** 2, axis=2).mean(axis=1))
                rank = len(local) + 1
                value = min(minima)
                energy = float(
                    re.search(r"REMARK VINA RESULT:\s+([-\d.]+)", block).group(1)
                )
                assert np.isfinite(energy) and np.isfinite(value)
                local.append(value)
                modes.append(
                    dict(
                        query=q,
                        receptor=r,
                        seed=seed,
                        arm=arm,
                        rank=rank,
                        rmsd=value,
                        score=energy,
                    )
                )
            assert 1 <= len(local) <= 9
            cells.append(
                dict(
                    query=q,
                    receptor=r,
                    seed=seed,
                    arm=arm,
                    modes=len(local),
                    aligned_ca=n_ca,
                    ca_rmsd=res,
                    top1=local[0],
                    best5=min(local[:5]),
                    oracle=min(local),
                    regret=local[0] - min(local),
                    output_sha256=sha(out),
                    search_gradient=sum(int(a) + int(b) for a, b in gradients)
                    if gradients
                    else "",
                    explicit_gradient=sum(map(int, explicit)) if explicit else "",
                )
            )
    if bins:
        write_csv(work / "cost.csv", costs)
    write_csv(work / "modes.csv", modes)
    write_csv(work / "cells.csv", cells)
    summary = []
    for arm in ARMS:
        obs = [r for r in cells if r["arm"] == arm]
        top = np.array([r["top1"] for r in obs])
        row = dict(
            arm=arm, cells=len(obs), top1_mean=top.mean(), top1_median=np.median(top)
        )
        for q in [25, 75, 90, 95]:
            row[f"top1_p{q}"] = np.percentile(top, q)
        for k in ("best5", "oracle", "regret"):
            row[f"{k}_mean"] = np.mean([r[k] for r in obs])
        summary.append(row)
    write_csv(work / "summary.csv", summary)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("work", type=Path)
    for b in ["official", "parent", "fixed"]:
        ap.add_argument("--" + b, type=Path)
    opt = ap.parse_args()
    bins = None
    if any([opt.official, opt.parent, opt.fixed]):
        if not all([opt.official, opt.parent, opt.fixed]):
            ap.error("supply all three binaries to run docking")
        bins = {x: getattr(opt, x).resolve() for x in ["official", "parent", "fixed"]}
        bins["parent_double"] = bins["parent"]
    run(opt.work.resolve(), bins)

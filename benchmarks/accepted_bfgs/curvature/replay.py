"""Replay both native builds on every condition in the archived glucose experiment."""

import argparse
import hashlib
import itertools
import json
import re
import shlex
import sys
import tarfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[1]
GLUCOSE = BENCH / "ligand_hydroxyl/glucose"
sys.path[:0] = [str(GLUCOSE), str(BENCH / "hydration_transfer")]
from independent import independent_mapping, xyz
from panel import timed_process, write_csv

CAPS = (3000, 12000)
SEEDS = (43103, 43117)


def prepare(work):
    archive = GLUCOSE / "evidence.tar.xz"
    with tarfile.open(archive) as tar:
        files = [
            m
            for m in tar.getmembers()
            if m.isfile()
            and (not m.name.startswith("runs/") or m.name.endswith("_upstream.pdbqt"))
        ]
        tar.extractall(work, members=files, filter="data")
    return json.loads((work / "metadata.json").read_text())


def dock_args(case, cap, seed, entry, out):
    pdb, compound = case.split("-")
    args = [
        "--receptor",
        f"{pdb}/receptor.pdbqt",
        "--ligand",
        f"{pdb}/{compound}.pdbqt",
        "--cpu",
        "1",
        "--seed",
        str(seed),
        "--exhaustiveness",
        "4",
        "--max_evals",
        str(cap),
        "--num_modes",
        "9",
        "--energy_range",
        "100",
        "--out",
        str(out),
        "--search_counters",
    ]
    for axis, value in zip("xyz", entry["center"]):
        args += [f"--center_{axis}", str(value), f"--size_{axis}", "24"]
    return args


def evaluate(pool, ids, ref):
    rows = []
    for block in pool.read_text().split("ENDMDL"):
        if not block.strip():
            continue
        lines = block.splitlines()
        coords = {
            int(l[6:11]): xyz(l) for l in lines if l.startswith(("ATOM", "HETATM"))
        }
        d = np.array([coords[i] for i in ids]) - ref
        score = float(
            next(l.split()[3] for l in lines if l.startswith("REMARK VINA RESULT:"))
        )
        rows.append(
            dict(
                rank=len(rows) + 1,
                score=score,
                rmsd_A=float(np.sqrt((d * d).sum() / len(ids))),
            )
        )
    if not rows:
        raise ValueError(f"No output models in {pool}")
    return rows


def run(parent, fixed, work):
    work.mkdir(parents=True, exist_ok=False)
    meta = prepare(work)
    poses, cells, costs = [], [], []
    for case, m in sorted(meta.items()):
        ids, ref = independent_mapping(work, *case.split("-"))
        for j, (cap, seed) in enumerate(itertools.product(CAPS, SEEDS)):
            arms = [("parent", parent), ("fixed", fixed)]
            if j % 2:
                arms.reverse()
            for arm, binary in [("upstream", None), *arms]:
                tag = f"{case}_{cap}_{seed}_{arm}"
                keys = dict(case=case, cap=cap, seed=seed, arm=arm)
                out = Path("runs") / (tag + ".pdbqt")
                if binary:
                    log = work / "runs" / (tag + ".log")
                    args = dock_args(case, cap, seed, m, out)
                    cost = timed_process([str(binary), *args], work, log)
                    txt = log.read_text()
                    pairs = re.findall(r"hunt_evals=(\d+).*?second_evals=(\d+)", txt)
                    ex = re.findall(r"explicit_evals=(\d+)", txt)
                    if len(pairs) != 4 or len(ex) != 1:
                        raise ValueError(f"Incomplete search counters in {log}")
                    costs.append(
                        {
                            **keys,
                            **cost,
                            "search_FG": sum(int(a) + int(b) for a, b in pairs),
                            "explicit_FG": int(ex[0]),
                            "args": shlex.join(args),
                        }
                    )
                rows = evaluate(work / out, ids, ref)
                poses.extend({**keys, **r} for r in rows)
                rmsd = [r["rmsd_A"] for r in rows]
                cells.append(
                    dict(
                        **keys,
                        modes=len(rows),
                        top1_A=rmsd[0],
                        best5_A=min(rmsd[:5]),
                        oracle_A=min(rmsd),
                        regret_A=rmsd[0] - min(rmsd),
                        pool_sha256=hashlib.sha256(
                            (work / out).read_bytes()
                        ).hexdigest(),
                    )
                )
    for name, rows in [("cells", cells), ("poses", poses), ("costs", costs)]:
        write_csv(work / f"{name}.csv", rows)
    summary = {}
    for cap in CAPS:
        for arm in ["upstream", "parent", "fixed"]:
            rows = [r for r in cells if r["cap"] == cap and r["arm"] == arm]
            values = [r["top1_A"] for r in rows]
            timing = [r for r in costs if r["cap"] == cap and r["arm"] == arm]
            summary[f"{cap}/{arm}"] = dict(
                n=len(rows),
                mean=np.mean(values),
                percentiles=dict(
                    zip(
                        ("q25", "median", "q75", "p90", "p95"),
                        np.percentile(values, [25, 50, 75, 90, 95]),
                    )
                ),
                mean_best5=np.mean([r["best5_A"] for r in rows]),
                mean_oracle=np.mean([r["oracle_A"] for r in rows]),
                mean_regret=np.mean([r["regret_A"] for r in rows]),
                top1_gt5=sum(v > 5 for v in values),
                total_cost=(
                    {
                        k: sum(r[k] for r in timing)
                        for k in ("search_FG", "explicit_FG", "wall_s", "cpu_s")
                    }
                    if timing
                    else None
                ),
            )
    (work / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (work / "binaries.json").write_text(
        json.dumps(
            {
                a: hashlib.sha256(p.read_bytes()).hexdigest()
                for a, p in [("parent", parent), ("fixed", fixed)]
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    for key in ("parent", "fixed", "work"):
        cli.add_argument(key, type=Path)
    cfg = cli.parse_args()
    run(cfg.parent.resolve(), cfg.fixed.resolve(), cfg.work.resolve())

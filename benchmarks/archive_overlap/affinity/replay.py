"""Run or independently audit the matched, known-site Gal3 affinity panel."""

import argparse
import hashlib
import itertools
import json
import os
import random
import re
import resource
import subprocess
import sys
from decimal import Decimal
from pathlib import Path

from scipy.stats import kendalltau, spearmanr

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[1]
sys.path.insert(0, str(BENCH / "accepted_retention/affinity"))
from run import COUNT, cells, dock, read, write

BUDGETS = (3000, 12000)
SEED = 391781
ARMS = ("upstream", "parent", "candidate")
RECEPTORS = ("1KJL", "3T1L")


def digest(data):
    return hashlib.sha256(data).hexdigest()


def cpu():
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    return usage.ru_utime + usage.ru_stime


def budgets(trajectories):
    return BUDGETS if trajectories == 1 else (3000,)


def execute(args):
    work = args.work
    (work / "raw").mkdir(parents=True, exist_ok=True)
    os.sched_setaffinity(0, {args.cpu})
    entries = list(itertools.product(budgets(args.trajectories), cells(work, "Gal3")))
    random.Random(617).shuffle(entries)
    rows = []
    for index, (budget, cell) in enumerate(entries):
        shift = index % len(ARMS)
        for arm in ARMS[shift:] + ARMS[:shift]:
            before = cpu()
            row = dock(
                cell,
                work,
                getattr(args, arm),
                arm,
                SEED,
                budget,
                arm != "upstream",
                trajectories=args.trajectories,
            )
            row["cpu_s"] = round(cpu() - before, 6)
            if arm == "upstream":
                counted = dock(
                    cell,
                    work,
                    args.upstream_counts,
                    arm,
                    SEED,
                    budget,
                    False,
                    diag=True,
                    trajectories=args.trajectories,
                )
                prefix = work / "raw" / row["key"]
                if (
                    prefix.with_suffix(".replay").read_bytes()
                    != prefix.with_suffix(".pdbqt").read_bytes()
                ):
                    raise ValueError(row["key"] + " upstream counting replay differs")
                prefix.with_suffix(".replay").unlink()
                row.update({k: counted[k] for k in ("joint", *COUNT)})
            rows.append(row)
    write(work / "runs.csv", rows)
    meta = {
        "seed": SEED,
        "cpu": args.cpu,
        "trajectories": args.trajectories,
        "max_child_rss_KiB": resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
        "sha256": {
            name: digest(getattr(args, name).read_bytes())
            for name in (*ARMS, "upstream_counts")
        },
    }
    (work / "engines.json").write_text(json.dumps(meta, indent=2) + "\n")


def fixed(args):
    entries = {(c[1], c[2]): c for c in cells(args.work, "Gal3")}
    pose = args.work / "fixed.pdbqt"
    fixed_rows = []
    for row in read(args.work / "runs.csv"):
        if row["arm"] != "parent":
            continue
        _, _, _, protein, box, ligand = entries[row["receptor"], row["compound"]]
        prefix = args.work / "raw" / row["key"]
        pose.write_bytes(
            prefix.with_suffix(".pdbqt").read_bytes().split(b"ENDMDL", 1)[0]
            + b"ENDMDL\n"
        )
        check_box = box
        for attempt in range(2):
            if attempt:
                check_box = args.work / "fixed.box"
                check_box.write_text(
                    re.sub(
                        r"(size_[xyz]\s*=\s*)([.\d]+)",
                        lambda m: m[1] + str(float(m[2]) + 0.75),
                        box.read_text(),
                    )
                )
            terms = []
            for arm in ("parent", "candidate"):
                out = subprocess.run(
                    [
                        getattr(args, arm),
                        "--receptor",
                        protein,
                        "--ligand",
                        ligand,
                        "--pose_coordinates",
                        pose,
                        "--config",
                        check_box,
                        "--score_only",
                        "--cpu",
                        "1",
                        "--seed",
                        str(SEED),
                    ],
                    capture_output=True,
                    check=False,
                    text=True,
                )
                if out.returncode and "outside the grid box" not in out.stderr:
                    raise ValueError(out.stderr)
                terms.append(re.findall(r":\s*([-+.\deE]+)\s+\(kcal/mol\)", out.stdout))
            if terms[0] != terms[1]:
                raise ValueError(row["key"] + " scoring mismatch")
            if len(terms[0]) == 10:
                break
            if terms[0] or attempt:
                raise ValueError(row["key"] + " incomplete scoring terms")
        fixed_rows.append(
            {
                "key": row["key"],
                "pose_sha256": digest(pose.read_bytes()),
                "expanded_box": bool(attempt),
                "terms": " ".join(terms[0]),
            }
        )
    write(args.work / "fixed.csv", fixed_rows)
    pose.unlink()
    (args.work / "fixed.box").unlink(missing_ok=True)


def audit(work):
    manifest = read(BENCH / "pair_gradient/galectin_external/manifest.csv")
    labels = {r["compound"]: r for r in manifest}
    kd = {i: Decimal(r["kd_uM"]) for i, r in labels.items()}
    compounds = sorted(labels, key=kd.get)
    runs = read(work / "runs.csv")
    trajectories = json.loads((work / "engines.json").read_text()).get(
        "trajectories", 1
    )
    caps = budgets(trajectories)
    records = {}
    for row in runs:
        prefix = work / "raw" / row["key"]
        raw = prefix.with_suffix(".pdbqt").read_bytes()
        scores = re.findall(rb"^REMARK VINA RESULT:\s+(\S+)", raw, re.MULTILINE)
        log = prefix.with_suffix(
            ".counts" if row["arm"] == "upstream" else ".log"
        ).read_text()
        for k in COUNT:
            actual = list(map(int, re.findall(r"\b" + k + r"=(\d+)", log)))
            assert len(actual) == (1 if k.startswith("explicit") else trajectories)
            assert sum(actual) == int(row[k])
        assert sum(int(row[k]) for k in COUNT[:3]) == int(row["joint"])
        assert (
            len(scores) == int(row["modes"])
            and scores[0].decode() == row["first_score"]
        )
        coords = b"\n".join(
            x
            for x in raw.split(b"ENDMDL", 1)[0].splitlines()
            if x.startswith((b"ATOM", b"HETATM"))
        )
        row["pool_sha256"], row["first_coords_sha256"] = digest(raw), digest(coords)
        records[int(row["cap"]), row["arm"], row["receptor"], row["compound"]] = row
    assert len(records) == len(runs) == len(caps) * 3 * 2 * 47
    write(work / "runs.csv", runs)
    cost_rows = []
    for budget, receiver, arm in itertools.product(caps, (*RECEPTORS, "both"), ARMS):
        recs = RECEPTORS if receiver == "both" else (receiver,)
        pairs = [
            (records[budget, arm, rec, c], records[budget, "parent", rec, c])
            for rec, c in itertools.product(recs, compounds)
        ]
        cost_rows.append(
            {
                "cap": budget,
                "receptor": receiver,
                "arm": arm,
                "n": len(pairs),
                "wall_s": sum(int(r["elapsed_ns"]) for r, _ in pairs) / 1e9,
                "cpu_s": round(sum(float(r["cpu_s"]) for r, _ in pairs), 6),
                "joint": sum(int(r["joint"]) for r, _ in pairs),
                "local_calls": sum(sum(int(r[k]) for k in COUNT[3:]) for r, _ in pairs),
                "pool_changes": sum(
                    r["pool_sha256"] != p["pool_sha256"] for r, p in pairs
                ),
                "first_coord_changes": sum(
                    r["first_coords_sha256"] != p["first_coords_sha256"]
                    for r, p in pairs
                ),
                "first_score_changes": sum(
                    r["first_score"] != p["first_score"] for r, p in pairs
                ),
            }
        )
    write(work / "costs.csv", cost_rows)

    controls = read(work / "fixed.csv")
    assert {r["key"] for r in controls} == {
        r["key"] for r in runs if r["arm"] == "parent"
    }
    for r in controls:
        raw = (work / "raw" / (r["key"] + ".pdbqt")).read_bytes()
        assert len(r["terms"].split()) == 10 and r["pose_sha256"] == digest(
            raw.split(b"ENDMDL", 1)[0] + b"ENDMDL\n"
        )

    summary, details, scaffolds = [], [], []
    for budget, receptor in itertools.product(caps, (*RECEPTORS, "mean")):
        recs = RECEPTORS if receptor == "mean" else (receptor,)
        scores = {"size": {c: -int(labels[c]["heavy"]) for c in compounds}}
        for arm in ARMS:
            scores[arm] = {
                c: sum(Decimal(records[budget, arm, r, c]["first_score"]) for r in recs)
                / len(recs)
                for c in compounds
            }
        edges = []
        for i, j in itertools.combinations(compounds, 2):
            if kd[j] < 3 * kd[i]:
                continue
            edge = {
                "cap": budget,
                "receptor": receptor,
                "strong": i,
                "weak": j,
                "fold10": kd[j] >= 10 * kd[i],
                "matched": all(
                    labels[i][k] == labels[j][k] for k in ("heavy", "charge")
                ),
                "gap_ln_kd": float((kd[j] / kd[i]).ln()),
            }
            for arm, values in scores.items():
                edge[arm] = (values[i] < values[j]) + 0.5 * (values[i] == values[j])
            edges.append(edge)
        assert len(edges) == 403 and sum(e["matched"] for e in edges) == 35
        details.extend(edges)
        for group in sorted({r["scaffold"] for r in manifest}) + ["between_scaffolds"]:
            if group == "between_scaffolds":
                selected = [
                    e
                    for e in edges
                    if labels[e["strong"]]["scaffold"] != labels[e["weak"]]["scaffold"]
                ]
            else:
                selected = [
                    e
                    for e in edges
                    if labels[e["strong"]]["scaffold"]
                    == labels[e["weak"]]["scaffold"]
                    == group
                ]
            scaffolds.append(
                dict(
                    cap=budget,
                    receptor=receptor,
                    scaffold=group,
                    pairs=len(selected),
                    **{arm: sum(e[arm] for e in selected) for arm in scores},
                )
            )
        for arm, values in scores.items():
            preds = [float(values[c]) for c in compounds]
            measured = [float(kd[c]) for c in compounds]
            deltas = [
                sum(
                    e[arm] - e["parent"]
                    for e in edges
                    if c not in (e["weak"], e["strong"])
                )
                for c in compounds
            ]
            fold10 = [e for e in edges if e["fold10"]]
            summary.append(
                {
                    "cap": budget,
                    "receptor": receptor,
                    "arm": arm,
                    "n_pairs": len(edges),
                    "credit": sum(e[arm] for e in edges),
                    "matched_credit": sum(e[arm] for e in edges if e["matched"]),
                    "tenfold_pairs": len(fold10),
                    "tenfold_credit": sum(e[arm] for e in fold10),
                    "spearman": spearmanr(preds, measured).statistic,
                    "kendall": kendalltau(preds, measured).statistic,
                    "leave_one_compound_delta_min": min(deltas),
                    "leave_one_compound_delta_max": max(deltas),
                }
            )
    write(work / "ranks.csv", summary)
    write(work / "pairs.csv", details)
    write(work / "scaffolds.csv", scaffolds)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=("run", "audit"))
    p.add_argument("work", type=Path)
    p.add_argument("--trajectories", choices=(1, 4), type=int, default=1)
    for name in (*ARMS, "upstream-counts"):
        p.add_argument("--" + name, type=Path)
    p.add_argument("--cpu", type=int, default=min(os.sched_getaffinity(0)))
    args = p.parse_args()
    if args.action == "run":
        if any(getattr(args, k) is None for k in (*ARMS, "upstream_counts")):
            p.error("run requires all four executables")
        execute(args)
        fixed(args)
    audit(args.work)

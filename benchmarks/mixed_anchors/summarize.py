"""Compare released poses, candidate coverage, and matched Ki/Kd order within each panel."""

import argparse
import itertools
import sys
from collections import defaultdict
from pathlib import Path
from statistics import median
from dock import ROOT, SEEDS, cases, read, write

sys.path.insert(0, str(ROOT / "benchmarks/directional_ion_pair"))
from evaluate import panel as ion_panel
from reference import header
from contrasts import unique_resampling


def ion_rows(base, runs):
    work = base / "fdv-ionpair"
    bench = ROOT / "benchmarks/directional_ion_pair/results"
    labels = {
        "fxa": {
            r["ligand"]: float(r["ki_nm"])
            for r in read(ROOT / "benchmarks/affinity_partition/fxa_series.csv")
            if r["included"] == "1"
        },
        "trypsin": {
            r["ligand"]: float(r["kd_uM"]) for r in read(work / "trypsin/series.csv")
        },
    }
    sizes = {
        r["ligand"]: (r["heavy"], r["charge"])
        for r in read(ROOT / "benchmarks/affinity_partition/fxa_series.csv")
        if r["included"] == "1"
    }
    geom = {}

    def append(rows):
        for row in rows:
            if row["method"] == "release":
                key = (
                    row["target"],
                    row["receptor"],
                    row["ligand"],
                    int(row["seed"]),
                    row["arm"],
                )
                geom[key] = row

    append(read(bench / "cells.csv"))
    append(read(bench / "inverted_cells.csv"))
    for name in ("fxa", "trypsin", "inverted"):
        target = "trypsin" if name == "trypsin" else "fxa"
        folders = (
            {"3PTB": work / "trypsin", "1BTY": work / "trypsin"}
            if target == "trypsin"
            else {"1FJS": base / "fdv-fxa1", "2P16": base / "fdv-fxa"}
        )
        _, rows = ion_panel(
            target,
            labels[target],
            runs / name,
            work / "raw",
            folders,
            SEEDS[:2] if name == "inverted" else SEEDS,
            ("mixed",),
        )
        if name == "inverted":
            for row in rows:
                row["arm"] = "mixed_inverted"
        append(rows)
    result = []
    for name, rec, lig, seed, _, _, _, _, native, guided in cases(base):
        if name == "cross":
            continue
        target = "trypsin" if name == "trypsin" else "fxa"
        is_wrong = name == "inverted"
        pots = {}
        for arm, path in (
            ("native", native.with_suffix(".pdbqt")),
            ("inverted" if is_wrong else "shell", guided.with_suffix(".pdbqt")),
            (
                "mixed_inverted" if is_wrong else "mixed",
                runs / name / f"{rec}_{lig}_{seed}_mixed.pdbqt",
            ),
        ):
            pots[arm] = header(path)
        ratios = [(total - unbound) / score for score, total, unbound in pots.values()]
        rot = round((median(ratios) - 1) / 0.05846 * 2) / 2
        denominator = 1 + 0.05846 * rot
        for arm, (score, total, unbound) in pots.items():
            assert abs((total - unbound) / denominator - score) < 0.0015
            g = geom[target, rec, lig, seed, arm]
            result.append(
                dict(
                    panel=name,
                    receptor=rec,
                    ligand=lig,
                    seed=seed,
                    arm=arm,
                    measurement=labels[target][lig],
                    reported=score,
                    fixed_reference=(total - pots["native"][2]) / denominator,
                    top_bidentate=int(g["bidentate"]),
                    oracle_bidentate=int(g["oracle_bidentate"]),
                    rms_top=g["heavy_rms"],
                    rms_best=g["oracle_rms"],
                )
            )
    return result, sizes


def scaffold_rows(base, runs):
    output = []
    for panel in ("brd4", "t4l", "hcaii"):
        for folder, arms in (
            (base / ("fdv-map-" + panel), ("native", "multi")),
            (runs / panel, ("mixed",)),
        ):
            output.extend(
                dict(
                    panel=panel,
                    receptor=r["receptor"],
                    ligand=r["ligand"],
                    seed=r["seed"],
                    arm=r["arm"],
                    rms_top=r["rmsd"],
                    rms_best=r["best_rmsd"],
                    top=int(r["strict"]),
                    oracle=int(r["strict_oracle"]),
                )
                for r in read(folder / "cells.csv")
                if r["arm"] in arms and r["method"] == "released"
            )
    return output


def group_by(rs, *keys):
    groups = defaultdict(list)
    for r in rs:
        groups[tuple(r[k] for k in keys)].append(r)
    return sorted(groups.items())


def affinity(rows, sizes):
    edges = []
    for (panel, arm), group in group_by(rows, "panel", "arm"):
        measurements = {r["ligand"]: r["measurement"] for r in group}
        for rec in ["all", *sorted({r["receptor"] for r in group})]:
            blocks = [
                b
                for _, b in group_by(
                    [r for r in group if rec == "all" or r["receptor"] == rec],
                    "receptor",
                    "seed",
                )
            ]
            for method in ("reported", "fixed_reference"):
                keys = defaultdict(list)
                for block in blocks:
                    for a, b in itertools.combinations(
                        sorted(block, key=lambda r: r["measurement"]), 2
                    ):
                        if b["measurement"] >= 3 * a["measurement"]:
                            keys[a["ligand"], b["ligand"]].append(b[method] - a[method])
                for (a, b), v in keys.items():
                    gap = round(median(v), 7)
                    edges.append(
                        dict(
                            panel=panel,
                            arm=arm,
                            receptor=rec,
                            reference=method,
                            strong=a,
                            weak=b,
                            matched=int(a in sizes and sizes[a] == sizes[b]),
                            measured_ratio=measurements[b] / measurements[a],
                            correct_in_every_block=int(min(v) > 0),
                            gap=gap,
                            credit=1 if gap > 0 else 0.5 if gap == 0 else 0,
                        )
                    )
    summary = []
    for (panel, rec, ref), g in group_by(edges, "panel", "receptor", "reference"):
        native = {
            (r["strong"], r["weak"]): r["credit"] for r in g if r["arm"] == "native"
        }
        for (arm,), rs in group_by(g, "arm"):
            ds = [
                (r["strong"], r["weak"], r["credit"] - native[r["strong"], r["weak"]])
                for r in rs
            ]
            ligs = sorted({lig for a, b, _ in ds for lig in (a, b)})
            ci = (
                unique_resampling(ds, ligs)[0]
                if panel != "trypsin" and rec == "all"
                else ("", "")
            )
            summary.append(
                dict(
                    panel=panel,
                    receptor=rec,
                    reference=ref,
                    arm=arm,
                    correct=sum(r["credit"] for r in rs),
                    pairs=len(rs),
                    changed_to_right=sum(d > 0 for a, b, d in ds),
                    changed_to_wrong=sum(d < 0 for a, b, d in ds),
                    matched_credit=sum(r["credit"] for r in rs if r["matched"]),
                    matched=sum(r["matched"] for r in rs),
                    all_blocks_correct=sum(r["correct_in_every_block"] for r in rs),
                    tenfold_credit=sum(
                        r["credit"] for r in rs if r["measured_ratio"] >= 10
                    ),
                    tenfold_pairs=sum(r["measured_ratio"] >= 10 for r in rs),
                    ligand_boot_low=ci[0],
                    ligand_boot_high=ci[1],
                )
            )
    return edges, summary


def scaffold_ranks(base, runs):
    ranks = []
    for p in ("brd4", "t4l", "hcaii"):
        for source, arms in (
            (base / ("fdv-map-" + p), ("native", "multi")),
            (runs / p, ("mixed",)),
        ):
            pairs = read(source / "pairs.csv")
            for arm in arms:
                raw = [r for r in pairs if r["arm"] == arm]
                for rec in ("all", *sorted({r["receptor"] for r in raw})):
                    groups = group_by(
                        [r for r in raw if rec == "all" or r["receptor"] == rec],
                        "strong",
                        "weak",
                    )
                    gaps = [
                        round(median(float(row["gap"]) for row in values), 7)
                        for _, values in groups
                    ]
                    ranks.append(
                        dict(
                            panel=p,
                            arm=arm,
                            receptor=rec,
                            correct=sum(
                                1 if g > 0 else 0.5 if g == 0 else 0 for g in gaps
                            ),
                            pairs=len(gaps),
                        )
                    )
    return ranks


def evaluate(base, runs, out):
    out.mkdir(exist_ok=True, parents=True)
    ion, sizes = ion_rows(base, runs)
    edges, summary = affinity(ion, sizes)
    write(out / "ion_cells.csv", ion)
    write(out / "edges.csv", edges)
    write(out / "affinity.csv", summary)
    structure = scaffold_rows(base, runs)
    structure.extend(
        dict(
            panel="cross",
            receptor=r["receptor"],
            ligand="DX9",
            seed=r["seed"],
            arm=r["arm"],
            rms_top=r["top"],
            rms_best=r["oracle"],
            top=int(float(r["top"]) <= 2),
            oracle=int(float(r["oracle"]) <= 2),
        )
        for r in read(runs / "cross_cells.csv")
    )
    for r in ion:
        if r["panel"] == "trypsin" and r["rms_top"] != "":
            structure.append(
                dict(
                    panel="benzamidine",
                    receptor=r["receptor"],
                    ligand=r["ligand"],
                    seed=r["seed"],
                    arm=r["arm"],
                    rms_top=r["rms_top"],
                    rms_best=r["rms_best"],
                    top=int(float(r["rms_top"]) <= 2),
                    oracle=int(float(r["rms_best"]) <= 2),
                )
            )
    write(out / "structure_cells.csv", structure)
    table = []
    for (panel, arm), g in group_by(structure, "panel", "arm"):
        for rec in ["all", *sorted({r["receptor"] for r in g})]:
            r = [v for v in g if rec == "all" or v["receptor"] == rec]
            table.append(
                dict(
                    panel=panel,
                    arm=arm,
                    receptor=rec,
                    cells=len(r),
                    top=sum(v["top"] for v in r),
                    oracle=sum(v["oracle"] for v in r),
                    median_rms=median(float(v["rms_top"]) for v in r),
                    top_2A=sum(float(v["rms_top"]) <= 2 for v in r),
                    pool_2A=sum(float(v["rms_best"]) <= 2 for v in r),
                )
            )
    write(out / "structures.csv", table)
    write(out / "scaffold_affinity.csv", scaffold_ranks(base, runs))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for n in ("base", "runs", "out"):
        p.add_argument(n, type=Path)
    a = p.parse_args()
    evaluate(a.base, a.runs, a.out)

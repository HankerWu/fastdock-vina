"""Paired centre-energy and restricted-volume comparisons on the same PDBQT pools."""

import argparse
import itertools
import statistics
from collections import defaultdict
from pathlib import Path
from run import table, write


def median(values):
    return statistics.median(values)


def credit(x):
    return 1 if x > 0 else 0.5 if x == 0 else 0


def evaluate(folder):
    cells, pairs, checks = [], [], []
    for panel in (
        "hcaii",
        "brd4",
        "t4l",
        "alkyl",
        "fxa",
        "fxa_pair",
        "trypsin",
        "bovine_ca",
    ):
        if not (folder / (panel + ".csv")).is_file():
            continue
        groups = defaultdict(list)
        rows = table(folder / (panel + ".csv"))
        err = [abs(float(r["window_delta"]) - float(r["half_delta"])) for r in rows]
        checks.append(
            dict(
                panel=panel,
                poses=len(rows),
                median_half_error=median(err),
                p95_half_error=sorted(err)[int(0.95 * (len(err) - 1))],
                median_effective=median([float(r["effective_samples"]) for r in rows]),
                max_inter_delta=max(
                    abs(float(r["source_inter"]) - float(r["inter"])) for r in rows
                ),
                max_intra_delta=max(
                    abs(float(r["source_intra"]) - float(r["intra"])) for r in rows
                ),
            )
        )
        for r in rows:
            groups[r["receptor"], r["ligand"], r["seed"], r["arm"]].append(r)
        for (rec, lig, seed, arm), g in groups.items():

            def order(r):
                return r["slot"] != "common", int(r["model"])

            if panel == "hcaii":
                raw = lambda r: float(r["source_inter"]) + float(r["source_intra"])
                winner = min(g, key=lambda r: (raw(r), order(r)))
                unbound = float(winner["source_intra"])
                for r in g:
                    r["center"] = (raw(r) - unbound) * float(r["slope"])
            else:
                for r in g:
                    r["center"] = float(r["source_affinity"])
            center = min(g, key=lambda r: (r["center"], order(r)))
            for method, key, fixed in (
                ("native", None, False),
                ("window", "window_delta", False),
                ("fixed_winner", "window_delta", True),
                ("half", "half_delta", False),
                ("uniform", "uniform_delta", False),
            ):

                def energy(r):
                    return r["center"] + (float(r[key]) if key else 0.0)

                top = center if fixed else min(g, key=lambda r: (energy(r), order(r)))
                hit = int(float(top["rmsd"]) <= 2) if top["rmsd"] else ""
                strict = (hit * int(top["sign_ok"])) if hit != "" else ""
                cells.append(
                    dict(
                        panel=panel,
                        receptor=rec,
                        ligand=lig,
                        seed=seed,
                        arm=arm,
                        method=method,
                        slot=top["slot"],
                        model=top["model"],
                        score=energy(top),
                        rmsd=top["rmsd"],
                        hit=hit,
                        strict=strict,
                        measurement=top["measurement"],
                        same_choice=int(order(top) == order(center)),
                    )
                )
    blocks = defaultdict(list)
    for r in cells:
        if r["measurement"] != "":
            blocks[r["panel"], r["receptor"], r["seed"], r["arm"], r["method"]].append(
                r
            )
    for (panel, rec, seed, arm, method), g in blocks.items():
        for s, w in itertools.combinations(
            sorted(g, key=lambda r: float(r["measurement"])), 2
        ):
            if float(w["measurement"]) < 3 * float(s["measurement"]):
                continue
            gap = round(w["score"] - s["score"], 9)
            pairs.append(
                dict(
                    panel=panel,
                    receptor=rec,
                    seed=seed,
                    arm=arm,
                    method=method,
                    strong=s["ligand"],
                    weak=w["ligand"],
                    gap=gap,
                    credit=credit(gap),
                )
            )
    summ = []
    for panel, arm, method in sorted(
        {(r["panel"], r["arm"], r["method"]) for r in cells}
    ):
        g = [
            r
            for r in cells
            if (r["panel"], r["arm"], r["method"]) == (panel, arm, method)
        ]
        edges = [
            r
            for r in pairs
            if (r["panel"], r["arm"], r["method"]) == (panel, arm, method)
        ]
        chemical = defaultdict(list)
        for edge in edges:
            chemical[edge["strong"], edge["weak"]].append(edge["gap"])
        known = [r for r in g if r["hit"] != ""]
        summ.append(
            dict(
                panel=panel,
                arm=arm,
                method=method,
                cells=len(g),
                known_pose_cells=len(known),
                identical=sum(r["same_choice"] for r in g),
                top_hits=sum(r["hit"] for r in known) if known else "",
                joint_hits=sum(r["strict"] for r in known) if known else "",
                median_rmsd=median([float(r["rmsd"]) for r in known]) if known else "",
                repeated=sum(r["credit"] for r in edges),
                comparisons=len(edges),
                unique=sum(credit(median(v)) for v in chemical.values()),
                edges=len(chemical),
            )
        )
    for n, rs in (
        ("cells", cells),
        ("pairs", pairs),
        ("summary", summ),
        ("numerics", checks),
    ):
        write(folder / (n + ".csv"), rs)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("folder", type=Path)
    evaluate(p.parse_args().folder)

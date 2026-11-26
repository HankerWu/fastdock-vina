"""Paired continuous-rotor docking with and without a soft search-only amide term."""

import argparse
import csv
import itertools
import json
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from statistics import median

import numpy as np
from evaluate import (
    ROOT,
    ccd_vectors,
    error_components,
    explicit_scores,
    handed,
    parse_atoms,
    read_models,
    write,
)
from run import dock, mapped_ligand_bonds, sources

# isort: split
# run installs the repository tools directory on the import path.
from amide_rotamers import active_amide_sites, search_amide_records
from rank_amide_pool import rank

SEEDS = (1618033, 1414213)
SEARCH_BARRIER = 5.0
RANK_BARRIERS = (0, 5, 15)


def pools(case, work, *, seeds=SEEDS):
    panel, rec, lig, states, receptor, box, anchor = case
    for state, path in states.items():
        text = path.read_text()
        eligible = bool(active_amide_sites(text, "all-arms"))
        prefix = work / "inputs" / f"{panel}_{lig}_{state}"
        bonds = prefix.with_suffix(".bonds")
        sitefile = prefix.with_suffix(".amide")
        for seed, slot in itertools.product(seeds, (0, 1)):
            for search in ("baseline", "guided") if eligible else ("baseline",):
                name = f"{rec}_{lig}_{state}_{seed}_{slot}.pdbqt"
                yield (
                    (state, seed, slot, search),
                    work / search / name,
                    (
                        receptor,
                        path,
                        box,
                        anchor,
                        bonds,
                        seed + 1009 * slot,
                    ),
                    sitefile if search == "guided" else None,
                )


def generate(binary, cases, work, workers, which):
    for subdir in ("inputs", "baseline", "guided"):
        (work / subdir).mkdir(parents=True, exist_ok=True)
    tasks = []
    for case in cases:
        panel, _, lig, states, *_ = case
        for state, path in states.items():
            text = path.read_text()
            prefix = work / "inputs" / f"{panel}_{lig}_{state}"
            prefix.with_suffix(".bonds").write_text(
                "".join(f"{i} {j}\n" for i, j in mapped_ligand_bonds(text))
            )
            if active_amide_sites(text, "all-arms"):
                prefix.with_suffix(".amide").write_text(
                    search_amide_records(text, SEARCH_BARRIER)
                )
        for (_, _, _, search), pool, args, sites in pools(case, work):
            if which in ("both", search):
                tasks.append((pool, args, sites))

    def run_one(task):
        pool, args, sites = task
        dock(binary, (pool, args), search_amide=sites)

    with ThreadPoolExecutor(workers) as executor:
        list(executor.map(run_one, tasks))
    print(f"Completed {len(tasks)} docking processes")


def score_case(binary, work, case, affinities, *, seeds=SEEDS):
    panel, rec, lig, states, receptor, _, _ = case
    original = states["orig"]
    orig_text = original.read_text()
    if panel == "ca":
        prep = original.parent
        meta = json.loads((prep / "metadata.json").read_text())[lig]
        ref = meta["references"][rec]["coords"]
        sign_ids, sign = ccd_vectors(prep.parent, lig, meta)
        ki = meta["Ki_nM"]
    else:
        ki = affinities[lig]
    rows = []
    for (state, seed, slot, search), path, _, _ in pools(case, work, seeds=seeds):
        text = path.read_text()
        models = read_models(text.splitlines())
        terms = {
            r["model"]: r
            for r in rank(states[state].read_text(), text, 1.0, "all-arms")
        }
        energies = explicit_scores(
            binary,
            original,
            receptor,
            work / "inputs" / f"{panel}_{lig}_orig.bonds",
            orig_text,
            models,
        )
        for (number, lines), energy in zip(models.items(), energies, strict=True):
            if int(energy["model"]) != number:
                raise ValueError(f"Energy/model mismatch in {path}")
            rmsd = correct_hand = None
            if panel == "ca":
                atoms = parse_atoms(lines)
                rmsd = error_components(prep, lig, atoms, ref)["rmsd"]
                correct_hand = int(
                    handed([np.array(atoms[i].xyz) for i in sign_ids]) == sign
                )
            rows.append(
                {
                    "panel": panel,
                    "rec": rec,
                    "ligand": lig,
                    "state": state,
                    "seed": seed,
                    "slot": slot,
                    "search": search,
                    "model": number,
                    "measurement": ki,
                    "inter": float(energy["inter"]),
                    "intra": float(energy["intra"]),
                    "slope": float(energy["slope"]),
                    "twist": terms[number]["amide_penalty"],
                    "rmsd": rmsd,
                    "hand": correct_hand,
                }
            )
    return rows


def hit(row):
    return row["rmsd"] <= 2 and row["hand"] == 1


def ranking_energy(row, unbound, barrier):
    native = (row["inter"] + row["intra"] - unbound) * row["slope"]
    return native + barrier * row["twist"]


def summarize(rows):
    groups = defaultdict(list)
    for r in rows:
        groups[r["panel"], r["rec"], r["ligand"], r["seed"]].append(r)
    cells = []
    for (panel, rec, lig, seed), group in groups.items():
        baseline = [r for r in group if r["search"] == "baseline"]
        guided = [r for r in group if r["search"] == "guided"] or baseline
        unbound = min(baseline, key=lambda r: r["inter"] + r["intra"])["intra"]
        arms = {"baseline": baseline, "guided": guided}
        for slot in (0, 1):
            arms[f"guided_slot{slot}"] = [r for r in baseline if r["slot"] != slot] + [
                r for r in guided if r["slot"] == slot
            ]
        for arm, candidates in arms.items():
            for barrier in RANK_BARRIERS:
                # A common unbound energy matters when comparing different search arms.
                def sort_key(row, unbound=unbound, barrier=barrier):
                    return (
                        ranking_energy(row, unbound, barrier),
                        row["state"],
                        row["slot"],
                        row["model"],
                    )

                top = min(candidates, key=sort_key)
                cells.append(
                    {
                        "panel": panel,
                        "receptor": rec,
                        "ligand": lig,
                        "seed": seed,
                        "arm": arm,
                        "barrier": barrier,
                        "score": ranking_energy(top, unbound, barrier),
                        "measurement": top["measurement"],
                        "unbound": unbound,
                        "state": top["state"],
                        "slot": top["slot"],
                        "model": top["model"],
                        "rmsd": top["rmsd"],
                        "strict": int(hit(top)) if panel == "ca" else None,
                        "oracle": int(any(map(hit, candidates)))
                        if panel == "ca"
                        else None,
                    }
                )
    return cells


def endpoints(cells):
    totals = []
    for (panel, arm, barrier), group in itertools.groupby(
        sorted(cells, key=lambda r: (r["panel"], r["arm"], r["barrier"])),
        key=lambda r: (r["panel"], r["arm"], r["barrier"]),
    ):
        rows = list(group)
        blocks = defaultdict(list)
        for r in rows:
            blocks[r["receptor"], r["seed"]].append(r)
        gaps = defaultdict(list)
        for block in blocks.values():
            for strong, weak in itertools.permutations(block, 2):
                if weak["measurement"] >= 3 * strong["measurement"]:
                    gaps[strong["ligand"], weak["ligand"]].append(
                        weak["score"] - strong["score"]
                    )
        cross = [r for r in rows if r["ligand"] != r["receptor"]]
        totals.append(
            {
                "panel": panel,
                "arm": arm,
                "barrier": barrier,
                "cells": len(rows),
                "top": sum(r["strict"] for r in rows) if panel == "ca" else None,
                "oracle": sum(r["oracle"] for r in rows) if panel == "ca" else None,
                "cross_top": sum(r["strict"] for r in cross) if panel == "ca" else None,
                "cross_oracle": sum(r["oracle"] for r in cross)
                if panel == "ca"
                else None,
                "unique_score": sum(
                    (1 + np.sign(median(v))) / 2 for v in gaps.values()
                ),
                "pairs": len(gaps),
            }
        )
    return totals


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("dock", "score"))
    for name in ("binary", "base", "work"):
        parser.add_argument(name, type=Path)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--panel", choices=("ca", "fxa"))
    parser.add_argument(
        "--search", choices=("both", "baseline", "guided"), default="both"
    )
    a = parser.parse_args()
    binary = a.binary.resolve()
    cases = [
        case
        for case in sources(a.base)
        if case[0] in ("ca", "fxa") and (a.panel is None or case[0] == a.panel)
    ]
    if a.action == "dock":
        generate(binary, cases, a.work, a.jobs, a.search)
    else:
        if a.search != "both":
            parser.error("paired scoring requires both search arms")
        with (ROOT / "benchmarks/affinity_partition/fxa_series.csv").open() as f:
            affinities = {
                r["ligand"]: float(r["ki_nm"])
                for r in csv.DictReader(f)
                if r["included"] == "1"
            }
        with ThreadPoolExecutor(a.jobs) as executor:
            rows = list(
                itertools.chain.from_iterable(
                    executor.map(
                        lambda c: score_case(binary, a.work, c, affinities),
                        cases,
                    )
                )
            )
        out = a.work / ("results-" + a.panel if a.panel else "results")
        out.mkdir(exist_ok=True)
        cells = summarize(rows)
        write(out / "modes.csv", rows)
        write(out / "cells.csv", cells)
        write(out / "summary.csv", endpoints(cells))
        print(f"Scored {len(rows)} poses")


if __name__ == "__main__":
    main()

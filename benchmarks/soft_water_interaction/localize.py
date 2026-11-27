"""For every changed winner, resolve the water penalty advantage by deposited site."""

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from evaluate import (
    Site,
    aligned_receivers,
    alignment,
    cross_sites,
    models,
    pdb_lines,
    xyz,
)


def localize(work, cells, panel):
    meta = json.loads((work / "metadata.json").read_text())
    rows = []
    for q, r in itertools.permutations(sorted(meta), 2):
        a, b = meta[q], meta[r]
        if a["target"] != b["target"]:
            continue
        pdb = pdb_lines(work / "raw" / f"{r}.pdb")
        protein = cKDTree(
            [xyz(l) for l in pdb if l.startswith("ATOM") and l[21] == "A"]
        )
        sites = cross_sites(pdb, protein, np.array(b["center"]))
        pair = aligned_receivers(a, b) if panel == "transfer" else (a, b)
        rot, trans, _, _ = alignment(*pair)
        references = [cKDTree(x @ rot + trans) for x in np.asarray(a["mapped"])]
        prepared = []
        for site, center, b_iso in sites:
            field = Site.prepare(center, b_iso, protein)
            reference = [field.evaluate(tree)[1] for tree in references]
            prepared.append((site, field, min(reference), max(reference)))
        pool_cells = cells[(cells["query"] == q) & (cells.receiver == r)]
        for cell in pool_cells.itertuples():
            p = work / f"runs_{cell.cap}/{q}_{r}_{cell.seed}_{cell.search}.pdbqt"
            pool = {model: (dry, tree) for model, dry, tree in models(p)}
            first, winner = pool[1], pool[cell.model]
            for name, field, low, high in prepared:
                original = field.evaluate(first[1])[1]
                selected = field.evaluate(winner[1])[1]
                rows.append(
                    dict(
                        panel=panel,
                        target=a["target"],
                        cap=cell.cap,
                        query=q,
                        receiver=r,
                        seed=cell.seed,
                        search=cell.search,
                        winner=cell.model,
                        dry_gap=winner[0] - first[0],
                        site=name,
                        penalty_original=original,
                        penalty_winner=selected,
                        penalty_difference=selected - original,
                        min_reference_penalty=low,
                        max_reference_penalty=high,
                    )
                )
    return rows


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("historical", type=Path)
    p.add_argument("transfer", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    cells = pd.read_csv(a.out / "cells.csv")
    changes = cells.loc[(cells.ranker == "water") & (cells.model != 1)]
    panels = [(a.historical, "historical"), (a.transfer, "transfer")]
    rows = [
        row
        for work, panel in panels
        for row in localize(work, changes[changes.panel == panel], panel)
    ]
    pd.DataFrame(rows).to_csv(a.out / "changed_winner_sites.csv", index=False)

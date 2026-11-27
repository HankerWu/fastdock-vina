"""Resolve each GAA pool's dry and corrected winners into aromatic residue contributions."""

from argparse import ArgumentParser
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from study import CHPi, WEIGHT, cross
from model import RING_NAMES, kernel


def decompose(inputs, results, out):
    rows = []
    cells = pd.read_csv(results / "cells.csv")
    scores = pd.read_csv(results / "modes.csv")
    models = {}
    for key, group in cells.groupby(
        ["panel", "query", "receptor", "seed", "cap", "search"]
    ):
        panel, query, receptor, seed, cap, search = key
        label = dict(zip(group.columns[:5], key[:5]), search=search)
        protein = (inputs / receptor / "receptor.pdbqt").read_text()
        if (query, receptor) not in models:
            rings = defaultdict(list)
            for line in protein.splitlines():
                res = line[17:20]
                if (
                    line.startswith(("ATOM  ", "HETATM"))
                    and line.split()[-1] in ("C", "A")
                    and line[12:16].strip() in RING_NAMES.get(res, ())
                ):
                    rings[line[17:27].strip()].append(cross.coords(line))
            model = CHPi(protein, (inputs / query / "protonated.pdbqt").read_text())
            models[query, receptor] = model, rings
        model, rings = models[query, receptor]
        poses = (
            (results / "pools" / ("_".join(map(str, key)) + ".pdbqt"))
            .read_text()
            .split("ENDMDL")
        )
        native = group[group.rerank == "native"].iloc[0]
        corrected = group[group.rerank == "field"].iloc[0]
        if native.winner == corrected.winner:
            continue
        terms = []
        for winner in (native.winner, corrected.winner):
            mode = scores[
                (scores.panel == panel)
                & (scores["query"] == query)
                & (scores.receptor == receptor)
                & (scores.cap == cap)
                & (scores.seed == seed)
                & (scores.search == search)
                & (scores["rank"] == winner)
            ].iloc[0]
            atoms = {
                int(l[6:11]): cross.coords(l)
                for l in poses[winner - 1].splitlines()
                if l.startswith(("ATOM  ", "HETATM"))
            }
            donors = [atoms[i] for i in model.l_ch]
            by_res = {
                res: WEIGHT * mode.factor * float(kernel(cdist(donors, coords)).sum())
                for res, coords in rings.items()
            }
            if not np.isclose(
                sum(by_res.values()), mode.corrected - mode.native, atol=1e-10, rtol=0
            ):
                raise ValueError("Residue contributions do not reconstruct the score")
            terms.append(by_res)
        for res in sorted(rings):
            if terms[0][res] or terms[1][res]:
                rows.append(
                    dict(
                        **label,
                        residue=res,
                        dry_winner=native.winner,
                        new_winner=corrected.winner,
                        dry_term=terms[0][res],
                        new_term=terms[1][res],
                        delta=terms[1][res] - terms[0][res],
                    )
                )
    pd.DataFrame(rows).to_csv(out, index=False)


if __name__ == "__main__":
    p = ArgumentParser(description=__doc__)
    for name in ("inputs", "results", "out"):
        p.add_argument(name, type=Path)
    decompose(**vars(p.parse_args()))

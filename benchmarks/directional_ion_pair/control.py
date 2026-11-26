"""Compare the reversed-head shell search to the declared pharmacophore, at unchanged budget."""

import argparse
from pathlib import Path

from evaluate import compare, panel
from experiment import ROOT, SEEDS, table, write


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "fjs", "p16", "out"):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    meta = [
        r
        for r in table(ROOT / "benchmarks/affinity_partition/fxa_series.csv")
        if int(r["included"])
    ]
    labels = {"fxa": {r["ligand"]: float(r["ki_nm"]) for r in meta}}
    sizes = {r["ligand"]: (r["heavy"], r["charge"]) for r in meta}
    _, cells = panel(
        "fxa",
        labels["fxa"],
        a.work / "inverted",
        a.work / "raw",
        {"1FJS": a.fjs, "2P16": a.p16},
        SEEDS,
        ("inverted",),
    )
    pairs = compare(cells, labels)
    write(a.out / "inverted_cells.csv", cells)
    group = [r for r in cells if r["method"] == "release"]
    eq = [r for r in pairs if r["method"] == "release" and r["ratio"] >= 3]
    matched = [r for r in eq if sizes[r["strong"]] == sizes[r["weak"]]]
    print(
        "Reversed-head initial seeds:",
        sum(r["bidentate"] for r in group),
        "of",
        len(group),
        "native top aryl contacts;",
        sum(r["credit"] for r in eq),
        "/",
        len(eq),
        "unique pairs;",
        sum(r["credit"] for r in matched),
        "/",
        len(matched),
        "matched",
    )
    print(
        next(
            r
            for r in eq
            if r["strong"] == "CHEMBL177142" and r["weak"] == "CHEMBL369161"
        )
    )


if __name__ == "__main__":
    main()

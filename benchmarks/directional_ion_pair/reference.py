"""Hold the unbound Vina reference fixed across search arms on the identical top poses."""

import argparse
import itertools
from pathlib import Path
from statistics import median
from experiment import ARMS, ROOT, SEEDS, table, write
from evaluate import compare
from contrasts import EXTRA
from trypsin import SERIES, RECEPTORS


def header(path):
    lines = path.read_text().split("ENDMDL", 1)[0].splitlines()

    def value(prefix):
        return float(
            next(l[len(prefix) :].split()[0] for l in lines if l.startswith(prefix))
        )

    return (
        value("REMARK VINA RESULT:"),
        value("REMARK INTER + INTRA:"),
        value("REMARK UNBOUND:"),
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("work", "fjs", "p16", "out"):
        p.add_argument(k, type=Path)
    a = p.parse_args()
    r = [
        r
        for r in table(ROOT / "benchmarks/affinity_partition/fxa_series.csv")
        if int(r["included"])
    ]
    labels = {
        "fxa": {k["ligand"]: float(k["ki_nm"]) for k in r},
        "trypsin": {k: 100 / v for k, (s, v) in SERIES.items()},
    }
    sizes = {v["ligand"]: (v["heavy"], v["charge"]) for v in r}
    raw = []
    for target, recs in (("fxa", ("1FJS", "2P16")), ("trypsin", RECEPTORS)):
        for rec, seed, lig in itertools.product(recs, (*SEEDS, *EXTRA), labels[target]):
            source = (a.fjs if rec == "1FJS" else a.p16) / f"{lig}_{seed}.pdbqt"
            for arm in ARMS:
                path = (
                    source
                    if arm == "native" and target == "fxa" and seed in SEEDS
                    else a.work / target / f"{rec}_{lig}_{seed}_{arm}.pdbqt"
                )
                score, energy, unbound = header(path)
                raw.append(
                    dict(
                        target=target,
                        receptor=rec,
                        seed=seed,
                        ligand=lig,
                        arm=arm,
                        reported=score,
                        total=energy,
                        unbound=unbound,
                    )
                )
    by = {
        (v["target"], v["receptor"], v["seed"], v["ligand"]): v
        for v in raw
        if v["arm"] == "native"
    }
    denominator = {}
    for t in labels:
        for lig in labels[t]:
            ratios = [
                (v["total"] - v["unbound"]) / v["reported"]
                for v in raw
                if v["target"] == t and v["ligand"] == lig
            ]
            # The original scoring-function constructor uses the rounded internal coefficient.
            rrot = round((median(ratios) - 1) / 0.05846 * 2) / 2
            factor = 1 + 0.05846 * rrot
            if max(abs(v - factor) for v in ratios) > 0.001:
                raise ValueError(f"Unknown torsion factor {t} {lig} {rrot}")
            denominator[t, lig] = factor
    cells = []
    for v in raw:
        key = (v["target"], v["receptor"], v["seed"], v["ligand"])
        base = by[key]
        for mode, ref in (("arm", v["unbound"]), ("fixed_native", base["unbound"])):
            score = (v["total"] - ref) / denominator[v["target"], v["ligand"]]
            if mode == "arm" and abs(score - v["reported"]) > 0.0015:
                raise ValueError("Cannot reproduce native score")
            cells.append(dict(**v, method=mode, score=score))
    summary = []
    for split, ss in (("initial", SEEDS), ("new", EXTRA), ("all", (*SEEDS, *EXTRA))):
        edges = compare([r for r in cells if r["seed"] in ss], labels)
        for t, arm, method in itertools.product(labels, ARMS, ("arm", "fixed_native")):
            rows = [
                r
                for r in edges
                if (r["target"], r["arm"], r["method"]) == (t, arm, method)
                and r["ratio"] >= 3
            ]
            matched = [
                r for r in rows if t == "fxa" and sizes[r["strong"]] == sizes[r["weak"]]
            ]
            diag = next(
                (
                    r["predicted_gap"]
                    for r in rows
                    if (r["strong"], r["weak"]) == ("CHEMBL177142", "CHEMBL369161")
                ),
                "",
            )
            summary.append(
                dict(
                    split=split,
                    target=t,
                    arm=arm,
                    reference=method,
                    credit=sum(r["credit"] for r in rows),
                    pairs=len(rows),
                    matched_credit=sum(r["credit"] for r in matched),
                    matched=len(matched),
                    diagnostic_gap=diag,
                )
            )
    write(a.out / "reference_cells.csv", cells)
    write(a.out / "reference_summary.csv", summary)


if __name__ == "__main__":
    main()

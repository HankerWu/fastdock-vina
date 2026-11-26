"""Post-hoc sensitivity to molecule size, charge, target weighting, and removal of one ligand."""

from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


def describe(pairs, labels, dims):
    key = dims + ["stronger", "weaker"]
    wide = pairs.pivot(index=key, columns="score", values="credit").reset_index()
    same_size = wide.stronger.map(labels.heavy).eq(wide.weaker.map(labels.heavy))
    same_charge = wide.stronger.map(labels.charge).eq(wide.weaker.map(labels.charge))
    wide["matched"] = same_size & same_charge
    if "target" in wide:
        wide["panel"] = np.where(
            wide.receptor.eq("mean_receptors"), "mean_receptors", "receptor_strata"
        )
        dims = [c for c in dims if c not in ("target", "assay", "receptor")] + ["panel"]
    rows = []
    for group_key, group in wide.groupby(dims):
        for subset in ("all", "equal_size_charge"):
            s = group if subset == "all" else group[group.matched]
            if s.empty:
                continue
            for baseline in ("native", "uniform"):
                delta = s.contact - s[baseline]
                changed = delta != 0
                ligset = set(s.stronger) | set(s.weaker)
                delete_one = [
                    delta[(s.stronger != c) & (s.weaker != c)].mean() for c in ligset
                ]
                rows.append(
                    dict(
                        **dict(zip(dims, group_key)),
                        subset=subset,
                        baseline=baseline,
                        pairs=len(s),
                        base=s[baseline].sum(),
                        contact=s.contact.sum(),
                        gained_credit=delta[delta > 0].sum(),
                        lost_credit=-delta[delta < 0].sum(),
                        switched_pairs=changed.sum(),
                        delete_ligand_delta_min=np.nanmin(delete_one),
                        delete_ligand_delta_max=np.nanmax(delete_one),
                    )
                )
    return pd.DataFrame(rows)


def main():
    source = HERE.parents[1] / "geometry_comparison" / "affinity"
    main = (
        pd.read_csv(source / "panel.csv", dtype={"compound": str})
        .drop_duplicates("compound")
        .set_index("compound")
    )
    main["charge"] = main.formal_charge
    main["heavy"] = (
        pd.read_csv(HERE / "native_rotors.csv", dtype={"compound": str})
        .set_index("compound")
        .heavy_atoms
    )
    pairs = pd.read_csv(HERE / "pairs.csv.xz", dtype={"stronger": str, "weaker": str})
    results = describe(
        pairs, main, ["budget", "optimizer", "policy", "target", "assay", "receptor"]
    )
    results.to_csv(HERE / "sensitivity.csv", index=False)
    fxa = pd.read_csv(
        HERE.parents[1] / "affinity_partition" / "fxa_series.csv"
    ).set_index("ligand")
    results = describe(
        pd.read_csv(HERE / "fxa_pairs.csv.xz"), fxa, ["policy", "receptor"]
    )
    results.to_csv(HERE / "fxa_sensitivity.csv", index=False)


if __name__ == "__main__":
    main()

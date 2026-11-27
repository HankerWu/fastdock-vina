"""Rebuild absolute scores directly, verify mode/pair ordering and tabulate correlations."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, spearmanr

SOURCE = Path(__file__).resolve().parents[1]


def run(out):
    modes = pd.read_csv(SOURCE / "modes.csv.xz", dtype={"compound": str})
    corrected = pd.read_csv(SOURCE / "corrected_modes.csv", dtype={"compound": str})
    modes = pd.concat(
        [modes[~modes.compound.isin(corrected.compound)], corrected], ignore_index=True
    )

    def minima(path):
        return (
            pd.read_csv(path, dtype={"compound": str})
            .groupby("compound")[["vacuum", "solution"]]
            .min()
        )

    old = minima(SOURCE / "references.csv.xz")
    new = minima(out / "independent_references.csv.xz")
    native_inter = modes.full - modes.scale * (
        modes.bound_mmff - modes.compound.map(old.vacuum)
    )
    water = modes.bound_mmff + modes.bound_gb
    np.testing.assert_allclose(
        native_inter + modes.scale * (water - modes.compound.map(old.solution)),
        modes.gb_reselected,
        atol=1e-8,
        rtol=0,
    )
    modes["fresh_full"] = native_inter + modes.scale * (
        modes.bound_mmff - modes.compound.map(new.vacuum)
    )
    modes["fresh_gb_reselected"] = native_inter + modes.scale * (
        water - modes.compound.map(new.solution)
    )
    for _, pool in modes.groupby(["compound", "pdb", "seed"]):
        for arm in ("full", "gb_reselected"):
            before = pool.sort_values([arm, "model"]).model.to_list()
            after = pool.sort_values(["fresh_" + arm, "model"]).model.to_list()
            if before != after:
                raise ValueError("Free ligand reference changed same-ligand pose order")

    first = modes[modes.model == 1]
    edges = pd.read_csv(out / "pairs.csv.xz", dtype={"stronger": str, "weaker": str})
    bace_fxa = (
        pd.read_csv(SOURCE.parent / "strain_components/scores.csv")
        .set_index("compound")
        .ki_nm
    )
    bace_fxa = bace_fxa.groupby(level=0).first()
    lectins = pd.read_csv(
        SOURCE.parents[1] / "galectin_affinity/table1.csv", dtype={"compound": str}
    )
    stats = []
    for (assay, receptor, arm), group in edges.groupby(["assay", "receptor", "arm"]):
        per_assay = edges[edges.assay == assay]
        recs = set(per_assay.receptor) - {"mean", "both"}
        if receptor not in ("mean", "both"):
            recs = {receptor}
        lab = (
            lectins.set_index("compound")[assay + "_Kd_uM"]
            if assay in ("Gal1", "Gal3")
            else bace_fxa
        )
        # Some measured compounds have no >=3-fold neighbor. Include them in correlations.
        if assay in ("Gal1", "Gal3"):
            selected = first[first.pdb.isin(recs)]
        else:
            selected = first[first.pdb.isin(recs) & (first.assay == assay)]

        for score in ("native", arm, "fresh_" + arm):
            wide = selected.pivot(
                index="compound", columns=["pdb", "seed"], values=score
            )
            if wide.isna().any().any():
                raise ValueError("Incomplete paired receptor/seed data")
            predicted = wide.T.groupby(level="pdb").median().mean(axis=0)
            if receptor != "both":
                stats.append(
                    {
                        "assay": assay,
                        "receptor": receptor,
                        "score": score,
                        "n": len(wide),
                        "spearman": spearmanr(lab.loc[wide.index], predicted).statistic,
                        "kendall_tau_b": kendalltau(
                            lab.loc[wide.index], predicted
                        ).statistic,
                    }
                )
            if score.startswith("fresh_"):
                if assay == "CHEMBL661804":
                    gaps = np.median(
                        wide.loc[group.weaker].values - wide.loc[group.stronger].values,
                        axis=1,
                    )
                else:
                    gaps = (
                        predicted.loc[group.weaker].values
                        - predicted.loc[group.stronger].values
                    )
                np.testing.assert_allclose(gaps, group.fresh_gap, atol=1e-9, rtol=0)
    pd.DataFrame(stats).drop_duplicates().to_csv(out / "correlations.csv", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    run(parser.parse_args().out)

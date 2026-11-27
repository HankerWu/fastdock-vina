"""Compare -docking score with measured -log10(K in mol/L), one row per ligand.

Bootstrap draws sample ligands with replacement within each assay. The two arms
use identical draws and rerank all observations (including duplicates) inside
each draw. Use percentile limits as descriptive within-series uncertainty.
The threefold measured-gap filter defines a sensitivity subset and does not
imply that the underlying experiments resolve individual differences.
"""

import argparse
from collections import defaultdict
from decimal import Decimal
from itertools import combinations
from math import log10
from pathlib import Path
import re
import sys

import numpy as np
from scipy.stats import kendalltau, rankdata, spearmanr

BENCH = Path(__file__).resolve().parents[1]
SOURCE = BENCH / "ligand_hydroxyl"
sys.path.insert(0, str(BENCH / "terminal_rotor_score"))
from assay import aggregate, read, save

ARMS = ("native", "bonded")
ASSAYS = {
    "FXa": ("CHEMBL661804", "Ki", "nM"),
    "Gal3": ("Sorme_Fig4_FP", "Kd", "uM"),
    "S6K1": ("Davis_binding", "Kd", "nM"),
    "ERRg": ("Thouennon_Table1_ITC", "Kd", "nM"),
}
UNIT_M = {
    "M": 1.0,
    "mM": 1e-3,
    "uM": 1e-6,
    "µM": 1e-6,
    "μM": 1e-6,
    "nM": 1e-9,
    "pM": 1e-12,
}


def ligand_rows(target, budget, pose_rule, cells):
    """One prediction per ligand. Take median over seeds, then mean over receptors."""
    assay, endpoint, unit = ASSAYS[target]
    scores = {arm: aggregate(cells, arm) for arm in ARMS}
    labels = defaultdict(set)
    runs = defaultdict(set)
    for r in cells:
        key = r["ligand"]
        run = (r["receptor"], r["seed"])
        if run in runs[key]:
            raise ValueError(f"Duplicate ligand/receptor/seed in {target}: {key} {run}")
        labels[key].add(Decimal(r["measurement"]))
        runs[key].add(run)
    if len(set(map(frozenset, runs.values()))) != 1:
        raise ValueError(f"Unbalanced receptor/seed panel in {target}")
    rows = []
    for key in sorted(labels):
        if len(labels[key]) != 1:
            raise ValueError(f"Conflicting measurement for {key}")
        concentration = next(iter(labels[key]))
        rows.append(
            dict(
                target=target,
                assay=assay,
                budget=budget,
                pose_rule=pose_rule,
                ligand=key,
                endpoint=endpoint,
                K=concentration,
                unit=unit,
                pK=-log10(float(concentration) * UNIT_M[unit]),
                runs=len(runs[key]),
                receptors=len({r for r, _ in runs[key]}),
                native=scores["native"][key],
                bonded=scores["bonded"][key],
            )
        )
    return rows


def load_panels():
    grouped = defaultdict(list)
    for r in read(SOURCE / "results/affinity_poses.csv"):
        grouped[r["target"]].append({**r, **{a: float(r[a]) for a in ARMS}})
    panels = [
        ligand_rows(t, "archived", "fixed_native_top1", rows)
        for t, rows in sorted(grouped.items())
    ]

    # The frozen protocol is the original ITC label table, with explicit nM units.
    text = (SOURCE / "transfer/PROTOCOL.md").read_text()
    labels = dict(re.findall(r"^\| (6I\d\d) \| [^|]+ \| ([\d.]+) ±", text, re.M))
    del labels["6I66"]  # Unresolved stereochemical identity; primary panel excludes it.
    cells = defaultdict(dict)
    for r in read(SOURCE / "transfer/results/cells.csv"):
        if (
            r["engine"] == "pro"
            and r["receiver"] == "6I63"
            and r["query"] in labels
            and r["arm"] in ARMS
        ):
            key = (int(r["cap"]), r["query"], r["seed"])
            if r["arm"] in cells[key]:
                raise ValueError(f"Duplicate selected pose for {key}")
            cells[key][r["arm"]] = float(r["score"])
    by_budget = defaultdict(list)
    for (cap, ligand, seed), arms in sorted(cells.items()):
        by_budget[cap].append(
            dict(
                ligand=ligand,
                seed=seed,
                receptor="6I63",
                measurement=labels[ligand],
                **arms,
            )
        )
    for cap, rows in sorted(by_budget.items()):
        if {r["ligand"] for r in rows} != set(labels) or len(labels) != 6:
            raise ValueError("Incomplete registered ERRg panel")
        panels.append(ligand_rows("ERRg", cap, "arm_selected_top1", rows))
    return panels


def read_ligand_table(path):
    records = read(path)
    required = {
        "target",
        "assay",
        "endpoint",
        "ligand",
        "K",
        "unit",
        "baseline_score",
        "candidate_score",
    }
    if not records or required - records[0].keys():
        raise ValueError(
            "Empty CSV or missing required columns: " + ", ".join(sorted(required))
        )
    panels = defaultdict(dict)
    for number, r in enumerate(records, start=2):
        try:
            target, assay, endpoint, ligand, unit = (
                r[key].strip()
                for key in ("target", "assay", "endpoint", "ligand", "unit")
            )
            if not (
                target
                and assay
                and ligand
                and endpoint in {"Ki", "Kd"}
                and unit in UNIT_M
            ):
                raise ValueError(
                    "need target, assay, ligand, endpoint Ki/Kd and a supported concentration unit"
                )
            k = Decimal(r["K"])
            if not k.is_finite() or k <= 0:
                raise ValueError("K must be an exact positive finite concentration")
            potency = -float((k * Decimal(str(UNIT_M[unit]))).log10())
            base, candidate = float(r["baseline_score"]), float(r["candidate_score"])
            if not np.isfinite([potency, base, candidate]).all():
                raise ValueError("potency and scores must be finite")
            budget, rule = r.get("budget", "").strip(), r.get("pose_rule", "").strip()
            panel = panels[target, assay, endpoint, budget, rule]
            if ligand in panel:
                raise ValueError(
                    f"duplicate ligand {ligand!r} in one assay and condition"
                )
            panel[ligand] = dict(
                target=target,
                assay=assay,
                budget=budget,
                pose_rule=rule,
                ligand=ligand,
                endpoint=endpoint,
                K=k,
                unit=unit,
                pK=potency,
                baseline=base,
                candidate=candidate,
            )
        except (ArithmeticError, ValueError, AttributeError) as exc:
            raise ValueError(f"{path}, row {number}: {exc}") from exc
    return [[panel[k] for k in sorted(panel)] for _, panel in sorted(panels.items())]


def rank_correlations(x, y):
    if len(x) < 2 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return np.nan, np.nan
    return spearmanr(x, y).statistic, kendalltau(x, y, variant="b").statistic


def accuracy(credits):
    return credits.mean() if credits.size else np.nan


def bootstrap_rho(potency, predicted, indices):
    """Ranks are recomputed in every ligand resample, assigning average ties."""

    def centered_ranks(x):
        r = rankdata(x[indices], axis=1)
        return r - r.mean(axis=1, keepdims=True)

    x = centered_ranks(potency)
    y = centered_ranks(predicted)
    denom = np.linalg.norm(x, axis=1) * np.linalg.norm(y, axis=1)
    return np.divide(
        (x * y).sum(axis=1), denom, out=np.full(denom.shape, np.nan), where=denom > 0
    )


def intervals(values):
    valid = values[np.isfinite(values)]
    if not len(valid):
        return np.nan, np.nan, 0
    return (*np.percentile(valid, [2.5, 97.5]), len(valid))


def evaluate(rows, rng, replicates, arms=ARMS):
    first = rows[0]
    base = {k: first[k] for k in ("target", "assay", "budget", "pose_rule", "endpoint")}
    K = [r["K"] * Decimal(str(UNIT_M[r["unit"]])) for r in rows]
    pK = np.array([r["pK"] for r in rows])
    n = len(K)
    # A two-ligand nonconstant draw only recovers its original rank correlation.
    replicates = replicates if n >= 3 else 0
    indices = rng.integers(n, size=(replicates, n))
    pairs = [
        (i, j) if K[i] < K[j] else (j, i)
        for i, j in combinations(range(n), 2)
        if K[i] != K[j]
    ]
    strong, weak = np.asarray(pairs, dtype=int).reshape(-1, 2).T
    # Decimal comparison includes ratios exactly equal to three in the source units.
    threefold = np.array([K[j] >= 3 * K[i] for i, j in pairs], dtype=bool)
    scores = {}
    boot = {}
    credits = {}
    summaries = []
    for arm in arms:
        predicted = -np.asarray([r[arm] for r in rows])
        scores[arm] = predicted
        boot[arm] = (
            bootstrap_rho(pK, predicted, indices) if replicates else np.array([])
        )
        diff = predicted[strong] - predicted[weak]
        credits[arm] = (np.sign(diff) + 1) / 2
        lo, hi, valid = intervals(boot[arm])
        rho, tau = rank_correlations(pK, predicted)
        summaries.append(
            dict(
                **base,
                arm=arm,
                ligands=n,
                pK_min=min(pK),
                pK_max=max(pK),
                spearman=rho,
                kendall_tau_b=tau,
                pairs_total=n * (n - 1) // 2,
                pairs_experimental_ties=n * (n - 1) // 2 - len(pairs),
                pairs_unequal=len(pairs),
                prediction_ties=int(np.sum(diff == 0)),
                correct_all=credits[arm].sum(),
                accuracy_all=accuracy(credits[arm]),
                pairs_ge3fold=threefold.sum(),
                correct_ge3fold=credits[arm][threefold].sum(),
                accuracy_ge3fold=accuracy(credits[arm][threefold]),
                rho_boot_p025=lo,
                rho_boot_p975=hi,
                rho_boot_valid=valid,
                bootstrap_draws=replicates,
            )
        )
    baseline, candidate = arms
    lo, hi, valid = intervals(boot[candidate] - boot[baseline])
    all_gains = credits[candidate] - credits[baseline]
    # The whole ligand (and every pair containing it) is removed from both arms.
    loo_rho, loo_ge3fold = [], []
    for k in range(n):
        keep = np.arange(n) != k
        rhos = [rank_correlations(pK[keep], scores[a][keep])[0] for a in arms]
        if np.all(np.isfinite(rhos)):
            loo_rho.append(rhos[1] - rhos[0])
        subset = threefold & (strong != k) & (weak != k)
        if subset.any():
            loo_ge3fold.append(all_gains[subset].mean())
    paired = dict(
        **base,
        ligands=n,
        delta_spearman=summaries[1]["spearman"] - summaries[0]["spearman"],
        delta_tau_b=summaries[1]["kendall_tau_b"] - summaries[0]["kendall_tau_b"],
        delta_accuracy_all=accuracy(all_gains),
        delta_accuracy_ge3fold=accuracy(all_gains[threefold]),
        gain_credit_all=np.maximum(all_gains, 0).sum(),
        lost_credit_all=np.maximum(-all_gains, 0).sum(),
        gain_credit_ge3fold=np.maximum(all_gains[threefold], 0).sum(),
        lost_credit_ge3fold=np.maximum(-all_gains[threefold], 0).sum(),
        delta_rho_boot_p025=lo,
        delta_rho_boot_p975=hi,
        delta_rho_boot_valid=valid,
        bootstrap_draws=replicates,
        delta_rho_leave1_min=min(loo_rho, default=np.nan),
        delta_rho_leave1_max=max(loo_rho, default=np.nan),
        delta_accuracy_ge3fold_leave1_min=min(loo_ge3fold, default=np.nan),
        delta_accuracy_ge3fold_leave1_max=max(loo_ge3fold, default=np.nan),
    )
    return summaries, paired


def verify(summaries):
    original = {
        (r["target"], "archived", r["arm"]): r
        for r in read(SOURCE / "results/affinity_summary.csv")
        if r["receptor"] == "mean"
    }
    for r in read(SOURCE / "transfer/results/affinity.csv"):
        if r["engine"] == "pro" and r["population"] == "single_microstate":
            original["ERRg", int(r["cap"]), r["arm"]] = r
    for s in summaries:
        old = original[s["target"], s["budget"], s["arm"]]
        prior = [
            old["spearman"],
            old["kendall"],
            old["correct"],
            old.get("edges", old.get("eligible")),
        ]
        new = [
            s["spearman"],
            s["kendall_tau_b"],
            s["correct_ge3fold"],
            s["pairs_ge3fold"],
        ]
        np.testing.assert_allclose(
            new, np.array(prior, dtype=float), rtol=0, atol=1e-12
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "output",
        nargs="?",
        type=Path,
        help="Output directory; defaults to results/ beside this script for archived inputs, or affinity_results/ beside a custom input.",
    )
    parser.add_argument(
        "--ligands",
        type=Path,
        metavar="CSV",
        help="Evaluate a new CSV instead of the frozen panels. Required columns: target, assay, endpoint (Ki or Kd), "
        "ligand, K (exact positive concentration), unit (M/mM/uM/µM/μM/nM/pM), baseline_score, candidate_score "
        "(both lower-is-stronger). Optional budget and pose_rule separate conditions. Supply one row per ligand per "
        "assay and condition after combining technical repeats. Undefined metrics are written as nan. Bootstrap limits are descriptive and omitted for fewer than three ligands.",
    )
    args = parser.parse_args()
    output = args.output
    if output is None:
        output = (
            args.ligands.parent / "affinity_results"
            if args.ligands
            else Path(__file__).resolve().parent / "results"
        )
    try:
        panels = read_ligand_table(args.ligands) if args.ligands else load_panels()
    except ValueError as exc:
        parser.error(str(exc))
    arms = ("baseline", "candidate") if args.ligands else ARMS
    rng = np.random.default_rng(43103)
    summaries, paired = [], []
    for rows in panels:
        result, contrast = evaluate(rows, rng, replicates=10000, arms=arms)
        summaries.extend(result)
        paired.append(contrast)
    if not args.ligands:
        verify(summaries)
    output.mkdir(parents=True, exist_ok=True)
    save(output / "affinity_ligands.csv", [r for panel in panels for r in panel])
    save(output / "affinity_summary.csv", summaries)
    save(output / "affinity_paired.csv", paired)


if __name__ == "__main__":
    main()

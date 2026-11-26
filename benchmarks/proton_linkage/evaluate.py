"""Evaluate matched fluorosulfonamide docks, CAIX pH transfer and buffer calorimetry."""

import argparse
import csv
import itertools
import math
import statistics as st
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "tools"))
from ionization import GAS_CONSTANT_KCAL_MOL_K as R
from ionization import (
    anion_pair_offset,
    anion_population_penalty,
    log_kd_ratio,
    net_proton_uptake,
)

RECEPTORS = ("1V9E", "6SKT")
SEEDS = (314159, 271828)
ARMS = ("native", "point")


def read(path):
    with path.open() as source:
        return list(csv.DictReader(source))


def write(path, rows):
    with path.open("w", newline="") as target:
        writer = csv.DictWriter(target, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {k: round(v, 8) if isinstance(v, float) else v for k, v in row.items()}
            )


def raw_scores(runs, ligands):
    rows = []
    for code, seed, arm, ligand in itertools.product(RECEPTORS, SEEDS, ARMS, ligands):
        path = runs / f"{code}_{ligand}_{seed}_{arm}.pdbqt"
        models = [
            line
            for line in path.read_text().splitlines()
            if line.startswith("REMARK VINA RESULT:")
        ]
        rows.append(
            {
                "receptor": code,
                "seed": seed,
                "arm": arm,
                "ligand": ligand,
                "top1_kcal_mol": float(models[0].split()[3]),
                "modes": len(models),
            }
        )
    return rows


def fluorinated(scores, ligands):
    indices = {
        (s["receptor"], int(s["seed"]), s["arm"], s["ligand"]): float(
            s["top1_kcal_mol"]
        )
        for s in scores
    }
    pairs = []
    for first, second in itertools.combinations(ligands, 2):
        strong, weak = sorted(
            (first, second), key=lambda lig: float(lig["measured_kd_nm"])
        )
        ks, kw = [float(lig["measured_kd_nm"]) for lig in (strong, weak)]
        if kw < 3 * ks:
            continue
        ps, pw = [float(lig["pka"]) for lig in (strong, weak)]
        ds, dw = [float(lig["kd_deviation_nm"]) for lig in (strong, weak)]
        # All four paired differences are taken within the same arm, receptor and seed.
        for arm in ARMS:
            gaps = [
                indices[code, seed, arm, weak["ligand"]]
                - indices[code, seed, arm, strong["ligand"]]
                for code, seed in itertools.product(RECEPTORS, SEEDS)
            ]
            raw = st.median(gaps)
            offset = anion_pair_offset(pw, ps, 7.5, 298.15)
            pairs.append(
                {
                    "strong": strong["ligand"],
                    "weak": weak["ligand"],
                    "arm": arm,
                    "same_fluorine_count": strong["fluorines"] == weak["fluorines"],
                    "label_robust_3fold": kw - dw >= 3 * (ks + ds),
                    "observed": R * 298.15 * math.log(kw / ks),
                    "native": raw,
                    "min_native": min(gaps),
                    "max_native": max(gaps),
                    "population_only": offset,
                    "corrected": raw + offset,
                    "zinc_surrogate": raw + offset - 0.40 * (pw - ps),
                    "corrected_lo": raw
                    + anion_pair_offset(pw - 0.2, ps + 0.2, 7.5, 298.15),
                    "corrected_hi": raw
                    + anion_pair_offset(pw + 0.2, ps - 0.2, 7.5, 298.15),
                }
            )
    return pairs


def ph_transfer():
    ligands = read(HERE / "caix.csv")
    shifts = {
        d["ligand"]: math.log10(float(d["kd_ph5_nm"]) / float(d["kd_ph7_nm"]))
        for d in ligands
    }
    rows = []
    for d in ligands:
        compound = d["ligand"]
        pka = float(d["ligand_pka"])
        observed = shifts[compound]
        # Receptor pKa 7.01 comes from an external catalytic-domain kinetic titration.
        predicted = log_kd_ratio(pka, 7.01, 7, 5, 310.15) / math.log(10)
        ligand_only = (
            anion_population_penalty(pka, 5, 310.15)
            - anion_population_penalty(pka, 7, 310.15)
        ) / (R * 310.15 * math.log(10))
        loo = st.median(v for name, v in shifts.items() if name != compound)
        rows.append(
            {
                "ligand": compound,
                "observed_log10_shift": observed,
                "linked_log10_shift": predicted,
                "observed_kd5_nm": float(d["kd_ph5_nm"]),
                "predicted_kd5_nm": float(d["kd_ph7_nm"]) * 10**predicted,
                "ligand_only_log10_shift": ligand_only,
                "flat_log10_shift": 0.0,
                "other_ligands_median_log10_shift": loo,
                "literature_6p6_log10_shift": log_kd_ratio(pka, 6.6, 7, 5, 310.15)
                / math.log(10),
                "pg_domain_log10_shift": log_kd_ratio(pka, 6.49, 7, 5, 310.15)
                / math.log(10),
            }
        )
    return rows


def alkaline_transfer():
    # Fig S25 panel C, its caption and the methods give three different pHs.
    # Carry all reported values through; the high point in panel D is unambiguous.
    inputs = read(HERE / "caix_itc.csv")
    out = []
    for r in inputs:
        old, new = float(r["reference_ph"]), float(r["test_ph"])
        kd0, kd1 = float(r["reference_kd_nm"]), float(r["test_kd_nm"])
        linked = kd0 * math.exp(log_kd_ratio(8.29, 7.01, old, new, 310.15))
        # The receptor fraction can decrease by at most 10**(new-old).
        f_ligand_ratio = (1 + 10 ** (8.29 - new)) / (1 + 10 ** (8.29 - old))
        worst_ratio = ""
        if r["reference_fit_error_nm"]:
            worst_ratio = (kd1 - float(r["test_fit_error_nm"])) / (
                kd0 + float(r["reference_fit_error_nm"])
            )
        out.append(
            {
                "reference": r["reference"],
                "reference_ph": old,
                "measured_nm": kd1,
                "predicted_nm": linked,
                "observed_fold": kd1 / kd0,
                "predicted_fold": linked / kd0,
                "max_fold_any_pkas": 10 ** (new - old),
                "max_fold_measured_ligand_pka": 10 ** (new - old) * f_ligand_ratio,
                "low_fold_one_reported_error": worst_ratio,
            }
        )
    return out


def calorimetry():
    rows = []
    for r in read(HERE / "buffer_enthalpy.csv"):
        ph = float(r["ph"])
        pka = float(r["receptor_pka"])
        uptake = (float(r["tris_dh_kj_mol"]) - float(r["phosphate_dh_kj_mol"])) / (
            47.5 - 5.1
        )
        rows.append(
            {
                "receptor": r["receptor"],
                "ph": ph,
                "observed_uptake": uptake,
                "predicted_uptake": net_proton_uptake(8, pka, ph, 298.15),
            }
        )
    return rows


def summarize(pairs, transfers, enthalpies):
    for arm in ARMS:
        arm_rows = [r for r in pairs if r["arm"] == arm]
        for label, subset in [
            ("all", arm_rows),
            ("same-F", [r for r in arm_rows if r["same_fluorine_count"]]),
            ("label-robust", [r for r in arm_rows if r["label_robust_3fold"]]),
        ]:
            for field in ("native", "corrected", "population_only", "zinc_surrogate"):
                print(
                    f"{arm} {label} {field}: {sum(r[field] > 0 for r in subset)}/{len(subset)}, "
                    f"MAE {st.mean(abs(r[field] - r['observed']) for r in subset):.4f} kcal/mol"
                )
        print(
            f"Independent ±0.2 pKa ranges: {sum(r['corrected_lo'] > 0 for r in arm_rows)} always right, "
            f"{sum(r['corrected_hi'] < 0 for r in arm_rows)} always wrong."
        )
    for method in (
        "flat",
        "ligand_only",
        "linked",
        "pg_domain",
        "literature_6p6",
        "other_ligands_median",
    ):
        field = f"{method}_log10_shift"
        err = [r[field] - r["observed_log10_shift"] for r in transfers]
        print(
            f"CAIX {method}: mean absolute {st.mean(map(abs, err)):.4f} log10, "
            f"max fold error {10 ** max(map(abs, err)):.3f}"
        )
    for field in ("linked_log10_shift", "flat_log10_shift"):
        errors = [
            abs(
                (a[field] - b[field])
                - (a["observed_log10_shift"] - b["observed_log10_shift"])
            )
            for a, b in itertools.combinations(transfers, 2)
        ]
        print(f"CAIX differential {field}: {st.mean(errors):.4f} log10")
    err = [r["predicted_uptake"] - r["observed_uptake"] for r in enthalpies]
    print(
        f"Buffer-derived proton uptake: {sum(r['predicted_uptake'] * r['observed_uptake'] > 0 for r in enthalpies)}/"
        f"{len(err)} signs; MAE {st.mean(map(abs, err)):.4f} protons "
        f"(zero uptake {st.mean(abs(r['observed_uptake']) for r in enthalpies):.4f})"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--runs",
        type=Path,
        help="Folder with 56 docked PDBQTs. Otherwise use frozen top scores.",
    )
    parser.add_argument("--out", type=Path, default=HERE / "results")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    ligands = read(HERE / "fluorinated.csv")
    if args.runs:
        scores = raw_scores(args.runs, (d["ligand"] for d in ligands))
        write(args.out / "docking.csv", scores)
    else:
        scores = read(HERE / "results/docking.csv")
    pairs = fluorinated(scores, ligands)
    transfers = ph_transfer()
    enthalpies = calorimetry()
    write(args.out / "pairs.csv", pairs)
    write(args.out / "transfer.csv", transfers)
    write(args.out / "uptake.csv", enthalpies)
    write(args.out / "alkaline.csv", alkaline_transfer())
    summarize(pairs, transfers, enthalpies)


if __name__ == "__main__":
    main()

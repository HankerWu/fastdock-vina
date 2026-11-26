"""Check a published-result replay against direct mass action and raw docking output."""

import argparse
import csv
import itertools
import math
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "tools"))
from ionization import log_kd_ratio, net_proton_uptake


def read(name):
    with (HERE / name).open() as source:
        return list(csv.DictReader(source))


def weight(pka_l, pka_r, ph):
    return 1 / ((1 + 10 ** (pka_l - ph)) * (1 + 10 ** (ph - pka_r)))


def verify_mass_action():
    for pl, pr, ph in itertools.product(
        (5.8, 8.29, 10.1), (6.49, 7.01, 8.4), (4, 7, 10)
    ):
        ratio = log_kd_ratio(pl, pr, 7, ph, 310.15)
        assert abs(ratio - math.log(weight(pl, pr, 7) / weight(pl, pr, ph))) < 1e-12
        h = 1e-5
        diff = math.log10(weight(pl, pr, ph - h) / weight(pl, pr, ph + h)) / (2 * h)
        assert abs(net_proton_uptake(pl, pr, ph, 310.15) - diff) < 1e-9
        assert abs(diff) <= 1
        assert abs(ratio + log_kd_ratio(pl, pr, ph, 7, 310.15)) < 1e-12
    for temp in (0, -1, math.inf, math.nan):
        try:
            log_kd_ratio(7, 8, 6, 9, temp)
        except ValueError:
            continue
        raise AssertionError(temp)


def verify_experiments():
    acidic = read("caix.csv")
    empirical = []
    for row in acidic:
        pka = float(row["ligand_pka"])
        reference = float(row["kd_ph7_nm"])
        expected = reference * weight(pka, 7.01, 7) / weight(pka, 7.01, 5)
        empirical.append(abs(math.log10(expected / float(row["kd_ph5_nm"]))))
    assert round(statistics.mean(empirical), 4) == 0.0938
    limits = read("results/alkaline.csv")
    for row in limits:
        assert float(row["observed_fold"]) > float(row["max_fold_any_pkas"])
    for row in limits[1:]:
        assert float(row["low_fold_one_reported_error"]) > float(
            row["max_fold_any_pkas"]
        )


def verify_scores(runs):
    frozen = read("results/docking.csv")
    ligands = {r["ligand"]: r for r in read("fluorinated.csv")}
    expected = set(
        itertools.product(
            ("1V9E", "6SKT"), ("native", "point"), ("271828", "314159"), ligands
        )
    )
    scores = {}
    for row in frozen:
        key = (row["receptor"], row["arm"], row["seed"], row["ligand"])
        assert key not in scores
        scores[key] = float(row["top1_kcal_mol"])
        if runs:
            code, arm, seed, lig = key
            contents = (
                (runs / f"{code}_{lig}_{seed}_{arm}.pdbqt").read_text().splitlines()
            )
            assert contents[0].split() == ["MODEL", "1"]
            assert float(contents[1].split()[3]) == scores[key]
            assert sum(line == "ENDMDL" for line in contents) == int(row["modes"])
    assert set(scores) == expected
    for arm in ("native", "point"):
        all_pairs = 0
        raw_wins = 0
        corrected_wins = 0
        for weak, strong in itertools.permutations(ligands, 2):
            w, s = ligands[weak], ligands[strong]
            if float(w["measured_kd_nm"]) < 3 * float(s["measured_kd_nm"]):
                continue
            all_pairs += 1
            gap = statistics.median(
                scores[c, arm, seed, weak] - scores[c, arm, seed, strong]
                for c, seed in itertools.product(("1V9E", "6SKT"), ("271828", "314159"))
            )
            raw_wins += gap > 0
            # A common receptor partition cancels even when its pKa is changed.
            corr = (
                0.00198720425864083
                * 298.15
                * math.log(
                    weight(float(s["pka"]), 8.4, 7.5)
                    / weight(float(w["pka"]), 8.4, 7.5)
                )
            )
            corrected_wins += gap + corr > 0
        assert (all_pairs, raw_wins, corrected_wins) == (15, 14, 13)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", type=Path)
    args = parser.parse_args()
    verify_mass_action()
    verify_experiments()
    verify_scores(args.runs)
    print(
        "PASS: mass action, uptake derivatives, range bound, pH reversal, external observations and matched docking ranks"
    )

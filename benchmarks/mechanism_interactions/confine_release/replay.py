"""Rebuild the published T4L Val111 thermodynamic cycle from Figure 2."""

import csv
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
RT = 0.00198720425864083 * 300  # kcal/mol at the simulated temperature


def load_rows(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def reconstruct():
    source = load_rows(HERE / "source.csv")
    archive = load_rows(HERE.parents[1] / "affinity_partition" / "published.csv")
    historical = next(row for row in archive if row["pdb"] == "187L")
    experimental = float(historical["exp"])
    endpoints = {
        "apo_like": ("cm2_multi", "cm2_val"),
        "holo_like": ("holo_cm2_multi", "holo_cm2_val"),
    }
    result = {"experimental_kcal_mol": experimental, "states": {}}
    contributions = []
    confined_energies = {}
    for row in source:
        state = row["rotamer"]
        no_ligand = float(row["apo_confine"])
        with_ligand = float(row["bound_confine"])
        p_free = math.exp(-no_ligand / RT)
        p_bound = math.exp(-with_ligand / RT)
        values = {"p_free_from_figure": p_free, "p_bound_from_figure": p_bound}
        result["states"][state] = values
        if not row["confined_binding"]:
            continue
        confined = float(row["confined_binding"])
        cycle = no_ligand + confined - with_ligand
        old_confined, old_released = endpoints[state]
        # The related JMB paper shares calculations and rounds its endpoints.
        if confined != float(historical[old_confined]) or not math.isclose(
            cycle, float(historical[old_released]), rel_tol=0, abs_tol=0.011
        ):
            raise ValueError(f"Cycle disagrees with the archived endpoint for {state}")
        values.update(
            cycle_kcal_mol=cycle,
            signed_error_before=confined - experimental,
            signed_error_after=cycle - experimental,
        )
        contributions.append(p_free * math.exp(-confined / RT))
        confined_energies[state] = confined

    states = result["states"]
    result["raw_p_free_sum"] = sum(s["p_free_from_figure"] for s in states.values())
    result["raw_p_bound_sum"] = sum(s["p_bound_from_figure"] for s in states.values())
    result["confined_holo_minus_apo"] = (
        confined_energies["holo_like"] - confined_energies["apo_like"]
    )
    result["cycle_holo_minus_apo"] = (
        states["holo_like"]["cycle_kcal_mol"] - states["apo_like"]["cycle_kcal_mol"]
    )
    result["two_state_partial_kcal_mol"] = -RT * math.log(sum(contributions))
    result["two_state_equal_weights_kcal_mol"] = -RT * math.log(
        sum(math.exp(-energy / RT) for energy in confined_energies.values())
        / len(confined_energies)
    )
    return result


if __name__ == "__main__":
    print(json.dumps(reconstruct(), indent=2, sort_keys=True, allow_nan=False))

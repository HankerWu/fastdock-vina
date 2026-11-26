"""Solution-state population penalties for a monoacid binding as its anion."""

import math

GAS_CONSTANT_KCAL_MOL_K = 0.00198720425864083


def anion_population_penalty(pka, ph, kelvin):
    """Return −RT ln(f_anion) in kcal/mol for the supplied solution conditions."""
    if not all(math.isfinite(v) for v in (pka, ph, kelvin)) or kelvin <= 0:
        raise ValueError(
            "pKa and pH must be finite and temperature must be positive and finite"
        )
    exponent = math.log(10.0) * (pka - ph)
    # Stable evaluation of ln(1 + exp(exponent)) for rare and common anions.
    log_inverse_fraction = max(exponent, 0.0) + math.log1p(math.exp(-abs(exponent)))
    return GAS_CONSTANT_KCAL_MOL_K * kelvin * log_inverse_fraction


def anion_pair_offset(first_pka, second_pka, ph, kelvin):
    """Return the additive correction to G(first) − G(second) at one pH and T."""
    return anion_population_penalty(first_pka, ph, kelvin) - anion_population_penalty(
        second_pka, ph, kelvin
    )


def linkage_penalty(ligand_pka, receptor_pka, ph, kelvin):
    """Cost to select an anionic monoacid and the protonated receptor at this pH."""
    ligand = anion_population_penalty(ligand_pka, ph, kelvin)
    receptor = anion_population_penalty(-receptor_pka, -ph, kelvin)
    return ligand + receptor


def log_kd_ratio(ligand_pka, receptor_pka, reference_ph, new_ph, kelvin):
    """Predict ln[Kd(new pH)/Kd(reference pH)] at fixed intrinsic affinity and T.

    Only the ligand anion and protonated receptor lead to the unique bound state.
    Both pKas must apply at the common temperature and solution conditions.
    """
    old = linkage_penalty(ligand_pka, receptor_pka, reference_ph, kelvin)
    new = linkage_penalty(ligand_pka, receptor_pka, new_ph, kelvin)
    return (new - old) / (GAS_CONSTANT_KCAL_MOL_K * kelvin)


def net_proton_uptake(ligand_pka, receptor_pka, ph, kelvin):
    """Mean protons taken from solution on binding, equal to d log10(Kd)/d pH."""
    thermal = GAS_CONSTANT_KCAL_MOL_K * kelvin
    ligand_anion = math.exp(-anion_population_penalty(ligand_pka, ph, kelvin) / thermal)
    receptor_acid = math.exp(
        -anion_population_penalty(-receptor_pka, -ph, kelvin) / thermal
    )
    return ligand_anion - receptor_acid


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=log_kd_ratio.__doc__)
    parser.add_argument("--ligand-pka", type=float, required=True)
    parser.add_argument("--receptor-pka", type=float, required=True)
    parser.add_argument("--reference-ph", type=float, required=True)
    parser.add_argument("--new-ph", type=float, required=True)
    parser.add_argument("--kelvin", type=float, required=True)
    args = parser.parse_args()
    print(
        f"predicted log10 Kd ratio (new/reference) = {log_kd_ratio(**vars(args)) / math.log(10):.8g}"
    )

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

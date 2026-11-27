"""Closed-form regression for the two final FKBP AD4 rescore operators."""
import itertools
import unittest

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from analyze import calculate_components


def atom(charge):
    return (f"ATOM      1  C   UNL     1    {0.5:8.3f}{0.5:8.3f}{0.5:8.3f}"
            f"  1.00  0.00    {charge:+6.3f} C")


class MapOperators(unittest.TestCase):
    def test_two_caps_and_printed_charge(self):
        axes = [[0., 1.]] * 3

        def field(value):
            return RegularGridInterpolator(axes, np.full((2, 2, 2), value))

        def cap(value):
            return value * 1000 / (1000 + value) if value > 0 else value

        for a, e, d, q in itertools.product((-200., 200.), (-140., 140.),
                                              (-60., 60.), (-.5, 0., .5)):
            fields = {"C": field(a), "e": field(e), "d": field(d)}
            lines = [atom(q), "TORSDOF 3"]
            separate = cap(a) + q * cap(e) + abs(q) * cap(d)
            combined = cap(a + q * e + abs(q) * d)
            for name, expected in (("separate", separate), ("combined", combined)):
                for raw in (True, False):
                    terms = calculate_components(lines, fields, operator=name, uncapped_terms=raw)
                    self.assertAlmostEqual(sum(terms.values()), expected + 0.8949, places=12)
                    if raw:
                        self.assertAlmostEqual(terms["type_maps"], a)
                        self.assertAlmostEqual(terms["electrostatic"], q * e)
                        self.assertAlmostEqual(terms["charge_desolvation"], abs(q) * d)
                self.assertAlmostEqual(terms["torsion"], 0.8949)


if __name__ == "__main__":
    unittest.main()

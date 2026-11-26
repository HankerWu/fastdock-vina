"""Check monoacid population penalties independently of the docking outcomes."""

import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
from ionization import (
    GAS_CONSTANT_KCAL_MOL_K,
    anion_pair_offset,
    anion_population_penalty,
)


class PopulationTests(unittest.TestCase):
    def test_equal_and_limiting_populations(self):
        rt = GAS_CONSTANT_KCAL_MOL_K * 298.15
        self.assertAlmostEqual(anion_population_penalty(7, 7, 298.15), rt * math.log(2))
        self.assertEqual(anion_population_penalty(-400, 0, 298.15), 0.0)
        self.assertAlmostEqual(
            anion_population_penalty(400, 0, 298.15), rt * 400 * math.log(10)
        )
        self.assertEqual(anion_pair_offset(7.3, 7.3, 7, 298.15), 0)

    def test_measured_pair_and_sign(self):
        offset = anion_pair_offset(7.3, 7.1, 7.0, 298.15)
        self.assertAlmostEqual(offset, 0.16716524022058915, places=12)
        self.assertEqual(offset, -anion_pair_offset(7.1, 7.3, 7.0, 298.15))

    def test_invalid_conditions(self):
        for args in (
            (7, 7, 0),
            (7, 7, -1),
            (float("nan"), 7, 298),
            (7, float("inf"), 298),
            (7, 7, float("nan")),
        ):
            with self.subTest(args=args), self.assertRaises(ValueError):
                anion_population_penalty(*args)


if __name__ == "__main__":
    unittest.main()

"""Closure carbons and construction markers in the geometry evaluator."""

import unittest

from paired import atoms


class ClosureGeometryTest(unittest.TestCase):
    def test_physical_carbon_retained_for_every_closure_pair(self):
        types = ["C", "H", "HD", "CG0", "G0", "CG1", "G1", "CG2", "G2", "CG3", "G3"]
        lines = [
            f"ATOM  {i:5d}  C   UNL A   1       1.000   2.000   3.000  1.00  0.00     0.000 {typ}"
            for i, typ in enumerate(types, 1)
        ]
        parsed = atoms(lines)
        self.assertEqual(list(parsed), [1, 4, 6, 8, 10])
        self.assertTrue(all(value == parsed[1] for value in parsed.values()))


if __name__ == "__main__":
    unittest.main()

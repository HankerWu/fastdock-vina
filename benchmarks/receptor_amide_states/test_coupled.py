"""Test the binary partition and its native-coordinate coupling against direct sums."""

from itertools import product
import unittest

import numpy as np
from scipy.special import logsumexp

from coupled import Component, coupling
from model import RT, pair_energy


class JointPartitionTest(unittest.TestCase):
    def test_mixed_difference_from_four_explicit_contact_tables(self):
        left = np.array([[0.0, 0.0, 0.0], [1.3, 1.7, 0.0]])
        right = np.array([[3.1, 0.0, 0.0], [3.5, 2.2, 1.0]])

        def direct(a, b):
            energy = 0.0
            for i, j in product(range(2), repeat=2):
                u, v = i ^ a, j ^ b
                energy += pair_energy(
                    np.linalg.norm(left[i] - right[j]),
                    (1.7, 1.8)[u],
                    (1.7, 1.8)[v],
                    u != v,
                )
            return energy

        mixed = direct(1, 1) - direct(1, 0) - direct(0, 1) + direct(0, 0)
        self.assertAlmostEqual(coupling(left, right), mixed)
        self.assertAlmostEqual(coupling(left[::-1], right), -mixed)
        self.assertEqual(coupling(left, right + 20.0), 0.0)

    def test_partition_matches_enumeration_and_coordinate_relabelling(self):
        rng = np.random.default_rng(314159)
        for n in range(1, 9):
            ids = list(range(n))
            edges = {
                (i, j): rng.normal()
                for i in ids
                for j in range(i + 1, n)
                if (i + j) % 3
            }
            h, dl = rng.normal(size=(2, n))
            states = np.array(list(product((0.0, 1.0), repeat=n)))
            e = states @ h + sum(
                v * states[:, i] * states[:, j] for (i, j), v in edges.items()
            )
            graph = Component(ids, edges, h)
            expected = RT * (logsumexp(-e / RT) - logsumexp(-(e + states @ dl) / RT))
            result = RT * (graph.apo - graph.log_partition(h + dl))
            self.assertAlmostEqual(expected, result, places=12)
            self.assertAlmostEqual(graph.apo, logsumexp(-e / RT), places=12)
            self.assertEqual(graph.log_partition(h), graph.apo)
            # Redefine the label of vertex zero. Its background becomes a constant
            # energy; ligand scoring at the new reference also gains dl[0].
            hr, dnew, er = h.copy(), dl.copy(), edges.copy()
            hr[0], dnew[0] = -h[0], -dl[0]
            for (i, j), v in edges.items():
                if i == 0:
                    hr[j] += v
                    er[i, j] = -v
            new = Component(ids, er, hr)
            self.assertAlmostEqual(
                result, dl[0] + RT * (new.apo - new.log_partition(hr + dnew)), places=12
            )


if __name__ == "__main__":
    unittest.main()

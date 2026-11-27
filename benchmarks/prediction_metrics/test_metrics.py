"""Regression checks for the public pose and assay metric inputs."""

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np
import pandas as pd

import affinity
import pose


class AffinityMetricsTest(unittest.TestCase):
    def evaluate(self, records):
        with TemporaryDirectory() as folder:
            path = Path(folder) / "ligands.csv"
            path.write_text(
                "target,assay,endpoint,ligand,K,unit,baseline_score,candidate_score\n"
                + records,
                encoding="utf-8",
            )
            (rows,) = affinity.read_ligand_table(path)
        summaries, _ = affinity.evaluate(
            rows, np.random.default_rng(1), 10, arms=("baseline", "candidate")
        )
        return rows, summaries

    def test_units_score_orientation_ties_and_inclusive_fold_cutoff(self):
        for endpoint in ("Ki", "Kd"):
            with self.subTest(endpoint=endpoint):
                records = "".join(
                    f"T,A,{endpoint},{ligand}\n"
                    for ligand in (
                        "a,1,nM,-8,-7",
                        "b,0.001,µM,-8,-7",
                        "c,3,nM,-6,-7",
                        "d,1.5,nM,-7,-7",
                    )
                )
                rows, (baseline, candidate) = self.evaluate(records)
                self.assertEqual([r["pK"] for r in rows[:2]], [9.0, 9.0])
                self.assertAlmostEqual(baseline["spearman"], 1)
                self.assertAlmostEqual(baseline["kendall_tau_b"], 1)
                self.assertEqual(baseline["accuracy_all"], 1)
                self.assertEqual(baseline["pairs_experimental_ties"], 1)
                self.assertEqual(baseline["pairs_unequal"], 5)
                self.assertEqual(baseline["pairs_ge3fold"], 2)
                self.assertEqual(candidate["prediction_ties"], 5)
                self.assertEqual(candidate["accuracy_ge3fold"], 0.5)
                self.assertTrue(np.isnan(candidate["spearman"]))

    def test_single_ligand_has_undefined_ranking(self):
        _, summaries = self.evaluate("T,A,Kd,a,1,nM,-7,-8\n")
        for result in summaries:
            self.assertEqual(result["pairs_total"], 0)
            self.assertEqual(result["bootstrap_draws"], 0)
            for metric in (
                "spearman",
                "kendall_tau_b",
                "accuracy_all",
                "accuracy_ge3fold",
            ):
                self.assertTrue(np.isnan(result[metric]), metric)


class PoseMetricsTest(unittest.TestCase):
    def test_equal_mass_per_target_query_receptor_then_seed(self):
        frame = pd.DataFrame(
            [
                ("A", "q", "r", 1),
                ("B", "p", "r0", 1),
                ("B", "p", "r1", 1),
                ("B", "p", "r1", 2),
                ("B", "t", "r", 1),
            ],
            columns=["target", "query", "receiver", "seed"],
        )
        weights = pose.balanced_weights(frame)
        np.testing.assert_allclose(weights, [1 / 2, 1 / 8, 1 / 16, 1 / 16, 1 / 4])
        self.assertEqual(pose.quantiles([0, 10, 10, 10, 10], weights, [0.5]), [5])

    def test_reranking_cannot_change_oracle_in_a_fixed_candidate_pool(self):
        baseline = dict(
            study="redock",
            task="pose",
            budget="full",
            cap=9,
            target="T",
            query="q",
            receiver="r",
            seed=1,
            pool="fixed",
            method="base",
            reference_pool="fixed",
            reference_method="base",
            n_poses=9,
            top1_A=3.0,
            best5_A=2.0,
            oracle_A=1.0,
        )
        changed = dict(baseline, method="candidate", oracle_A=1.5)
        with self.assertRaisesRegex(ValueError, "Pool coverage changed"):
            pose.pair_runs(pd.DataFrame([baseline, changed]))


if __name__ == "__main__":
    unittest.main()

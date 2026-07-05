#!/usr/bin/env python3
import importlib.util
import pathlib
import unittest

MODULE_PATH = pathlib.Path(__file__).resolve().parents[1] / "pipeline" / "scripts" / "rna" / "build_lamar_training_table.py"
spec = importlib.util.spec_from_file_location("build_lamar_training_table", MODULE_PATH)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class TrainingLabelStatisticsTests(unittest.TestCase):
    def test_median_abs_deviation(self):
        self.assertAlmostEqual(module.median_abs_deviation([0.1, 0.2, 0.3]), 0.1)

    def test_bh_fdr_is_monotone_in_rank(self):
        p_values = [0.01, 0.04, 0.03, 0.002]
        adjusted = module.bh_fdr(p_values)
        ranked = sorted(zip(p_values, adjusted))
        self.assertTrue(all(ranked[i][1] <= ranked[i + 1][1] for i in range(len(ranked) - 1)))
        self.assertTrue(all(0 <= value <= 1 for value in adjusted))

    def test_fisher_exact_detects_large_difference(self):
        self.assertLess(module.fisher_exact_two_sided(30, 70, 1, 99), 1e-5)

    def test_fisher_exact_no_difference(self):
        self.assertAlmostEqual(module.fisher_exact_two_sided(10, 90, 10, 90), 1.0)

    def test_reverse_complement(self):
        self.assertEqual(module.reverse_complement("ACGTN"), "NACGT")


if __name__ == "__main__":
    unittest.main()

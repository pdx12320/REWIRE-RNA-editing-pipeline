#!/usr/bin/env python3
import math
import sys
import unittest
from pathlib import Path

sys.path.insert(
    0,
    str(Path(__file__).resolve().parents[1] / "pipeline" / "scripts" / "rna"),
)

from run_audited_lamar_background_correction import (
    build_labels,
    benjamini_hochberg,
    cluster_positions,
    exact_fisher_two_sided,
    mad,
    parse_mpileup_bases,
    revcomp,
    validate_outputs,
)


class PipelineUnitTests(unittest.TestCase):
    def test_reverse_complement(self):
        self.assertEqual(revcomp("AACGTN"), "NACGTT")

    def test_mad(self):
        self.assertAlmostEqual(mad([0.1, 0.2, 0.3]), 0.1)
        self.assertIsNone(mad([]))

    def test_bh_is_monotonic_by_pvalue(self):
        adjusted = benjamini_hochberg([0.01, 0.04, None, 0.03])
        self.assertIsNone(adjusted[2])
        self.assertAlmostEqual(adjusted[0], 0.03)
        self.assertAlmostEqual(adjusted[1], 0.04)
        self.assertAlmostEqual(adjusted[3], 0.04)

    def test_fisher_exact(self):
        value = exact_fisher_two_sided(1, 9, 11, 3)
        self.assertTrue(math.isclose(value, 0.002759456, rel_tol=1e-5))

    def test_mpileup_parser_strand_and_indels(self):
        counts = parse_mpileup_bases(".,Tt^F.+2ag,-1cA$", "C")
        self.assertEqual(counts["C"], 4)
        self.assertEqual(counts["C_forward"], 2)
        self.assertEqual(counts["C_reverse"], 2)
        self.assertEqual(counts["T"], 2)
        self.assertEqual(counts["A"], 1)

    def test_position_clustering(self):
        clusters = cluster_positions([1, 3, 500, 501], maximum_gap=10)
        self.assertEqual([(a, b) for a, b, _ in clusters], [(1, 4), (500, 502)])

    def test_insufficient_control_keeps_descriptive_median_but_not_label(self):
        candidate = {"chrom": "chr1", "position": "101", "ref": "C", "alt": "T", "vep_strand": "1"}
        rows = []
        samples = (
            ("CU517_GC_T1", "treated", 1, 20, 4),
            ("CU517_GC_T2", "treated", 2, 20, 6),
            ("CU517_GC_T3", "treated", 3, 20, 8),
            ("CU517_GC_C1", "control", 1, 20, 2),
            ("CU517_GC_C2", "control", 2, 10, 0),
            ("CU517_GC_C3", "control", 3, 0, 0),
        )
        for sample, group, replicate, depth, alt in samples:
            rows.append(
                {
                    **candidate,
                    "transcript_strand": "1",
                    "sample": sample,
                    "group": group,
                    "replicate": replicate,
                    "usable_depth": depth,
                    "ref_count": depth - alt,
                    "alt_count": alt,
                    "edit_rate": alt / depth if depth else None,
                    "coverage_status": "adequate_with_alt" if depth >= 20 else "inadequate_low_depth",
                }
            )
        context = {
            ("chr1", 101, "C", "T"): {
                "orientation_qc": "pass",
            }
        }
        labels = build_labels([candidate], rows, 20, 2, context)
        self.assertEqual(labels[0]["control_median"], 0.1)
        self.assertIsNone(labels[0]["raw_edit_rate_difference"])
        self.assertIsNone(labels[0]["corrected_editing_efficiency"])
        self.assertFalse(labels[0]["training_eligible"])
        self.assertIn("low_control_coverage", labels[0]["exclusion_reason"])
        validate_outputs(labels, 1)


if __name__ == "__main__":
    unittest.main()

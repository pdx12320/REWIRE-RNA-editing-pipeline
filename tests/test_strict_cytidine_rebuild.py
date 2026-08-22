#!/usr/bin/env python3
import tempfile
import unittest
from pathlib import Path
import sys


RNA_SCRIPTS = Path(__file__).resolve().parents[1] / "pipeline" / "scripts" / "rna"
sys.path.insert(0, str(RNA_SCRIPTS))

import build_augmented_reference as augmented
import build_full_coverage_cytidine_labels as strict
import finalize_strict_cytidine_dataset as finalize
import pileup_candidate_base_counts as candidate_counts
import audit_sequence_mappability as mappability


class FakeRead:
    def __init__(self, base="C", nh=None):
        self.is_unmapped = False
        self.is_duplicate = False
        self.is_secondary = False
        self.is_supplementary = False
        self.is_qcfail = False
        self.is_reverse = False
        self.mapping_quality = 60
        self.query_qualities = [40]
        self.query_sequence = base
        self._nh = nh

    def has_tag(self, tag):
        return tag == "NH" and self._nh is not None

    def get_tag(self, tag):
        if not self.has_tag(tag):
            raise KeyError(tag)
        return self._nh


class FakePileupRead:
    def __init__(self, read):
        self.alignment = read
        self.is_del = False
        self.is_refskip = False
        self.query_position = 0


class FakeColumn:
    reference_pos = 9

    def __init__(self, reads):
        self.pileups = [FakePileupRead(read) for read in reads]


class FakeAlignment:
    def __init__(self, reads):
        self.column = FakeColumn(reads)

    def pileup(self, *args, **kwargs):
        return [self.column]


class StrictCytidineRebuildTests(unittest.TestCase):
    def test_mappability_requires_full_length_alignment(self):
        self.assertTrue(mappability.full_length_alignment("101M", 101))
        self.assertTrue(mappability.full_length_alignment("50M1I50M", 101))
        self.assertFalse(mappability.full_length_alignment("5S96M", 101))

    def test_candidate_count_requires_nh_equal_one(self):
        alignment = FakeAlignment([FakeRead("C", 1), FakeRead("T", 2), FakeRead("T", None)])
        result = candidate_counts.count_site(alignment, "chr1", 10, 30, 20, 1000)
        self.assertEqual(result["C_count"], 1)
        self.assertEqual(result["T_count"], 0)
        self.assertEqual(result["excluded_nh_multimapper"], 1)
        self.assertEqual(result["excluded_missing_nh"], 1)

    def test_reporter_context_is_verified_before_append(self):
        with tempfile.TemporaryDirectory() as temporary:
            temporary = Path(temporary)
            genome = temporary / "genome.fa"
            reporter = temporary / "reporter.fa"
            output = temporary / "augmented.fa"
            genome.write_text(">chr1\nACGTACGT\n", encoding="utf-8")
            sequence = list("A" * 500)
            sequence[457] = "G"
            sequence[458] = "C"
            reporter.write_text(">verified_construct\n{}\n".format("".join(sequence)), encoding="utf-8")
            report = augmented.build_reference(
                genome, reporter, output, "EGFP_GC_reporter", 459, "GC"
            )
            self.assertEqual(report["target_context"], "GC")
            records = list(augmented.read_fasta(output))
            self.assertEqual([name for name, _ in records], ["chr1", "EGFP_GC_reporter"])

    def test_wrong_reporter_context_stops(self):
        with self.assertRaisesRegex(ValueError, "validation failed"):
            augmented.validate_reporter("A" * 500, 459, "GC")

    def test_positive_and_strict_negative_rules(self):
        manifest = [
            {"sample": "T1", "group": "treated"},
            {"sample": "T2", "group": "treated"},
            {"sample": "T3", "group": "treated"},
            {"sample": "C1", "group": "control"},
            {"sample": "C2", "group": "control"},
            {"sample": "C3", "group": "control"},
        ]
        positive = strict.Site(100, "C", "T", "+")
        rates = {"T1": (90, 10), "T2": (80, 20), "T3": (85, 15), "C1": (100, 0), "C2": (100, 0), "C3": (100, 0)}
        for sample, (ref, alt) in rates.items():
            positive.sample_counts[sample] = {
                "ref_count": ref, "alt_count": alt, "allele_depth": ref + alt,
                "acgt_depth": ref + alt, "edit_rate": alt / (ref + alt),
            }
        label = strict.classify_site(positive, manifest, 0.10, 2, 3, 0.01)
        self.assertEqual(label["label_class"], "positive")

        negative = strict.Site(200, "C", "T", "+")
        for row in manifest:
            negative.sample_counts[row["sample"]] = {
                "ref_count": 100, "alt_count": 0, "allele_depth": 100,
                "acgt_depth": 100, "edit_rate": 0.0,
            }
        label = strict.classify_site(negative, manifest, 0.10, 2, 3, 0.01)
        self.assertEqual(label["label_class"], "strict_negative")

    def test_deterministic_one_to_two_selection(self):
        rows = []
        for index in range(2):
            rows.append({
                "chrom": "chr1", "position": str(index + 1), "ref": "C", "alt": "T",
                "training_eligible": "1", "mappability_pass": "1", "label_class": "positive",
            })
        for index in range(10):
            rows.append({
                "chrom": "chr2", "position": str(index + 1), "ref": "C", "alt": "T",
                "training_eligible": "1", "mappability_pass": "1", "label_class": "strict_negative",
            })
        selected, audit = finalize.select_rows(rows, negative_ratio=2, seed=7)
        self.assertEqual(audit["selected_positive"], 2)
        self.assertEqual(audit["selected_strict_negative"], 4)
        self.assertEqual(len(selected), 6)


if __name__ == "__main__":
    unittest.main()

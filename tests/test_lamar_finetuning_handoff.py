#!/usr/bin/env python3
import copy
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


RNA_SCRIPTS = Path(__file__).resolve().parents[1] / "pipeline" / "scripts" / "rna"
sys.path.insert(0, str(RNA_SCRIPTS))

import export_lamar_scalar_regression as exporter
import prepare_lamar_finetuning_handoff as handoff

BASELINE_PATH = Path(__file__).resolve().parents[1] / "examples" / "train_scalar_baseline.py"
BASELINE_SPEC = importlib.util.spec_from_file_location("train_scalar_baseline", BASELINE_PATH)
baseline = importlib.util.module_from_spec(BASELINE_SPEC)
BASELINE_SPEC.loader.exec_module(baseline)
SKLEARN_AVAILABLE = importlib.util.find_spec("sklearn") is not None


def unique_sequence(index):
    sequence = list(("ACGT" * 26)[:101])
    value = index
    for offset in range(8):
        sequence[offset] = "ACGT"[value % 4]
        value //= 4
    sequence[50] = "C"
    return "".join(sequence)


def synthetic_rows(count=18):
    labels = []
    metadata = []
    for index in range(count):
        position = 1000 + index * 300
        if index == 2:
            position = 2000
        elif index == 3:
            position = 2050
        sequence = unique_sequence(0 if index == 1 else index)
        eligible = index < count - 2
        confidence = "moderate" if eligible and index in {4, 5} else "high" if eligible else "low"
        if not eligible:
            treated = "0.05"
            control = "NA"
            raw = "NA"
            corrected = "NA"
            reason = "low_control_coverage"
        elif index % 3 == 0:
            treated = "0.01"
            control = "0.02"
            raw = "-0.01"
            corrected = "0"
            reason = "none"
        else:
            treated = "0.10"
            control = "0.02"
            raw = "0.08"
            corrected = "0.08"
            reason = "none"
        key = {"chrom": "chr1", "position": str(position), "ref": "C", "alt": "T"}
        label = {
            **key,
            "treated_median": treated,
            "control_median": control,
            "raw_edit_rate_difference": raw,
            "corrected_editing_efficiency": corrected,
            "training_eligible": "1" if eligible else "0",
            "label_confidence": confidence,
            "elevated_control_background": "1" if index == 6 else "0",
            "exclusion_reason": reason,
        }
        meta = {
            **key,
            "sequence_context": sequence,
            "sequence_length": "101",
            "center_index": "50",
            "center_base": "C",
            "orientation_qc": "pass",
            "transcript_oriented_ref": "C",
            "transcript_oriented_alt": "T",
            "corrected_editing_efficiency": corrected,
            "raw_edit_rate_difference": raw,
            "training_eligible": "1" if eligible else "0",
            "label_confidence": confidence,
            "exclusion_reason": reason,
        }
        labels.append(label)
        metadata.append(meta)
    return labels, metadata


def joined_eligible_rows():
    labels, metadata = synthetic_rows()
    label_fields = list(labels[0])
    metadata_fields = list(metadata[0])
    joined = handoff.merge_rows(
        label_fields,
        metadata_fields,
        handoff.index_unique(labels, "labels"),
        handoff.index_unique(metadata, "metadata"),
    )
    return handoff.construct_subsets(joined)["all_eligible"]


class LamarFinetuningHandoffTests(unittest.TestCase):
    def test_background_correction_formula(self):
        self.assertAlmostEqual(handoff.corrected_from_medians(0.30, 0.10), 0.20)

    def test_negative_raw_difference_becomes_corrected_zero(self):
        row = {
            "chrom": "chr1",
            "position": "10",
            "ref": "C",
            "alt": "T",
            "treated_median": "0.1",
            "control_median": "0.2",
            "raw_edit_rate_difference": "-0.1",
            "corrected_editing_efficiency": "0",
            "training_eligible": "1",
            "exclusion_reason": "none",
        }
        handoff.validate_label_row(row)
        self.assertEqual(handoff.corrected_from_medians(0.1, 0.2), 0)

    def test_missing_control_does_not_become_zero_label(self):
        row = {
            "chrom": "chr1",
            "position": "10",
            "ref": "C",
            "alt": "T",
            "treated_median": "0.1",
            "control_median": "NA",
            "raw_edit_rate_difference": "NA",
            "corrected_editing_efficiency": "NA",
            "training_eligible": "0",
            "exclusion_reason": "low_control_coverage",
        }
        handoff.validate_label_row(row)
        self.assertIsNone(handoff.corrected_from_medians(0.1, None))

    def test_sequence_length_and_center_c_validation(self):
        _, metadata = synthetic_rows()
        handoff.validate_sequence_row(metadata[0])
        broken = dict(metadata[0], center_base="A")
        with self.assertRaisesRegex(ValueError, "center base"):
            handoff.validate_sequence_row(broken)

    def test_overlap_cluster_construction(self):
        rows = joined_eligible_rows()
        clusters = handoff.build_overlap_clusters(rows)
        row_by_position = {int(row["position"]): row for row in rows}
        self.assertEqual(
            clusters[handoff.allele_key(row_by_position[2000])],
            clusters[handoff.allele_key(row_by_position[2050])],
        )
        self.assertNotEqual(
            clusters[handoff.allele_key(row_by_position[1000])],
            clusters[handoff.allele_key(row_by_position[1300])],
        )

    def test_duplicate_sequences_stay_in_one_split(self):
        split_rows = handoff.assign_splits(joined_eligible_rows(), seed=20260715)
        groups = {}
        for row in split_rows:
            groups.setdefault(row["duplicate_sequence_group_id"], set()).add(row["split"])
        self.assertTrue(all(len(splits) == 1 for splits in groups.values()))

    def test_deterministic_splitting(self):
        rows = joined_eligible_rows()
        first = {
            row["allele_key"]: row["split"] for row in handoff.assign_splits(rows, seed=20260715)
        }
        second = {
            row["allele_key"]: row["split"] for row in handoff.assign_splits(rows, seed=20260715)
        }
        self.assertEqual(first, second)

    def test_no_allele_key_leakage(self):
        rows = handoff.assign_splits(joined_eligible_rows(), seed=20260715)
        leaked = rows + [dict(rows[0], split="test")]
        with self.assertRaisesRegex(ValueError, "allele key"):
            handoff.validate_split_assignments(leaked)

    def test_no_overlap_cluster_leakage(self):
        rows = handoff.assign_splits(joined_eligible_rows(), seed=20260715)
        groups = {}
        for row in rows:
            groups.setdefault(row["overlap_cluster_id"], []).append(row)
        pair = next(group for group in groups.values() if len(group) > 1)
        leaked = copy.deepcopy(rows)
        target_key = pair[0]["allele_key"]
        for row in leaked:
            if row["allele_key"] == target_key:
                row["split"] = "test" if pair[1]["split"] != "test" else "validation"
        with self.assertRaisesRegex(ValueError, "overlap cluster"):
            handoff.validate_split_assignments(leaked)

    def test_no_duplicate_sequence_leakage(self):
        rows = handoff.assign_splits(joined_eligible_rows(), seed=20260715)
        groups = {}
        for row in rows:
            groups.setdefault(row["duplicate_sequence_group_id"], []).append(row)
        pair = next(group for group in groups.values() if len(group) > 1)
        leaked = copy.deepcopy(rows)
        target_key = pair[0]["allele_key"]
        for row in leaked:
            if row["allele_key"] == target_key:
                row["split"] = "test" if pair[1]["split"] != "test" else "validation"
                row["overlap_cluster_id"] = "independent_for_this_test"
        with self.assertRaisesRegex(ValueError, "identical sequence"):
            handoff.validate_split_assignments(leaked)

    def test_strict_and_all_eligible_subsets(self):
        labels, metadata = synthetic_rows()
        joined = handoff.merge_rows(
            list(labels[0]),
            list(metadata[0]),
            handoff.index_unique(labels, "labels"),
            handoff.index_unique(metadata, "metadata"),
        )
        subsets = handoff.construct_subsets(joined)
        self.assertEqual(len(subsets["all_eligible"]), 16)
        self.assertEqual(len(subsets["high_confidence"]), 14)
        self.assertEqual(len(subsets["high_confidence_low_control"]), 13)
        self.assertEqual(len(subsets["excluded"]), 2)

    def test_checksum_generation_and_full_synthetic_handoff(self):
        labels, metadata = synthetic_rows()
        with tempfile.TemporaryDirectory() as temporary:
            temporary = Path(temporary)
            labels_path = temporary / "background_corrected_labels.tsv.gz"
            metadata_path = temporary / "lamar_ready_metadata.tsv.gz"
            output = temporary / "handoff"
            handoff.write_tsv(labels_path, list(labels[0]), labels)
            handoff.write_tsv(metadata_path, list(metadata[0]), metadata)
            result = handoff.build_handoff(labels_path, metadata_path, output)
            self.assertEqual(result["counts"]["all_eligible"], 16)
            self.assertEqual(result["split_qc"]["validation_status"], "pass")
            self.assertTrue(handoff.verify_checksums(output))
            self.assertTrue(all((output / name).is_file() for name in handoff.EXPECTED_OUTPUTS))
            manifest = json.loads((output / "handoff_manifest.json").read_text())
            self.assertTrue(manifest["guardrails"]["missing_labels_are_not_zero"])
            scalar_path = temporary / "scalar.tsv.gz"
            exported = exporter.export_scalar(
                output / "CU5.17_lamar_splits.tsv.gz", scalar_path
            )
            self.assertEqual(exported, 14)
            scalar_fields, scalar_rows = handoff.read_tsv(scalar_path)
            self.assertEqual(tuple(scalar_fields), exporter.OUTPUT_FIELDS)
            self.assertTrue(all(row["center_index"] == "50" for row in scalar_rows))

    def test_puf_target_is_required_for_token_mask_export(self):
        with self.assertRaisesRegex(ValueError, "experimentally confirmed"):
            exporter.validate_puf_target_requirement(Path("token.tsv.gz"), None)
        self.assertEqual(
            exporter.validate_puf_target_requirement(Path("token.tsv.gz"), "ACGTU"), "ACGTU"
        )
        masked = exporter.token_mask_rows(
            [
                {
                    "sequence": unique_sequence(1),
                    "center_index": "50",
                    "corrected_editing_efficiency": "0.08",
                    "split": "train",
                    "chrom": "chr1",
                    "position": "100",
                    "ref": "C",
                    "alt": "T",
                }
            ],
            "ACGTU",
        )[0]
        mask = json.loads(masked["label_mask"])
        values = json.loads(masked["label_values"])
        self.assertEqual(sum(mask), 1)
        self.assertEqual(mask[50], 1)
        self.assertEqual(values[50], 0.08)
        self.assertTrue(all(value is None for index, value in enumerate(values) if index != 50))
        self.assertNotIn("label_total_count", masked)

    @unittest.skipUnless(SKLEARN_AVAILABLE, "optional scikit-learn baseline dependency is absent")
    def test_scalar_baseline_uses_generated_split(self):
        labels, metadata = synthetic_rows()
        with tempfile.TemporaryDirectory() as temporary:
            temporary = Path(temporary)
            labels_path = temporary / "labels.tsv.gz"
            metadata_path = temporary / "metadata.tsv.gz"
            output = temporary / "handoff"
            handoff.write_tsv(labels_path, list(labels[0]), labels)
            handoff.write_tsv(metadata_path, list(metadata[0]), metadata)
            handoff.build_handoff(labels_path, metadata_path, output)
            rows = baseline.read_split_rows(
                output / "CU5.17_lamar_splits.tsv.gz", "high_confidence"
            )
            report = baseline.run_baselines(rows, kmer_size=2, alpha=1.0)
            self.assertEqual(set(report["row_counts"]), {"train", "validation", "test"})
            self.assertIn("kmer_ridge", report["metrics"])


if __name__ == "__main__":
    unittest.main()

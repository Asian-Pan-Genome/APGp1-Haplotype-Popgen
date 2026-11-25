import json
import random
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from flanking_block_ld import background_block_ld as workflow
from flanking_block_ld.tests import test_background_block_ld as existing_tests


class BackgroundRegionPairPTests(unittest.TestCase):
    def test_cli_default_and_skip_flag(self):
        argv = ["target.bed", "background.bed", "cohort.vcf", "-o", "output"]
        defaults = workflow.build_parser().parse_args(argv)
        self.assertFalse(defaults.skip_background_region_pair_p)
        self.assertFalse(defaults.skip_pair_test)
        self.assertFalse(defaults.skip_region_test)
        skipped = workflow.build_parser().parse_args(
            argv + ["--skip-background-region-pair-p", "--skip-pair-test", "--skip-region-test"]
        )
        self.assertTrue(skipped.skip_background_region_pair_p)
        self.assertTrue(skipped.skip_pair_test)
        self.assertTrue(skipped.skip_region_test)

    def test_max_only_analysis_does_not_require_pair_p(self):
        row = existing_tests.BackgroundBlockLDTests.pair_row(1, 1, 0.9, 0.001)
        del row["pair_nominal_p"]
        stats = workflow.analyze_region([row], [(1, 1)], 0.01, calculate_breadth=False)
        self.assertEqual(stats["max_dprime"], 0.9)
        self.assertEqual(stats["n_valid_pairs"], 1)
        self.assertEqual(stats["n_eligible_left_blocks"], 1)
        self.assertEqual(stats["eligible_left_window_coverage_bp"], 10)
        for field in ("n_significant_pairs", "significant_pair_fraction", "peak_pair_nominal_p",
                      "n_supported_left_blocks", "supported_left_block_fraction",
                      "supported_left_window_coverage_bp", "supported_left_window_coverage_fraction"):
            self.assertIsNone(stats[field], field)
            self.assertEqual(workflow.serialized(stats)[field], "NA", field)
        with self.assertRaisesRegex(ValueError, "calibrated pair background"):
            workflow.analyze_region([row], [(1, 1)], 0.01)

    def test_max_only_analysis_still_enforces_qc_and_snp_counts(self):
        row = existing_tests.BackgroundBlockLDTests.pair_row(1, 1, 0.9, 0.001)
        row.update(block1_snp_number_pass_threshold="1", block2_snp_number_pass_threshold="2")
        required = {("left", 1): 1, ("right", 1): 2}
        workflow.analyze_region([row], [(1, 1)], 0.01, required, calculate_breadth=False)
        row["block2_snp_number_pass_threshold"] = "1"
        with self.assertRaisesRegex(ValueError, "SNP count matching"):
            workflow.analyze_region([row], [(1, 1)], 0.01, required, calculate_breadth=False)
        row["status"] = "insufficient_right_snps"
        with self.assertRaisesRegex(ValueError, "eligible pair mask"):
            workflow.analyze_region([row], [(1, 1)], 0.01, required, calculate_breadth=False)

    @staticmethod
    def write_inputs(directory):
        (directory / "target.bed").write_text("1\t300\t400\ttarget\n", encoding="utf-8")
        (directory / "background.bed").write_text("1\t0\t1000\tbackground\n", encoding="utf-8")
        rows = [
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4",
        ]
        rng = random.Random(812)
        for pos in range(5, 1000, 5):
            parity_mask = rng.randint(1, 7)
            if 280 <= pos - 1 < 300 or 400 <= pos - 1 < 420:
                parity_mask = 1
            alleles = [bin(haplotype & parity_mask).count("1") % 2 for haplotype in range(8)]
            genotypes = "\t".join(f"{alleles[i]}|{alleles[i + 1]}" for i in range(0, 8, 2))
            rows.append(f"1\t{pos}\tv{pos}\tA\tG\t.\tPASS\t.\tGT\t{genotypes}")
        (directory / "cohort.vcf").write_text("\n".join(rows) + "\n", encoding="utf-8")

    def test_skip_mode_preserves_target_and_max_test_but_disables_background_breadth(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            self.write_inputs(directory)
            for prefix, skip in (("enabled", False), ("disabled", True)):
                args = existing_tests.BackgroundBlockLDTests.args(directory, prefix)
                args.background_pairs_per_distance = 40
                args.background_regions = 20
                args.pair_significance_threshold = 0.3
                args.skip_background_region_pair_p = skip
                with mock.patch.object(workflow, "attach_pair_p", wraps=workflow.attach_pair_p) as attach:
                    workflow.run(args)
                if skip:
                    self.assertEqual(attach.call_count, 1)  # Target pairs only.
                    self.assertTrue(all(row["region_id"] == "target" for row in attach.call_args.args[0]))
                else:
                    self.assertGreater(attach.call_count, 1)
            for suffix in ("target_pairs.tsv", "background_pairs.tsv", "background_regions.bed"):
                self.assertEqual((directory / f"enabled.{suffix}").read_text(),
                                 (directory / f"disabled.{suffix}").read_text(), suffix)
            enabled = existing_tests.BackgroundBlockLDTests.read_rows(directory / "enabled.summary.tsv")[0]
            disabled = existing_tests.BackgroundBlockLDTests.read_rows(directory / "disabled.summary.tsv")[0]
            for field in ("max_dprime", "region_max_empirical_p", "peak_pair_nominal_p",
                          "n_significant_pairs", "significant_pair_fraction", "n_supported_left_blocks",
                          "supported_right_window_coverage_fraction"):
                self.assertEqual(disabled[field], enabled[field], field)
                self.assertNotEqual(disabled[field], "NA", field)
            self.assertGreater(int(disabled["n_significant_pairs"]), 0)
            self.assertEqual(disabled["region_breadth_empirical_p"], "NA")
            self.assertEqual(disabled["region_breadth_test_status"], "disabled_background_region_pair_p")
            self.assertEqual(disabled["status"], "ok")
            self.assertEqual(disabled["background_region_pair_p_enabled"], "false")
            self.assertNotEqual(enabled["region_breadth_empirical_p"], "NA")
            for row in existing_tests.BackgroundBlockLDTests.read_rows(directory / "disabled.background_regions.tsv"):
                for field in ("peak_pair_nominal_p", "n_significant_pairs", "significant_pair_fraction",
                              "n_supported_left_blocks", "n_supported_right_blocks",
                              "supported_left_window_coverage_fraction"):
                    self.assertEqual(row[field], "NA", field)
                self.assertEqual(row["background_region_pair_p_enabled"], "false")
                self.assertEqual(int(row["n_compared_pairs"]), 4)
                self.assertGreaterEqual(float(row["max_dprime"]), 0)
            metadata = json.loads((directory / "disabled.run.json").read_text())
            self.assertFalse(metadata["background_region_pair_p_enabled"])
            self.assertFalse(metadata["region_breadth_test_enabled"])
            self.assertEqual(metadata["status"], "complete")

    def test_pair_and_region_tests_can_be_skipped_independently(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            self.write_inputs(directory)

            pair_skipped = existing_tests.BackgroundBlockLDTests.args(directory, "pair_skipped")
            pair_skipped.skip_pair_test = True
            with mock.patch.object(
                    workflow, "random_pair_windows",
                    side_effect=AssertionError("pair background sampling must be skipped")):
                pair_result = workflow.run(pair_skipped)
            self.assertEqual(pair_result["background_pair_count"], 0)
            self.assertEqual(pair_result["background_region_count"], pair_skipped.background_regions)
            self.assertEqual(existing_tests.BackgroundBlockLDTests.read_rows(
                directory / "pair_skipped.background_pairs.tsv"), [])
            pair_target_rows = existing_tests.BackgroundBlockLDTests.read_rows(
                directory / "pair_skipped.target_pairs.tsv")
            self.assertTrue(all(row["multiallelic_dprime"] != "NA" for row in pair_target_rows))
            self.assertTrue(all(row["pair_nominal_p"] == "NA" for row in pair_target_rows))
            pair_summary = existing_tests.BackgroundBlockLDTests.read_rows(
                directory / "pair_skipped.summary.tsv")[0]
            self.assertEqual(pair_summary["pair_test_enabled"], "false")
            self.assertEqual(pair_summary["region_test_enabled"], "true")
            self.assertNotEqual(pair_summary["region_max_empirical_p"], "NA")
            self.assertEqual(pair_summary["region_breadth_empirical_p"], "NA")
            self.assertEqual(pair_summary["significant_pair_fraction"], "NA")
            self.assertEqual(pair_summary["region_breadth_test_status"], "disabled_pair_test")
            self.assertEqual(pair_summary["status"], "ok")

            region_skipped = existing_tests.BackgroundBlockLDTests.args(directory, "region_skipped")
            region_skipped.skip_region_test = True
            with mock.patch.object(
                    workflow, "random_region_intervals",
                    side_effect=AssertionError("region background sampling must be skipped")):
                region_result = workflow.run(region_skipped)
            self.assertEqual(region_result["background_pair_count"], 24)
            self.assertEqual(region_result["background_region_count"], 0)
            self.assertEqual(existing_tests.BackgroundBlockLDTests.read_rows(
                directory / "region_skipped.background_regions.tsv"), [])
            region_target_rows = existing_tests.BackgroundBlockLDTests.read_rows(
                directory / "region_skipped.target_pairs.tsv")
            self.assertTrue(all(row["pair_nominal_p"] != "NA" for row in region_target_rows))
            region_summary = existing_tests.BackgroundBlockLDTests.read_rows(
                directory / "region_skipped.summary.tsv")[0]
            self.assertEqual(region_summary["pair_test_enabled"], "true")
            self.assertEqual(region_summary["region_test_enabled"], "false")
            self.assertEqual(region_summary["region_max_empirical_p"], "NA")
            self.assertEqual(region_summary["region_breadth_empirical_p"], "NA")
            self.assertNotEqual(region_summary["significant_pair_fraction"], "NA")
            self.assertEqual(region_summary["status"], "ok_region_test_skipped")

            both_skipped = existing_tests.BackgroundBlockLDTests.args(directory, "both_skipped")
            both_skipped.skip_pair_test = True
            both_skipped.skip_region_test = True
            with mock.patch.object(
                    workflow, "random_pair_windows",
                    side_effect=AssertionError("pair background sampling must be skipped")), \
                    mock.patch.object(
                        workflow, "random_region_intervals",
                        side_effect=AssertionError("region background sampling must be skipped")):
                both_result = workflow.run(both_skipped)
            self.assertEqual(both_result["background_pair_count"], 0)
            self.assertEqual(both_result["background_region_count"], 0)
            both_summary = existing_tests.BackgroundBlockLDTests.read_rows(
                directory / "both_skipped.summary.tsv")[0]
            self.assertEqual(both_summary["max_dprime"], pair_summary["max_dprime"])
            self.assertEqual(both_summary["region_max_empirical_p"], "NA")
            self.assertEqual(both_summary["significant_pair_fraction"], "NA")
            self.assertEqual(both_summary["status"], "ok_region_test_skipped")
            metadata = json.loads((directory / "both_skipped.run.json").read_text())
            self.assertFalse(metadata["pair_test_enabled"])
            self.assertFalse(metadata["region_test_enabled"])
            self.assertEqual(metadata["pair_sampling"], "skipped_by_user")
            self.assertEqual(metadata["region_sampling"], "skipped_by_user")


if __name__ == "__main__":
    unittest.main()

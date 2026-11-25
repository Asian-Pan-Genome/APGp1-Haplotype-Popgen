import contextlib
import io
import json
import math
import tempfile
import unittest
from pathlib import Path

from flanking_block_ld import background_block_ld as workflow
from flanking_block_ld.flanking_block_ld import Region, common_haplotype_keys, read_vcf
from flanking_block_ld.tests import test_background_block_ld as existing_tests
from flanking_block_ld.tests import test_background_region_pair_p as skip_tests


class PrefilteredVCFTests(unittest.TestCase):
    def test_cli_requires_explicit_prefiltered_declaration_and_supports_alias(self):
        argv = ["target.bed", "background.bed", "cohort.vcf", "-o", "output"]
        self.assertFalse(workflow.build_parser().parse_args(argv).vcf_prefiltered)
        for flag in ("--vcf-prefiltered", "--skip-vcf-maf-missing-qc"):
            self.assertTrue(workflow.build_parser().parse_args(argv + [flag]).vcf_prefiltered)

    def test_reader_skips_rate_calculation_and_filtering_but_keeps_other_checks(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "cohort.vcf"
            path.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n"
                "1\t95\tlow_maf\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|0\t0|0\t0|0\n"
                "1\t96\tmissing\tA\tG\t.\tPASS\t.\tGT\t.|.\t0|1\t0|1\t0|1\n"
                "1\t97\tgood\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|1\t0|1\t0|1\n"
                "1\t98\tall_missing\tA\tG\t.\tPASS\t.\tGT\t.|.\t.|.\t.|.\t.|.\n"
                "1\t99\tmonomorphic\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|0\t0|0\t0|0\n"
                "1\t100\tunphased\tA\tG\t.\tPASS\t.\tGT\t0/1\t0|1\t0|1\t0|1\n"
                "1\t101\tindel\tA\tGG\t.\tPASS\t.\tGT\t0|1\t0|1\t0|1\t0|1\n"
                "1\t102\tfailed_filter\tA\tG\t.\tLowQual\t.\tGT\t0|1\t0|1\t0|1\t0|1\n",
                encoding="utf-8")
            params = dict(path=path, regions=[Region("1", "1", 0, 200, "all", 1)],
                          flank_bp=0, min_maf=0.2, min_call_rate=0.95,
                          include_filtered=False, allow_unphased_missing=False,
                          include_region_interiors=True)
            normal = read_vcf(**params)
            self.assertFalse(normal.maf_missing_qc_skipped)
            self.assertEqual([v.variant_id for v in normal.variants_by_chrom["1"]], ["good"])
            self.assertEqual(normal.filtered_maf_count, 2)
            self.assertEqual(normal.filtered_call_rate_count, 2)
            updates = []
            trusted = read_vcf(**params, skip_maf_missing_qc=True, progress_callback=updates.append)
            self.assertTrue(trusted.maf_missing_qc_skipped)
            self.assertEqual([v.variant_id for v in trusted.variants_by_chrom["1"]],
                             ["low_maf", "missing", "good", "all_missing", "monomorphic"])
            self.assertTrue(all(math.isnan(v.maf) and math.isnan(v.call_rate)
                                for v in trusted.variants_by_chrom["1"]))
            self.assertEqual(trusted.filtered_maf_count, 0)  # Not executed, not proof of zero failures.
            self.assertEqual(trusted.filtered_call_rate_count, 0)
            self.assertEqual(trusted.filtered_unphased_count, 1)
            self.assertEqual(trusted.filtered_non_snp_count, 1)
            self.assertEqual(trusted.filtered_vcf_filter_count, 1)
            self.assertEqual(trusted.pre_missing_snp_positions_by_chrom, {})
            self.assertEqual(trusted.missing_filtered_snp_positions_by_chrom, {})
            self.assertEqual(updates[-1]["records_scanned"], 8)
            self.assertEqual(updates[-1]["retained_snps"], 5)
            self.assertEqual(updates[-1]["maf_missing_qc_skipped"], 1)
            self.assertEqual(common_haplotype_keys((trusted.variants_by_chrom["1"][3],)), ())

    def test_prefiltered_reader_still_validates_gt(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "bad.vcf"
            path.write_text("##fileformat=VCFv4.2\n"
                            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
                            "1\t95\tbad\tA\tG\t.\tPASS\t.\tGT\t0|2\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "unsupported GT allele"):
                read_vcf(path, [Region("1", "1", 0, 200, "all", 1)], 0, 0.1, 0.95,
                         False, False, include_region_interiors=True, skip_maf_missing_qc=True)

    def test_progress_marks_skipped_missing_qc_as_na(self):
        with tempfile.TemporaryDirectory() as temp_dir, contextlib.redirect_stderr(io.StringIO()):
            path = Path(temp_dir) / "vcf.log"
            logger = workflow.RunLogger(path)
            try:
                logger.vcf_progress(dict(records_scanned=1000, candidate_variants=80, retained_snps=75,
                                         filtered_missing=0, maf_missing_qc_skipped=1, finished=1))
                log = path.read_text()
                self.assertIn("VCF DONE records_scanned=1000", log)
                self.assertIn("filtered_missing=NA maf_missing_qc=skipped", log)
                self.assertNotIn("filtered_missing=0", log)
            finally:
                logger.close()

    def test_prefiltered_workflow_preserves_ld_and_marks_unknown_missing_counts(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            skip_tests.BackgroundRegionPairPTests.write_inputs(directory)
            for prefix, prefiltered in (("normal", False), ("trusted", True)):
                args = existing_tests.BackgroundBlockLDTests.args(directory, prefix)
                args.vcf_prefiltered = prefiltered
                args.skip_region_entropy_matching = False
                with contextlib.redirect_stderr(io.StringIO()):
                    workflow.run(args)
            read_rows = existing_tests.BackgroundBlockLDTests.read_rows
            normal = read_rows(directory / "normal.summary.tsv")[0]
            trusted = read_rows(directory / "trusted.summary.tsv")[0]
            for field in workflow.STAT_FIELDS + ("region_max_empirical_p", "region_breadth_empirical_p"):
                self.assertEqual(normal[field], trusted[field], field)
            self.assertEqual(trusted["vcf_maf_missing_qc_skipped"], "true")
            self.assertEqual(trusted["window_missing_qc_status"], "not_assessable_prefiltered_vcf")
            self.assertEqual(trusted["region_entropy_matching_enabled"], "true")
            self.assertEqual((directory / "normal.background_regions.bed").read_text(),
                             (directory / "trusted.background_regions.bed").read_text())
            for suffix in ("target_pairs.tsv", "background_pairs.tsv"):
                for row in read_rows(directory / f"trusted.{suffix}"):
                    for field in ("left_pre_missing_snp_count", "left_missing_filtered_snp_count",
                                  "right_missing_filtered_fraction", "block1_missed_snp_number",
                                  "block2_missed_snp_number", "total_missed_snp_number"):
                        self.assertEqual(row[field], "NA", field)
                    self.assertEqual(row["status"], "ok")
                    self.assertEqual(row["window_missing_qc_status"], "not_assessable_prefiltered_vcf")
            metadata = json.loads((directory / "trusted.run.json").read_text())
            self.assertTrue(metadata["vcf_qc"]["maf_missing_qc_skipped"])
            self.assertFalse(metadata["vcf_qc"]["site_maf_filter_enabled"])
            self.assertFalse(metadata["vcf_qc"]["site_missing_filter_enabled"])
            self.assertIsNone(metadata["vcf_qc"]["filtered_maf_count"])
            self.assertIsNone(metadata["vcf_qc"]["filtered_site_missing_count"])
            self.assertEqual(metadata["vcf_qc"]["rate_calculations"], "skipped")
            log = (directory / "trusted.log").read_text()
            self.assertIn("VCF PREFILTERED", log)
            self.assertIn("filtered_missing=NA maf_missing_qc=skipped", log)
            self.assertIn("cannot be reconstructed", log)


if __name__ == "__main__":
    unittest.main()

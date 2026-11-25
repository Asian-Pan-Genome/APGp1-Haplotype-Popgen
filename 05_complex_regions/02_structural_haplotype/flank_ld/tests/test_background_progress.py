import contextlib
import io
import json
import random
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from flanking_block_ld import background_block_ld as workflow
from flanking_block_ld.flanking_block_ld import Region, read_vcf
from flanking_block_ld.tests import test_background_block_ld as existing_tests
from flanking_block_ld.tests import test_background_region_pair_p as skip_tests


class BackgroundProgressTests(unittest.TestCase):
    def test_cli_defaults_and_invalid_logging_options(self):
        argv = ["target.bed", "background.bed", "cohort.vcf", "-o", "output"]
        args = workflow.build_parser().parse_args(argv)
        self.assertIsNone(args.log_file)
        self.assertEqual(args.log_interval_seconds, 30.0)
        self.assertEqual(args.log_every_attempts, 1000)
        for flag, value in (("--log-interval-seconds", "0"), ("--log-interval-seconds", "nan"),
                            ("--log-interval-seconds", "inf"), ("--log-every-attempts", "0")):
            with self.subTest(flag=flag, value=value), contextlib.redirect_stderr(io.StringIO()), \
                    self.assertRaises(SystemExit) as error:
                workflow.main(argv + [flag, value])
            self.assertEqual(error.exception.code, 2)

    def test_sampling_reports_zero_acceptance_time_trigger_and_completion_and_flushes(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "progress.log"
            clock = [0.0]
            stderr = io.StringIO()
            with mock.patch.object(workflow.time, "monotonic", side_effect=lambda: clock[0]), \
                    contextlib.redirect_stderr(stderr):
                logger = workflow.RunLogger(path, interval_seconds=10, every_attempts=2)
                try:
                    tracker = workflow.SamplingProgress(logger, "PAIR stratum=1/1", 4, 100, 20)
                    tracker.update(0, 1, {"insufficient_snps": 1})
                    self.assertNotIn("PROGRESS", path.read_text())
                    clock[0] = 5
                    tracker.update(0, 2, {"insufficient_snps": 2})
                    self.assertIn("accepted=0/4 percent=0.0% attempts=2/20", path.read_text())
                    self.assertIn("eta=NA", path.read_text())
                    clock[0] = 20  # Only one more attempt, but time trigger fires.
                    reasons = {"insufficient_snps": 2, "selected_snp_count_mismatch": 1}
                    tracker.update(0, 3, reasons)
                    self.assertIn("attempts=3/20 rejected=3 qc_rejected=2 snp_count_rejected=1", path.read_text())
                    clock[0] = 25
                    tracker.update(1, 5, dict(reasons, insufficient_snps=3))
                    self.assertIn("eta=75.0s", path.read_text())
                    tracker.update(4, 8, dict(reasons, insufficient_snps=3))
                    self.assertIn("DONE accepted=4/4 percent=100.0%", path.read_text())
                    self.assertIn("eta=0.0s", path.read_text())
                    # Reading before close checks that each message was flushed.
                    self.assertEqual(stderr.getvalue(), path.read_text())
                    self.assertRegex(path.read_text(), r"\[\d{4}-\d{2}-\d{2}T.*\] \[elapsed=0.0s\] \[INFO\]")
                finally:
                    logger.close()

    def test_vcf_callback_counts_all_records_and_keeps_results_unchanged(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "cohort.vcf"
            rows = ["##fileformat=VCFv4.2",
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1"]
            for index in range(1, 2502):
                chrom, pos, gt = "2", index, "0|1"
                if index in (1, 1500, 2501):
                    chrom, pos = "1", {1: 1, 1500: 2, 2501: 3}[index]
                    if index == 1500:
                        gt = ".|."
                rows.append(f"{chrom}\t{pos}\tv{index}\tA\tG\t.\tPASS\t.\tGT\t{gt}")
            path.write_text("\n".join(rows) + "\n", encoding="utf-8")
            params = dict(path=path, regions=[Region("1", "1", 0, 10, "test", 1)], flank_bp=0,
                          min_maf=0.1, min_call_rate=0.95, include_filtered=False,
                          allow_unphased_missing=False, include_region_interiors=True)
            updates = []
            baseline = read_vcf(**params)
            logged = read_vcf(**params, progress_callback=updates.append)
            self.assertEqual(baseline, logged)
            self.assertEqual([item["records_scanned"] for item in updates], [1000, 2000, 2501])
            self.assertEqual([item["candidate_variants"] for item in updates], [1, 2, 3])
            self.assertEqual([item["retained_snps"] for item in updates], [1, 1, 2])
            self.assertEqual(updates[-1]["filtered_missing"], 1)
            self.assertEqual(updates[-1]["finished"], 1)

    def test_vcf_time_throttle_and_forced_completion(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "vcf.log"
            clock = [0.0]
            with mock.patch.object(workflow.time, "monotonic", side_effect=lambda: clock[0]), \
                    contextlib.redirect_stderr(io.StringIO()):
                logger = workflow.RunLogger(path)
                try:
                    counts = dict(records_scanned=1000, candidate_variants=20, retained_snps=10,
                                  filtered_missing=2, finished=0)
                    logger.vcf_progress(counts)
                    self.assertEqual(path.read_text(), "")
                    clock[0] = 30
                    logger.vcf_progress(counts)
                    logger.vcf_progress(dict(counts, records_scanned=1001, finished=1))
                    text = path.read_text()
                    self.assertEqual(text.count("VCF PROGRESS"), 1)
                    self.assertEqual(text.count("VCF DONE"), 1)
                    self.assertNotIn("percent=", text)  # No known total for streaming VCF.
                finally:
                    logger.close()

    @staticmethod
    def range_inputs(directory):
        (directory / "target.bed").write_text("1\t300\t400\ttarget\n", encoding="utf-8")
        (directory / "background.bed").write_text("1\t0\t1000\tbackground\n", encoding="utf-8")
        rows = ["##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4"]
        rng = random.Random(201)
        for start in range(0, 1000, 20):
            for offset in (4, 14, 17):
                pos = start + offset
                mask = rng.randint(1, 7)
                alleles = [bin(haplotype & mask).count("1") % 2 for haplotype in range(8)]
                gt = "\t".join(f"{alleles[i]}|{alleles[i + 1]}" for i in range(0, 8, 2))
                rows.append(f"1\t{pos + 1}\tv{pos}\tA\tG\t.\tPASS\t.\tGT\t{gt}")
        (directory / "cohort.vcf").write_text("\n".join(rows) + "\n", encoding="utf-8")

    def test_full_workflow_logs_sampling_rejections_and_frequency_does_not_change_statistics(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            self.range_inputs(directory)
            for prefix in ("normal", "verbose"):
                args = existing_tests.BackgroundBlockLDTests.args(directory, prefix)
                args.snps_per_block = None
                args.min_snps_per_block, args.max_snps_per_block = 1, 2
                args.max_attempt_factor = 100
                if prefix == "verbose":
                    args.log_every_attempts = 1
                    args.log_interval_seconds = 0.000001
                    args.log_file = directory / "custom" / "verbose.log"
                with contextlib.redirect_stderr(io.StringIO()) as stderr:
                    workflow.run(args)
                log = args.log_file if prefix == "verbose" else directory / "normal.log"
                self.assertEqual(stderr.getvalue(), log.read_text())
            for suffix in ("target_pairs.tsv", "background_pairs.tsv", "background_regions.tsv",
                           "background_regions.bed", "summary.tsv"):
                self.assertEqual((directory / f"normal.{suffix}").read_text(),
                                 (directory / f"verbose.{suffix}").read_text(), suffix)
            text = (directory / "custom" / "verbose.log").read_text()
            for stage in ("RUN START", "BED DONE", "VCF DONE", "TARGET 1/1", "PAIR BACKGROUND PLAN",
                          "PAIR BACKGROUND TOTAL", "REGION BACKGROUND TOTAL", "RUN COMPLETE"):
                self.assertIn(stage, text)
            self.assertRegex(text, r"PAIR stratum=.* PROGRESS accepted=\d+/8")
            self.assertRegex(text, r"REGION target=target .* PROGRESS accepted=\d+/5")
            self.assertIn('"selected_snp_count_mismatch":', text)
            self.assertIn("rank_pair_mask_matching=false pair_distance_matching=false", text)
            metadata = json.loads((directory / "verbose.run.json").read_text())
            self.assertEqual(metadata["outputs"]["log"], str(directory / "custom" / "verbose.log"))
            self.assertEqual(metadata["logging"]["every_attempts"], 1)

    def test_failure_keeps_zero_acceptance_progress_and_error_for_both_sampling_stages(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            skip_tests.BackgroundRegionPairPTests.write_inputs(directory)
            original_summarize = workflow.BlockStore.summarize
            original_screen_region = workflow.BlockStore.screen_region

            def reject_pair(store, region, first, second):
                row = original_summarize(store, region, first, second)
                if region.region_id.startswith("pair_"):
                    row["status"] = "below_min_complete_haplotypes"
                return row

            def reject_region(store, region):
                if "|bg_" in region.region_id:
                    raise ValueError("background region fails target available block count matching")
                return original_screen_region(store, region)

            for prefix, method, replacement, label in (
                ("pair_fail", "summarize", reject_pair, "PAIR stratum="),
                ("region_fail", "screen_region", reject_region, "REGION target="),
            ):
                with self.subTest(stage=prefix):
                    args = existing_tests.BackgroundBlockLDTests.args(directory, prefix)
                    args.background_pairs_per_distance = args.background_regions = 2
                    args.max_attempt_factor = 2
                    args.log_every_attempts = 1
                    with mock.patch.object(workflow.BlockStore, method, replacement), \
                            contextlib.redirect_stderr(io.StringIO()), self.assertRaises(ValueError):
                        workflow.run(args)
                    text = (directory / f"{prefix}.log").read_text()
                    updates = [line for line in text.splitlines() if label in line and "PROGRESS" in line]
                    self.assertTrue(updates)
                    self.assertTrue(all("accepted=0/2" in line and "eta=NA" in line for line in updates))
                    self.assertIn("EXHAUSTED accepted=0/2", text)
                    self.assertIn("[ERROR] RUN FAILED", text)
                    self.assertNotIn("RUN COMPLETE", text)
                    self.assertFalse((directory / f"{prefix}.run.json").exists())

    def test_log_paths_do_not_overwrite_inputs_outputs_or_existing_logs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            skip_tests.BackgroundRegionPairPTests.write_inputs(directory)
            args = existing_tests.BackgroundBlockLDTests.args(directory, "safe")
            args.log_file = args.vcf
            old_vcf = args.vcf.read_text()
            with self.assertRaisesRegex(ValueError, "cannot be input"):
                workflow.run(args)
            self.assertEqual(args.vcf.read_text(), old_vcf)
            args.log_file = directory / "safe.summary.tsv"
            with self.assertRaisesRegex(ValueError, "must be distinct"):
                workflow.run(args)
            args.log_file = directory / "existing.log"
            args.log_file.write_text("old log\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "log file already exists"):
                workflow.run(args)
            self.assertEqual(args.log_file.read_text(), "old log\n")


if __name__ == "__main__":
    unittest.main()

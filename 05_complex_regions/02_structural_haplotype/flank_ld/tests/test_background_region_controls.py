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


class BackgroundRegionControlTests(unittest.TestCase):
    def test_block_count_cli_defaults_and_strict_override(self):
        argv = ["target.bed", "background.bed", "cohort.vcf", "-o", "output"]
        args = workflow.build_parser().parse_args(argv)
        self.assertEqual(args.region_block_count_tolerance, 1)
        self.assertTrue(args.region_block_count_allow_swap)
        strict = workflow.build_parser().parse_args(
            argv + ["--region-block-count-tolerance", "0", "--no-region-block-count-swap"])
        self.assertEqual(strict.region_block_count_tolerance, 0)
        self.assertFalse(strict.region_block_count_allow_swap)
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as error:
            workflow.main(argv + ["--region-block-count-tolerance", "-1"])
        self.assertEqual(error.exception.code, 2)

        self.assertEqual(args.region_match_fallback, "nearest")
        self.assertIsNone(args.background_region_attempts)
        self.assertEqual(args.region_fallback_max_score, 2.0)
        configured = workflow.build_parser().parse_args(
            argv + ["--background-region-attempts", "4000",
                    "--region-match-fallback", "none",
                    "--region-fallback-max-score", "1.5"])
        self.assertEqual(configured.background_region_attempts, 4000)
        self.assertEqual(configured.region_match_fallback, "none")
        self.assertEqual(configured.region_fallback_max_score, 1.5)

    def test_count_tolerance_swap_boundaries_ties_and_zero_pair_guard(self):
        pair = existing_tests.BackgroundBlockLDTests.pair_row

        def rows_for(counts):
            return [pair(left, right, 0.5, 0.2)
                    for left in range(1, counts[0] + 1) for right in range(1, counts[1] + 1)]

        for actual, orientation in (((8, 9), "direct"), ((9, 8), "swapped"),
                                    ((7, 10), "direct"), ((10, 7), "swapped"),
                                    ((8, 10), "direct"), ((10, 8), "swapped"),
                                    ((9, 9), "direct"), ((7, 7), None),
                                    ((10, 10), None), ((6, 9), None)):
            with self.subTest(actual=actual):
                if orientation is None:
                    with self.assertRaisesRegex(ValueError, "available block count matching"):
                        workflow.check_available_block_counts(rows_for(actual), (8, 9))
                else:
                    self.assertEqual(workflow.check_available_block_counts(rows_for(actual), (8, 9)), orientation)
        self.assertEqual(workflow.check_available_block_counts(rows_for((9, 8)), (8, 9), tolerance=0), "swapped")
        for actual, tolerance, allow_swap in (((10, 8), 1, False), ((9, 8), 0, False), ((7, 10), 0, True)):
            with self.assertRaisesRegex(ValueError, "available block count matching"):
                workflow.check_available_block_counts(rows_for(actual), (8, 9), tolerance, allow_swap)
        with self.assertRaisesRegex(ValueError, "available block count matching"):
            workflow.check_available_block_counts([], (1, 1), tolerance=100)
        with self.assertRaisesRegex(ValueError, "tolerance must be >= 0"):
            workflow.check_available_block_counts(rows_for((8, 9)), (8, 9), tolerance=-1)

    def test_workflow_accepts_tolerance_and_swap_and_reports_aligned_deltas(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            skip_tests.BackgroundRegionPairPTests.write_inputs(directory)
            args = existing_tests.BackgroundBlockLDTests.args(directory, "tolerant")
            args.block_count = 3
            args.max_flank_bp = 30
            args.background_regions = 3
            args.region_block_count_tolerance = 1
            args.region_block_count_allow_swap = True
            args.skip_background_region_pair_p = True
            args.log_every_attempts = 1
            original = workflow.BlockStore.summarize
            original_eligibility = workflow.BlockStore.pair_is_eligible_without_dprime

            def controlled_rows(store, region, first, second):
                row = original(store, region, first, second)
                if region.region_id == "target" and first.rank == 3:
                    row["status"] = "insufficient_complete_haplotypes"  # Target counts = 2,3.
                elif "|bg_" in region.region_id:
                    if ((region.start == 30 and (first.rank != 1 or second.rank != 1)) or
                            (region.start == 40 and second.rank == 3) or
                            (region.start == 45 and first.rank != 1)):
                        row["status"] = "insufficient_complete_haplotypes"
                return row

            def controlled_eligibility(store, region, first, second):
                if not original_eligibility(store, region, first, second):
                    return False
                return not ("|bg_" in region.region_id and
                            ((region.start == 30 and (first.rank != 1 or second.rank != 1)) or
                             (region.start == 40 and second.rank == 3) or
                             (region.start == 45 and first.rank != 1)))

            def candidates(domains, limit, rng):
                for start in (30, 35, 40, 45):
                    yield domains[0], start, 100

            with mock.patch.object(workflow.BlockStore, "summarize", controlled_rows), \
                    mock.patch.object(workflow.BlockStore, "pair_is_eligible_without_dprime",
                                      controlled_eligibility), \
                    mock.patch.object(workflow, "random_region_intervals", candidates), \
                    contextlib.redirect_stderr(io.StringIO()):
                workflow.run(args)
            read_rows = existing_tests.BackgroundBlockLDTests.read_rows
            summary = read_rows(directory / "tolerant.summary.tsv")[0]
            self.assertEqual(summary["background_region_attempts"], "4")
            self.assertEqual(summary["background_region_block_count_rejections"], "1")
            self.assertEqual(summary["region_block_count_tolerance"], "1")
            self.assertEqual(summary["region_block_count_allow_swap"], "true")
            backgrounds = read_rows(directory / "tolerant.background_regions.tsv")
            self.assertEqual([(row["n_eligible_left_blocks"], row["n_eligible_right_blocks"]) for row in backgrounds],
                             [("3", "3"), ("3", "2"), ("1", "3")])
            self.assertEqual([row["block_count_match_orientation"] for row in backgrounds], ["direct", "swapped", "direct"])
            self.assertEqual([(row["matched_left_block_count_delta"], row["matched_right_block_count_delta"])
                              for row in backgrounds], [("1", "0"), ("0", "0"), ("-1", "0")])
            self.assertEqual([row["region_start"] for row in backgrounds], ["35", "40", "45"])
            self.assertIn("L3-R2", backgrounds[1]["background_valid_pair_labels"])
            self.assertNotIn("R3", backgrounds[1]["background_valid_pair_labels"])
            metadata = json.loads((directory / "tolerant.run.json").read_text())
            self.assertEqual(metadata["region_block_count_tolerance"], 1)
            self.assertTrue(metadata["region_block_count_allow_swap"])
            self.assertEqual(metadata["region_sampling_counts"]["target"]["accepted_block_count_match_orientations"],
                             {"direct": 2, "swapped": 1})
            self.assertEqual(metadata["region_sampling_counts"]["target"]["accepted_available_block_count_counts"],
                             {"1,3": 1, "3,2": 1, "3,3": 1})
            self.assertIn("block_count_tolerance=1 block_count_allow_swap=true", (directory / "tolerant.log").read_text())

    def test_nearest_fallback_fills_shortfall_and_keeps_left_right_symmetry(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            skip_tests.BackgroundRegionPairPTests.write_inputs(directory)
            args = existing_tests.BackgroundBlockLDTests.args(directory, "fallback")
            args.block_count = 3
            args.max_flank_bp = 30
            args.background_regions = 2
            args.background_region_attempts = 3
            args.region_block_count_tolerance = 0
            args.region_block_count_allow_swap = True
            args.region_match_fallback = "nearest"
            args.region_fallback_max_score = 1
            args.skip_background_region_pair_p = True
            original = workflow.BlockStore.summarize
            original_eligibility = workflow.BlockStore.pair_is_eligible_without_dprime

            def controlled_rows(store, region, first, second):
                row = original(store, region, first, second)
                if region.region_id == "target" and first.rank == 3:
                    row["status"] = "insufficient_complete_haplotypes"  # target = (2, 3)
                elif "|bg_" in region.region_id:
                    if region.start == 30 and second.rank == 3:
                        row["status"] = "insufficient_complete_haplotypes"  # (3, 2), exact after swap
                    elif region.start == 35 and first.rank != 1:
                        row["status"] = "insufficient_complete_haplotypes"  # (1, 3), score 1
                    # start 40 remains (3, 3), also score 1; stable sampling order breaks the tie.
                return row

            def controlled_eligibility(store, region, first, second):
                if not original_eligibility(store, region, first, second):
                    return False
                if "|bg_" not in region.region_id:
                    return True
                if region.start == 30 and second.rank == 3:
                    return False
                if region.start == 35 and first.rank != 1:
                    return False
                return True

            def candidates(domains, limit, rng):
                self.assertEqual(limit, 3)
                for start in (30, 35, 40):
                    yield domains[0], start, 100

            with mock.patch.object(workflow.BlockStore, "summarize", controlled_rows), \
                    mock.patch.object(workflow.BlockStore, "pair_is_eligible_without_dprime",
                                      controlled_eligibility), \
                    mock.patch.object(workflow, "random_region_intervals", candidates), \
                    contextlib.redirect_stderr(io.StringIO()):
                workflow.run(args)

            read_rows = existing_tests.BackgroundBlockLDTests.read_rows
            summary = read_rows(directory / "fallback.summary.tsv")[0]
            self.assertEqual(summary["background_region_attempts"], "3")
            self.assertEqual(summary["region_strict_match_count"], "1")
            self.assertEqual(summary["region_fallback_candidate_count"], "2")
            self.assertEqual(summary["region_fallback_selected_count"], "1")
            self.assertEqual(float(summary["selected_region_match_score_max"]), 1)

            backgrounds = read_rows(directory / "fallback.background_regions.tsv")
            self.assertEqual([row["region_start"] for row in backgrounds], ["30", "35"])
            self.assertEqual([row["region_match_selection"] for row in backgrounds],
                             ["strict", "nearest_fallback"])
            self.assertEqual([row["region_match_strict"] for row in backgrounds], ["true", "false"])
            self.assertEqual([row["block_count_match_orientation"] for row in backgrounds],
                             ["swapped", "direct"])
            self.assertEqual([float(row["region_match_score"]) for row in backgrounds], [0, 1])
            metadata = json.loads((directory / "fallback.run.json").read_text())
            info = metadata["region_sampling_counts"]["target"]
            self.assertEqual(info["fallback_selected_count"], 1)
            self.assertEqual(info["attempt_limit"], 3)
            log = (directory / "fallback.log").read_text()
            self.assertIn("REGION FALLBACK target=target", log)
            self.assertIn("fallback_selected=1", log)

    def test_length_bounds_cli_default_zero_and_negative(self):
        self.assertEqual(workflow.region_length_bounds(100_000, 10_000), (90_000, 110_000))
        self.assertEqual(workflow.region_length_bounds(100, 0), (100, 100))
        self.assertEqual(workflow.region_length_bounds(5, 10), (1, 15))
        argv = ["target.bed", "background.bed", "cohort.vcf", "-o", "output"]
        self.assertEqual(workflow.build_parser().parse_args(argv).region_length_tolerance_bp, 10_000)
        self.assertEqual(workflow.build_parser().parse_args(
            argv + ["--region-length-tolerance-bp", "0"]).region_length_tolerance_bp, 0)
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as error:
            workflow.main(argv + ["--region-length-tolerance-bp", "-1"])
        self.assertEqual(error.exception.code, 2)

    def test_configuration_domains_and_decode_match_exhaustive_containment(self):
        for component_bp in range(3, 22):
            for central_bp in (1, 3, 8):
                for tolerance in (0, 1, 5):
                    with self.subTest(component_bp=component_bp, length=central_bp, tolerance=tolerance):
                        background = [Region("1", "1", 10, 10 + component_bp, "bed", 1)]
                        domains = workflow.region_length_sampling_domains(background, central_bp, 2, tolerance)
                        lower, upper = workflow.region_length_bounds(central_bp, tolerance)
                        expected = {(start, length)
                                    for start in range(10, 10 + component_bp)
                                    for length in range(lower, upper + 1)
                                    if start - 2 >= 10 and start + length + 2 <= 10 + component_bp}
                        decoded = [domain.decode(index) for domain in domains for index in range(domain.size)]
                        self.assertEqual(set(decoded), expected)
                        self.assertEqual(len(decoded), len(expected))

    def test_sampling_is_unique_deterministic_and_fixed_length_matches_old_sampler(self):
        background = [Region("1", "1", 0, 18, "first", 1),
                      Region("2", "2", 40, 55, "second", 2)]
        domains = workflow.region_length_sampling_domains(background, 4, 3, 2)
        first = [(domain.chrom, start, length) for domain, start, length in
                 workflow.random_region_intervals(domains, 10000, random.Random(41))]
        second = [(domain.chrom, start, length) for domain, start, length in
                  workflow.random_region_intervals(domains, 10000, random.Random(41))]
        self.assertEqual(first, second)
        self.assertEqual(len(first), sum(domain.size for domain in domains))
        self.assertEqual(len(first), len(set(first)))
        fixed = workflow.region_length_sampling_domains(background, 4, 3, 0)
        old = workflow.region_sampling_domains(background, 4, 3)
        actual = [(domain.chrom, start) for domain, start, _length in
                  workflow.random_region_intervals(fixed, 10000, random.Random(41))]
        expected = [(domain.chrom, start) for domain, start in
                    workflow.random_starts(old, 10000, random.Random(41))]
        self.assertEqual(actual, expected)

    def test_available_blocks_count_unique_valid_endpoints_not_all_windows(self):
        pair = existing_tests.BackgroundBlockLDTests.pair_row
        rows = [pair(1, 1, 0.5, 0.2), pair(1, 3, 0.6, 0.2),
                pair(4, 1, 0.5, 0.2, status="insufficient_complete_haplotypes")]
        self.assertEqual(workflow.available_block_counts(rows), (1, 2))
        workflow.check_available_block_counts(rows, (1, 2))
        self.assertEqual(workflow.check_available_block_counts(rows, (2, 1)), "swapped")
        with self.assertRaisesRegex(ValueError, "available block count matching"):
            workflow.check_available_block_counts(rows, (2, 1), tolerance=0, allow_swap=False)
        self.assertEqual(workflow.available_block_counts([]), (0, 0))

    def test_region_accepts_extra_pairs_and_uses_its_own_maximum_not_target_mask(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            skip_tests.BackgroundRegionPairPTests.write_inputs(directory)
            args = existing_tests.BackgroundBlockLDTests.args(directory, "own_pairs")
            args.background_regions = 2
            args.skip_background_region_pair_p = True
            original = workflow.BlockStore.summarize

            def controlled_rows(store, region, first, second):
                row = original(store, region, first, second)
                if region.region_id == "target":
                    row["multiallelic_dprime"] = "0.4"
                    if (first.rank, second.rank) == (2, 1):
                        row["status"] = "insufficient_complete_haplotypes"
                elif "|bg_" in region.region_id:
                    row["multiallelic_dprime"] = "1" if (first.rank, second.rank) == (2, 1) else "0.2"
                return row

            with mock.patch.object(workflow.BlockStore, "summarize", controlled_rows), \
                    contextlib.redirect_stderr(io.StringIO()):
                workflow.run(args)
            read_rows = existing_tests.BackgroundBlockLDTests.read_rows
            summary = read_rows(directory / "own_pairs.summary.tsv")[0]
            self.assertEqual(summary["n_valid_pairs"], "3")
            for row in read_rows(directory / "own_pairs.background_regions.tsv"):
                self.assertEqual((row["n_eligible_left_blocks"], row["n_eligible_right_blocks"]), ("2", "2"))
                self.assertEqual(row["n_compared_pairs"], "4")
                self.assertEqual(row["max_dprime"], "1")
                self.assertEqual((row["peak_block1_label"], row["peak_block2_label"]), ("L2", "R1"))
                self.assertIn("L2-R1", row["background_valid_pair_labels"])
                self.assertEqual(row["target_pair_mask"], "NA")

    def test_counts_match_each_side_separately_and_background_ranks_can_differ(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            skip_tests.BackgroundRegionPairPTests.write_inputs(directory)
            args = existing_tests.BackgroundBlockLDTests.args(directory, "counts")
            args.background_regions = 1
            args.log_every_attempts = 1
            original = workflow.BlockStore.summarize
            original_eligibility = workflow.BlockStore.pair_is_eligible_without_dprime

            def controlled_rows(store, region, first, second):
                row = original(store, region, first, second)
                if region.region_id == "target" and first.rank == 2:
                    row["status"] = "insufficient_complete_haplotypes"
                elif "|bg_" in region.region_id:
                    # Same total count (3) but reversed sides must be rejected.
                    if ((region.start == 20 and second.rank == 2) or
                            (region.start != 20 and first.rank == 1)):
                        row["status"] = "insufficient_complete_haplotypes"
                return row

            def controlled_eligibility(store, region, first, second):
                if not original_eligibility(store, region, first, second):
                    return False
                if "|bg_" not in region.region_id:
                    return True
                return not ((region.start == 20 and second.rank == 2) or
                            (region.start != 20 and first.rank == 1))

            def two_candidates(domains, limit, rng):
                yield domains[0], 20, 100
                yield domains[0], 25, 100

            with mock.patch.object(workflow.BlockStore, "summarize", controlled_rows), \
                    mock.patch.object(workflow.BlockStore, "pair_is_eligible_without_dprime",
                                      controlled_eligibility), \
                    mock.patch.object(workflow, "random_region_intervals", two_candidates), \
                    contextlib.redirect_stderr(io.StringIO()):
                workflow.run(args)
            read_rows = existing_tests.BackgroundBlockLDTests.read_rows
            summary = read_rows(directory / "counts.summary.tsv")[0]
            self.assertEqual(summary["background_region_attempts"], "2")
            self.assertEqual(summary["background_region_block_count_rejections"], "1")
            self.assertEqual(summary["background_region_qc_rejections"], "0")
            row = read_rows(directory / "counts.background_regions.tsv")[0]
            self.assertEqual((row["n_eligible_left_blocks"], row["n_eligible_right_blocks"]), ("1", "2"))
            self.assertEqual(row["background_valid_pair_labels"], "L2-R1,L2-R2")
            metadata = json.loads((directory / "counts.run.json").read_text())
            self.assertEqual(metadata["pair_sampling_counts"]["130:1:1"]["source"], "background_region_pairs")
            self.assertIn("block_count_rejected=1", (directory / "counts.log").read_text())

    def test_variable_lengths_exact_pair_references_and_skip_mode_preserve_region_sampling(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            skip_tests.BackgroundRegionPairPTests.write_inputs(directory)
            for prefix, skip, tolerance in (("full", False, 3), ("skip", True, 3), ("strict", True, 0)):
                args = existing_tests.BackgroundBlockLDTests.args(directory, prefix)
                args.background_regions = 12
                args.region_length_tolerance_bp = tolerance
                args.skip_background_region_pair_p = skip
                args.pair_significance_threshold = 0.3
                with contextlib.redirect_stderr(io.StringIO()):
                    workflow.run(args)
            self.assertEqual((directory / "full.background_regions.bed").read_text(),
                             (directory / "skip.background_regions.bed").read_text())
            for prefix in ("skip", "strict"):
                self.assertEqual((directory / "full.target_pairs.tsv").read_text(),
                                 (directory / f"{prefix}.target_pairs.tsv").read_text())
            read_rows = existing_tests.BackgroundBlockLDTests.read_rows
            backgrounds = read_rows(directory / "full.background_regions.tsv")
            self.assertGreater(len({row["central_length_bp"] for row in backgrounds}), 1)
            banks = {}
            for row in read_rows(directory / "full.background_pairs.tsv"):
                banks.setdefault(workflow.pair_background_key(row, True), []).append(float(row["multiallelic_dprime"]))
            data = read_vcf(directory / "cohort.vcf", [Region("1", "1", 0, 1000, "all", 1)],
                            0, 0.1, 0.95, False, False, include_region_interiors=True)
            store = workflow.BlockStore(data, args)
            for bg in backgrounds:
                start, end = int(bg["region_start"]), int(bg["region_end"])
                self.assertLessEqual(abs(end - start - 100), 3)
                self.assertGreaterEqual(start - 20, 0)
                self.assertLessEqual(end + 20, 1000)
                self.assertTrue(end + 20 <= 280 or start - 20 >= 420)
                self.assertEqual(bg["target_pair_mask"], "NA")
                rows = [row for row in store.region_pairs(Region("1", "1", start, end, "check", 0))
                        if workflow.valid_row(row)]
                self.assertEqual(workflow.available_block_counts(rows), (2, 2))
                p_values = []
                for row in rows:
                    bank = banks[workflow.pair_background_key(row, True)]
                    self.assertEqual(len(bank), 8)
                    p_values.append((1 + sum(value >= float(row["multiallelic_dprime"]) - 1e-12
                                             for value in bank)) / 9)
                self.assertAlmostEqual(float(bg["significant_pair_fraction"]),
                                       sum(p < 0.3 for p in p_values) / len(rows))
                peak_index = max(range(len(rows)), key=lambda index: float(rows[index]["multiallelic_dprime"]))
                self.assertAlmostEqual(float(bg["peak_pair_nominal_p"]), p_values[peak_index])
            full_summary = read_rows(directory / "full.summary.tsv")[0]
            skip_summary = read_rows(directory / "skip.summary.tsv")[0]
            self.assertEqual(full_summary["region_max_empirical_p"], skip_summary["region_max_empirical_p"])
            metadata = json.loads((directory / "full.run.json").read_text())
            self.assertFalse(metadata["region_rank_pair_mask_matching"])
            self.assertEqual(metadata["region_pair_distance_matching"], "none")
            self.assertEqual(sum(metadata["region_sampling_counts"]["target"]["accepted_central_length_counts"].values()), 12)
            self.assertTrue(any(info["source"] == "background_region_pairs"
                                for info in metadata["pair_sampling_counts"].values()))
            skip_meta = json.loads((directory / "skip.run.json").read_text())
            self.assertTrue(all(info["source"] == "target_pairs"
                                for info in skip_meta["pair_sampling_counts"].values()))
            self.assertIn("on_demand=true", (directory / "full.log").read_text())

    def test_missing_on_demand_bank_fails_run_instead_of_rejecting_region(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            skip_tests.BackgroundRegionPairPTests.write_inputs(directory)
            args = existing_tests.BackgroundBlockLDTests.args(directory, "bank_fail")
            args.background_regions = 1
            args.region_length_tolerance_bp = 1
            args.max_attempt_factor = 2
            original = workflow.BlockStore.summarize

            def reject_new_bank(store, region, first, second):
                row = original(store, region, first, second)
                if region.region_id.startswith("pair_d111_"):
                    row["status"] = "below_min_complete_haplotypes"
                return row

            def one_candidate(domains, limit, rng):
                yield domains[0], 20, 101

            with mock.patch.object(workflow.BlockStore, "summarize", reject_new_bank), \
                    mock.patch.object(workflow, "random_region_intervals", one_candidate), \
                    contextlib.redirect_stderr(io.StringIO()), self.assertRaisesRegex(ValueError, "distance 111"):
                workflow.run(args)
            self.assertFalse((directory / "bank_fail.run.json").exists())
            self.assertIn("RUN FAILED", (directory / "bank_fail.log").read_text())
            self.assertEqual(existing_tests.BackgroundBlockLDTests.read_rows(
                directory / "bank_fail.background_regions.tsv"), [])


if __name__ == "__main__":
    unittest.main()

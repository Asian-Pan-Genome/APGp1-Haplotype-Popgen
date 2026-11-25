import contextlib
import io
import json
import math
import random
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from flanking_block_ld import background_block_ld as workflow
from flanking_block_ld.flanking_block_ld import Region, read_vcf
from flanking_block_ld.tests import test_background_block_ld as existing_tests
from flanking_block_ld.tests import test_background_region_pair_p as skip_tests


class BackgroundRegionEntropyTests(unittest.TestCase):
    @staticmethod
    def rows_for(left, right):
        return [existing_tests.BackgroundBlockLDTests.pair_row(i, j, 0.5, 0.2)
                for i in range(1, left + 1) for j in range(1, right + 1)]

    def test_cli_default_override_opt_out_and_invalid_thresholds(self):
        argv = ["target.bed", "background.bed", "cohort.vcf", "-o", "output"]
        args = workflow.build_parser().parse_args(argv)
        self.assertEqual(args.region_entropy_max_wasserstein, 0.05)
        self.assertFalse(args.skip_region_entropy_matching)
        changed = workflow.build_parser().parse_args(
            argv + ["--region-entropy-max-wasserstein", "0.1", "--skip-region-entropy-matching"])
        self.assertEqual(changed.region_entropy_max_wasserstein, 0.1)
        self.assertTrue(changed.skip_region_entropy_matching)
        for value in ("-0.01", "nan", "inf", "-inf"):
            with self.subTest(value=value), contextlib.redirect_stderr(io.StringIO()), \
                    self.assertRaises(SystemExit) as error:
                workflow.main(argv + ["--region-entropy-max-wasserstein", value])
            self.assertEqual(error.exception.code, 2)

    def test_w1_sorted_equal_counts_unequal_counts_and_equal_block_weights(self):
        distance = workflow.wasserstein_distance
        self.assertAlmostEqual(distance([1.5, 0.5, 1.0], [0.6, 1.2, 1.8]), 0.2)
        self.assertAlmostEqual(distance([0.5, 1.0, 1.5], [0.5, 1.5]), 1 / 6)
        self.assertEqual(distance([0.5, 1.5], [0.5, 0.5, 1.5, 1.5]), 0)
        self.assertEqual(distance([1, 1], [1]), 0)
        self.assertEqual(distance([1], [2, 2]), 1)

    def test_w1_matches_independent_quantile_integral_reference(self):
        rng = random.Random(182)
        for n in range(1, 7):
            for m in range(1, 7):
                for _ in range(5):
                    a = sorted(rng.randint(0, 10) / 3 for _ in range(n))
                    b = sorted(rng.randint(0, 10) / 3 for _ in range(m))
                    expected = math.fsum(abs(a[j // m] - b[j // n])
                                         for j in range(n * m)) / (n * m)
                    self.assertAlmostEqual(workflow.wasserstein_distance(a, b), expected)
                    self.assertAlmostEqual(workflow.wasserstein_distance(b, a), expected)

    def test_w1_rejects_empty_or_nonfinite_samples(self):
        for first, second in (([], [1]), ([1], []), ([math.nan], [1]), ([1], [math.inf])):
            with self.assertRaisesRegex(ValueError, "nonempty finite"):
                workflow.wasserstein_distance(first, second)

    def test_swap_is_chosen_by_entropy_not_count_only_preselection(self):
        rows = self.rows_for(2, 2)
        target = workflow.RegionEntropyProfile((0.2, 0.4), (1.0, 1.2))
        background = workflow.RegionEntropyProfile((1.2, 1.0), (0.4, 0.2))
        self.assertEqual(workflow.check_available_block_counts(rows, (2, 2), 0), "direct")
        match = workflow.match_region_entropy(rows, (2, 2), background, target, tolerance=0)
        self.assertEqual(match.orientation, "swapped")
        self.assertEqual(match.max_wasserstein, 0)
        self.assertGreater(match.direct_max_wasserstein, 0.05)
        with self.assertRaisesRegex(ValueError, "entropy distribution matching"):
            workflow.match_region_entropy(rows, (2, 2), background, target, 0, False)
        disabled = workflow.match_region_entropy(rows, (2, 2), background, target, enabled=False)
        self.assertEqual(disabled.orientation, "direct")
        self.assertGreater(disabled.max_wasserstein, 0.05)

    def test_joint_candidate_score_is_left_right_symmetric_and_uses_one_orientation(self):
        rows = self.rows_for(9, 8)
        target = workflow.RegionEntropyProfile((0.2,) * 8, (1.0,) * 9)

        # Counts favor an exact swap, but entropy only passes in the direct
        # orientation.  The scorer must not mix the count and entropy mappings.
        direct_entropy = workflow.RegionEntropyProfile((0.2,) * 9, (1.0,) * 8)
        direct = workflow.score_region_candidate(
            rows, (8, 9), direct_entropy, target,
            central_length_delta=0, length_tolerance=10_000,
            count_tolerance=1, allow_swap=True,
            max_wasserstein=0.05, entropy_enabled=True,
        )
        self.assertTrue(direct.strict)
        self.assertEqual(direct.orientation, "direct")
        self.assertEqual((direct.left_count_error, direct.right_count_error), (1, 1))

        # Reversing the two entropy distributions makes (left 9, right 8)
        # equivalent to target (left 8, right 9) under one global swap.
        swapped_entropy = workflow.RegionEntropyProfile((1.0,) * 9, (0.2,) * 8)
        swapped = workflow.score_region_candidate(
            rows, (8, 9), swapped_entropy, target,
            central_length_delta=0, length_tolerance=10_000,
            count_tolerance=0, allow_swap=True,
            max_wasserstein=0.05, entropy_enabled=True,
        )
        self.assertTrue(swapped.strict)
        self.assertEqual(swapped.orientation, "swapped")
        self.assertEqual(swapped.score, 0)
        self.assertEqual((swapped.left_count_error, swapped.right_count_error), (0, 0))

    def test_entropy_can_choose_less_exact_count_alignment_within_tolerance(self):
        rows = self.rows_for(3, 2)
        target = workflow.RegionEntropyProfile((0.2, 0.2), (1.0, 1.0, 1.0))
        background = workflow.RegionEntropyProfile((0.2, 0.2, 0.2), (1.0, 1.0))
        self.assertEqual(workflow.check_available_block_counts(rows, (2, 3)), "swapped")
        match = workflow.match_region_entropy(rows, (2, 3), background, target)
        self.assertEqual(match.orientation, "direct")
        self.assertEqual(match.max_wasserstein, 0)
        # Direct passes entropy but not exact counts; swapped passes exact
        # counts but not entropy. No single orientation passes both.
        with self.assertRaisesRegex(ValueError, "entropy distribution matching"):
            workflow.match_region_entropy(rows, (2, 3), background, target, tolerance=0)

    def test_per_side_matching_cannot_be_replaced_by_pooled_entropy(self):
        target = workflow.RegionEntropyProfile((0.2, 0.2), (1.0, 1.0))
        background = workflow.RegionEntropyProfile((0.2, 1.0), (0.2, 1.0))
        self.assertEqual(workflow.wasserstein_distance(target.left + target.right,
                                                      background.left + background.right), 0)
        with self.assertRaisesRegex(ValueError, "entropy distribution matching"):
            workflow.match_region_entropy(self.rows_for(2, 2), (2, 2), background, target)

    def test_threshold_boundary_and_profile_validation(self):
        rows = self.rows_for(1, 1)
        target = workflow.RegionEntropyProfile((0.1,), (0.1,))
        workflow.match_region_entropy(rows, (1, 1), workflow.RegionEntropyProfile((0.15,), (0.15,)), target)
        workflow.match_region_entropy(rows, (1, 1), target, target, max_wasserstein=0)
        with self.assertRaisesRegex(ValueError, "entropy distribution matching"):
            workflow.match_region_entropy(rows, (1, 1), workflow.RegionEntropyProfile((0.150001,), (0.1,)), target)
        with self.assertRaisesRegex(ValueError, "available blocks"):
            workflow.match_region_entropy(rows, (1, 1), workflow.RegionEntropyProfile((0.1, 0.1), (0.1,)), target)
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            workflow.match_region_entropy(rows, (1, 1), workflow.RegionEntropyProfile((-0.1,), (0.1,)), target)

    def test_block_entropy_is_independent_of_partner_calls_and_counted_once(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            vcf = directory / "cohort.vcf"
            vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n"
                "1\t95\tleft\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|1\t1|1\t1|1\n"
                "1\t205\tright1\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|1\t0|1\t.|.\n"
                "1\t215\tright2\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|1\t0|1\t0|1\n",
                encoding="utf-8")
            region = Region("1", "1", 100, 200, "test", 1)
            data = read_vcf(vcf, [region], flank_bp=20, min_maf=0.1, min_call_rate=0.5,
                            include_filtered=False, allow_unphased_missing=False)
            store = workflow.BlockStore(data, existing_tests.BackgroundBlockLDTests.args(directory, "test"))
            rows = store.region_pairs(region)
            valid = [row for row in rows if workflow.valid_row(row)]
            self.assertEqual([row["complete_haplotype_count"] for row in valid], ["6", "8"])
            self.assertNotEqual(valid[0]["left_haplotype_entropy"], valid[1]["left_haplotype_entropy"])
            with mock.patch.object(store, "window_entropy", wraps=store.window_entropy) as calls:
                profile = store.region_entropy_profile(region, rows)
                self.assertEqual(calls.call_count, 3)  # One left + two right blocks.
            expected = -(3 / 8) * math.log(3 / 8) - (5 / 8) * math.log(5 / 8)
            self.assertEqual(profile.left_complete_haplotypes, (8,))
            self.assertEqual(profile.right_complete_haplotypes, (6, 8))
            self.assertAlmostEqual(profile.left[0], expected)
            self.assertEqual(profile.left_ranks, (1,))
            self.assertEqual(profile.right_ranks, (1, 2))

    def test_candidate_screen_skips_dprime_and_matches_full_pair_validity(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            skip_tests.BackgroundRegionPairPTests.write_inputs(directory)
            args = existing_tests.BackgroundBlockLDTests.args(directory, "screen")
            data = read_vcf(
                directory / "cohort.vcf",
                [Region("1", "1", 0, 1000, "all", 1)],
                flank_bp=0, min_maf=0.1, min_call_rate=0.95,
                include_filtered=False, allow_unphased_missing=False,
                include_region_interiors=True,
            )
            store = workflow.BlockStore(data, args)
            region = Region("1", "1", 100, 200, "candidate", 1)
            args.skip_region_entropy_matching = False

            with mock.patch.object(
                    workflow, "summarize_block_pair",
                    side_effect=AssertionError("screen must not calculate pair Dprime/entropy")):
                screen = store.screen_region(region)

            rows = store.region_pairs_from_blocks(
                region, screen.left_blocks, screen.right_blocks,
            )
            self.assertEqual(screen.eligible_pairs, workflow.region_valid_pairs(rows))
            self.assertEqual(screen.available_counts, workflow.available_block_counts(rows))
            self.assertEqual(screen.entropy, store.region_entropy_profile(region, rows))

            # The complete-haplotype threshold remains part of screening even
            # though the numerical D' value is deferred.
            args.min_complete_haplotypes = data.sample_count * 2 + 1
            filtered = store.screen_region(region)
            self.assertEqual(filtered.eligible_pairs, ())
            self.assertEqual(filtered.available_counts, (0, 0))

    def test_skip_entropy_matching_avoids_target_and_candidate_region_entropy(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            skip_tests.BackgroundRegionPairPTests.write_inputs(directory)
            args = existing_tests.BackgroundBlockLDTests.args(directory, "skip_entropy_work")
            args.skip_region_entropy_matching = True
            with mock.patch.object(
                    workflow.BlockStore, "region_entropy_profile",
                    side_effect=AssertionError("region entropy must be skipped")), \
                    contextlib.redirect_stderr(io.StringIO()):
                workflow.run(args)

            read_rows = existing_tests.BackgroundBlockLDTests.read_rows
            for row in read_rows(directory / "skip_entropy_work.background_regions.tsv"):
                self.assertEqual(row["region_entropy_wasserstein"], "NA")
                self.assertEqual(row["left_block_haplotype_entropies"], "NA")
                self.assertEqual(row["right_block_haplotype_entropies"], "NA")
            metadata = json.loads((directory / "skip_entropy_work.run.json").read_text())
            info = metadata["region_sampling_counts"]["target"]
            self.assertIsNone(info["target_entropy_profile"])
            self.assertIsNone(info["accepted_entropy_wasserstein_summary"])
            self.assertIn(
                "marginal_entropy=skipped",
                (directory / "skip_entropy_work.log").read_text(),
            )

    def test_skip_entropy_candidate_screen_does_not_select_snp_subsets(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            skip_tests.BackgroundRegionPairPTests.write_inputs(directory)
            args = existing_tests.BackgroundBlockLDTests.args(directory, "deferred")
            args.skip_region_entropy_matching = True
            data = read_vcf(
                directory / "cohort.vcf",
                [Region("1", "1", 0, 1000, "all", 1)],
                flank_bp=0, min_maf=0.1, min_call_rate=0.95,
                include_filtered=False, allow_unphased_missing=False,
                include_region_interiors=True,
            )
            store = workflow.BlockStore(data, args)
            with mock.patch.object(
                    workflow, "cap_variants_by_entropy",
                    side_effect=AssertionError("candidate screen must defer SNP entropy selection")):
                screen = store.screen_region(Region("1", "1", 100, 200, "candidate", 1))
            self.assertTrue(screen.eligible_pairs)
            self.assertIsNone(screen.entropy)

    def test_workflow_defaults_apply_matching_without_changing_pair_background(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            skip_tests.BackgroundRegionPairPTests.write_inputs(directory)
            for prefix, entropy_enabled in (("enabled", True), ("disabled", False)):
                args = existing_tests.BackgroundBlockLDTests.args(directory, prefix)
                if entropy_enabled:
                    del args.region_entropy_max_wasserstein
                    del args.skip_region_entropy_matching
                with contextlib.redirect_stderr(io.StringIO()):
                    workflow.run(args)
            read_rows = existing_tests.BackgroundBlockLDTests.read_rows
            summary = read_rows(directory / "enabled.summary.tsv")[0]
            self.assertEqual(summary["region_entropy_matching_enabled"], "true")
            self.assertEqual(float(summary["region_entropy_max_wasserstein"]), 0.05)
            for row in read_rows(directory / "enabled.background_regions.tsv"):
                self.assertLessEqual(float(row["region_entropy_wasserstein"]), 0.05 + 1e-12)
                self.assertEqual(len(json.loads(row["left_block_haplotype_entropies"])), 2)
            for suffix in ("target_pairs.tsv", "background_pairs.tsv", "background_regions.bed"):
                self.assertEqual((directory / f"enabled.{suffix}").read_text(),
                                 (directory / f"disabled.{suffix}").read_text())

    def test_workflow_entropy_rejection_orientation_diagnostics_and_skip_pair_p(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            skip_tests.BackgroundRegionPairPTests.write_inputs(directory)
            original_profile = workflow.BlockStore.region_entropy_profile
            original_rows = workflow.BlockStore.region_pairs
            original_eligibility = workflow.BlockStore.pair_is_eligible_without_dprime

            def profiles(store, region, rows):
                original = original_profile(store, region, rows)
                if "|bg_" not in region.region_id or region.start == 30:
                    values = ((0.2, 0.4), (1.0, 1.2))
                elif region.start == 25:
                    values = ((1.0, 1.2), (0.4, 0.2))
                else:
                    values = ((0.6, 0.6), (1.0, 1.2))
                return workflow.replace(original, left=values[0], right=values[1])

            def rows_for(store, region, required_snp_counts=None):
                rows = original_rows(store, region, required_snp_counts)
                if "|bg_" in region.region_id and region.start == 35:
                    for row in rows:
                        if workflow.row_key(row) != (1, 1):
                            row["status"] = "insufficient_complete_haplotypes"
                return rows

            def eligible_without_dprime(store, region, first, second):
                if not original_eligibility(store, region, first, second):
                    return False
                return not ("|bg_" in region.region_id and region.start == 35 and
                            (first.rank, second.rank) != (1, 1))

            def candidates(domains, limit, rng):
                for start in (20, 35, 25, 30):
                    yield domains[0], start, 100

            for prefix, skip_pair_p in (("full", False), ("skip", True)):
                args = existing_tests.BackgroundBlockLDTests.args(directory, prefix)
                args.background_regions = 2
                args.skip_region_entropy_matching = False
                args.region_block_count_allow_swap = True
                args.skip_background_region_pair_p = skip_pair_p
                args.log_every_attempts = 1
                with mock.patch.object(workflow.BlockStore, "region_entropy_profile", profiles), \
                        mock.patch.object(workflow.BlockStore, "region_pairs", rows_for), \
                        mock.patch.object(workflow.BlockStore, "pair_is_eligible_without_dprime",
                                          eligible_without_dprime), \
                        mock.patch.object(workflow, "random_region_intervals", candidates), \
                        contextlib.redirect_stderr(io.StringIO()):
                    workflow.run(args)
                read_rows = existing_tests.BackgroundBlockLDTests.read_rows
                summary = read_rows(directory / f"{prefix}.summary.tsv")[0]
                self.assertEqual(summary["background_region_attempts"], "4")
                self.assertEqual(summary["background_region_entropy_rejections"], "1")
                self.assertEqual(summary["background_region_block_count_rejections"], "1")
                self.assertEqual(summary["background_region_qc_rejections"], "0")
                backgrounds = read_rows(directory / f"{prefix}.background_regions.tsv")
                self.assertEqual([row["region_start"] for row in backgrounds], ["25", "30"])
                self.assertEqual([row["block_count_match_orientation"] for row in backgrounds], ["swapped", "direct"])
                self.assertEqual(float(backgrounds[0]["swapped_entropy_max_wasserstein"]), 0)
                self.assertGreater(float(backgrounds[0]["direct_entropy_max_wasserstein"]), 0.05)
                self.assertEqual(json.loads(backgrounds[0]["left_block_haplotype_entropies"]), [1.0, 1.2])
                metadata = json.loads((directory / f"{prefix}.run.json").read_text())
                info = metadata["region_sampling_counts"]["target"]
                self.assertEqual(info["entropy_rejections"], 1)
                self.assertEqual(info["accepted_entropy_wasserstein_summary"]["max"], 0)
                self.assertIn("entropy_rejected=1", (directory / f"{prefix}.log").read_text())
            self.assertEqual((directory / "full.background_regions.bed").read_text(),
                             (directory / "skip.background_regions.bed").read_text())

    def test_entropy_exhaustion_is_logged_and_never_relaxes_threshold(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            skip_tests.BackgroundRegionPairPTests.write_inputs(directory)
            args = existing_tests.BackgroundBlockLDTests.args(directory, "fail")
            args.skip_region_entropy_matching = False
            args.background_regions = 2
            args.max_attempt_factor = 2
            args.log_every_attempts = 1
            original = workflow.BlockStore.region_entropy_profile

            def reject_background(store, region, rows):
                profile = original(store, region, rows)
                if "|bg_" in region.region_id:
                    return workflow.replace(profile, left=tuple(value + 0.1 for value in profile.left))
                return profile

            with mock.patch.object(workflow.BlockStore, "region_entropy_profile", reject_background), \
                    contextlib.redirect_stderr(io.StringIO()), \
                    self.assertRaisesRegex(ValueError, "never relaxed automatically"):
                workflow.run(args)
            log = (directory / "fail.log").read_text()
            self.assertIn("EXHAUSTED accepted=0/2", log)
            self.assertIn("entropy_rejected=4", log)
            self.assertIn("[ERROR] RUN FAILED", log)
            self.assertFalse((directory / "fail.run.json").exists())


if __name__ == "__main__":
    unittest.main()

import argparse
import contextlib
import io
import json
import random
import tempfile
import unittest
from pathlib import Path

from flanking_block_ld import background_block_ld as workflow
from flanking_block_ld.background_block_ld import (
    BlockStore,
    analyze_region,
    attach_pair_p,
    build_parser,
    empirical_tail_p,
    main,
    pair_background_key,
    region_snp_counts,
    resolve_snp_limits,
    run,
)
from flanking_block_ld.flanking_block_ld import Region, Variant, VCFData
from flanking_block_ld.tests import test_background_block_ld as existing_tests


class BackgroundSNPRangeTests(unittest.TestCase):
    @staticmethod
    def limits(**overrides):
        values = dict(snps_per_block=None, min_snps_per_block=None, max_snps_per_block=None)
        values.update(overrides)
        return argparse.Namespace(**values)

    @staticmethod
    def pair_row(left_rank=1, right_rank=1, left_n=1, right_n=2, dprime=0.8):
        row = existing_tests.BackgroundBlockLDTests.pair_row(left_rank, right_rank, dprime, 0.005)
        row.update(
            block1_snp_number_pass_threshold=str(left_n),
            block2_snp_number_pass_threshold=str(right_n),
        )
        return row

    def test_snp_limits_preserve_default_and_legacy_fixed_k(self):
        self.assertEqual(resolve_snp_limits(self.limits()), (10, 10))
        self.assertEqual(resolve_snp_limits(argparse.Namespace(snps_per_block=4)), (4, 4))
        self.assertEqual(resolve_snp_limits(self.limits(min_snps_per_block=3, max_snps_per_block=5)), (3, 5))
        self.assertEqual(resolve_snp_limits(self.limits(min_snps_per_block=3)), (3, 10))
        self.assertEqual(resolve_snp_limits(self.limits(max_snps_per_block=12)), (10, 12))

    def test_snp_limits_reject_conflicting_and_invalid_arguments(self):
        for values in (
            dict(snps_per_block=3, min_snps_per_block=3),
            dict(snps_per_block=3, max_snps_per_block=3),
            dict(snps_per_block=0),
            dict(min_snps_per_block=0, max_snps_per_block=3),
            dict(min_snps_per_block=5, max_snps_per_block=3),
            dict(max_snps_per_block=5),
        ):
            with self.subTest(values=values), self.assertRaises(ValueError):
                resolve_snp_limits(self.limits(**values))

    def test_cli_range_and_per_stratum_alias(self):
        args = build_parser().parse_args([
            "target.bed", "background.bed", "cohort.vcf", "-o", "results",
            "--min-snps-per-block", "3", "--max-snps-per-block", "5",
            "--background-pairs-per-stratum", "20",
        ])
        self.assertIsNone(args.snps_per_block)
        self.assertEqual(resolve_snp_limits(args), (3, 5))
        self.assertEqual(args.background_pairs_per_distance, 20)
        self.assertEqual(args.min_pair_complete_haplotype_fraction, 0.20)
        changed = build_parser().parse_args([
            "target.bed", "background.bed", "cohort.vcf", "-o", "results",
            "--min-pair-complete-haplotype-fraction", "0.35",
        ])
        self.assertEqual(changed.min_pair_complete_haplotype_fraction, 0.35)
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as error:
            main([
                "target.bed", "background.bed", "cohort.vcf", "-o", "results",
                "--snps-per-block", "3", "--min-snps-per-block", "3",
            ])
        self.assertEqual(error.exception.code, 2)
        for value in ("-0.01", "1.01", "nan", "inf"):
            with self.subTest(value=value), contextlib.redirect_stderr(io.StringIO()), \
                    self.assertRaises(SystemExit) as error:
                main([
                    "target.bed", "background.bed", "cohort.vcf", "-o", "results",
                    "--min-pair-complete-haplotype-fraction", value,
                ])
            self.assertEqual(error.exception.code, 2)

    def test_pair_p_background_is_conditioned_on_ordered_selected_counts(self):
        rows = [
            self.pair_row(left_n=1, right_n=2),
            self.pair_row(left_n=2, right_n=1),
            self.pair_row(left_n=2, right_n=2),
        ]
        distance = int(float(rows[0]["block_center_distance_bp"]))
        self.assertEqual(pair_background_key(rows[0], True), (distance, 1, 2))
        backgrounds = {(distance, 1, 2): [0.1, 0.2], (distance, 2, 1): [0.9, 1.0]}
        attach_pair_p(rows, backgrounds, 0.5, match_snp_counts=True, distance_tolerance=0)
        self.assertAlmostEqual(float(rows[0]["pair_nominal_p"]), 1 / 3)
        self.assertEqual(rows[0]["pair_significant"], "true")
        self.assertEqual(rows[0]["pair_snp_count_stratum"], "1,2")
        self.assertEqual(float(rows[1]["pair_nominal_p"]), 1.0)
        self.assertEqual(rows[1]["pair_significant"], "false")
        self.assertEqual(rows[2]["pair_nominal_p"], "NA")
        self.assertEqual(rows[2]["pair_distance_background_n"], "0")

    def test_region_exact_count_matching_applies_only_to_eligible_blocks(self):
        rows = [
            self.pair_row(1, 1, 1, 2),
            self.pair_row(1, 2, 1, 1),
            self.pair_row(2, 1, 3, 2),
        ]
        mask = [(1, 1), (1, 2)]
        required = {("left", 1): 1, ("right", 1): 2, ("right", 2): 1}
        self.assertEqual(region_snp_counts(rows, mask), required)
        summary = analyze_region(rows, mask, 0.01, required)
        self.assertEqual(summary["n_valid_pairs"], 2)
        self.assertEqual(summary["n_eligible_left_blocks"], 1)
        rows[0]["block2_snp_number_pass_threshold"] = "1"
        with self.assertRaisesRegex(ValueError, "SNP count matching"):
            analyze_region(rows, mask, 0.01, required)
        # A partial count map must not silently leave an eligible endpoint unmatched.
        with self.assertRaisesRegex(ValueError, "SNP count matching"):
            analyze_region(rows, mask, 0.01, {("left", 1): 1})

    def test_block_store_keeps_range_counts_caps_max_and_filters_below_min(self):
        positions = (2, 12, 14, 22, 24, 26, 28)
        alleles = {(sample, copy): copy for sample in range(4) for copy in range(2)}
        variants = tuple(
            Variant("1", pos + 1, f"v{pos}", "A", "G", 0.5, 1.0, alleles, index)
            for index, pos in enumerate(positions, 1)
        )
        data = VCFData(4, {"1": variants}, {"1": positions}, {}, len(variants), 0, 0, 0, 0, 0)
        args = existing_tests.BackgroundBlockLDTests.args(Path("."), "unused")
        args.snps_per_block = None
        args.min_snps_per_block, args.max_snps_per_block = 2, 3
        store = BlockStore(data, args)
        sparse = store.window("1", 0, "left")
        within = store.window("1", 10, "right")
        dense = store.window("1", 20, "left")
        self.assertEqual(len(sparse.variants), 1)
        self.assertEqual(len(within.variants), 2)
        self.assertEqual(within.selection_method, "deferred_pair_specific_fixed_cohort")
        self.assertEqual(dense.snp_count_before_cap, 4)
        self.assertEqual(len(dense.variants), 4)
        self.assertEqual(dense.selection_method, "deferred_pair_specific_fixed_cohort")
        region = Region("1", "chr1", 10, 20, "target", 1)
        failed = store.summarize(region, sparse, within)
        self.assertEqual(failed["status"], "insufficient_left_snps")
        self.assertEqual(failed["block1_snp_number_pass_threshold"], "1")
        passed = store.summarize(region, dense, within)
        self.assertEqual(passed["status"], "ok")
        self.assertEqual(passed["block1_snp_number_pass_threshold"], "3")
        self.assertEqual(passed["block2_snp_number_pass_threshold"], "2")
        self.assertEqual(passed["left_snp_selection_method"], "exact_max_entropy")
        self.assertEqual(passed["right_snp_selection_method"], "all_within_upper_limit")

    def test_pair_specific_selection_uses_one_fixed_complete_haplotype_cohort(self):
        keys = [(f"S{index // 2}", index % 2) for index in range(10)]

        def variant(pos, name, values, keep=10):
            return Variant(
                "1", pos, name, "A", "G", 0.5, 1.0,
                {key: value for key, value in zip(keys[:keep], values)}, pos,
            )

        # Marginally A has higher entropy, but it is monomorphic on the six
        # haplotypes retained by the partner. Fixed-cohort selection must pick B.
        left_a = variant(1, "A", [0, 0, 0, 0, 0, 0, 1, 1, 1, 1])
        left_b = variant(2, "B", [0, 0, 0, 1, 1, 1, 0, 0, 0, 0])
        right = variant(11, "R", [0, 1, 0, 1, 0, 1], keep=6)
        data = VCFData(5, {}, {}, {}, 3, 0, 0, 0, 0, 0)
        args = existing_tests.BackgroundBlockLDTests.args(Path("."), "unused")
        args.snps_per_block = 1
        store = BlockStore(data, args)
        left = workflow.FixedBlock("left", 1, 0, 10, 2, 0, 0.0, True, 2,
                                   "deferred_pair_specific_fixed_cohort", (left_a, left_b))
        right_block = workflow.FixedBlock("right", 1, 10, 20, 1, 0, 0.0, True, 1,
                                          "deferred_pair_specific_fixed_cohort", (right,))
        region = Region("1", "1", 10, 10, "pair", 1)
        self.assertTrue(store.pair_is_eligible_without_dprime(region, left, right_block))
        row = store.summarize(region, left, right_block)
        self.assertEqual(row["status"], "ok")
        self.assertEqual(row["left_snps"], "B@2")
        self.assertEqual(row["complete_haplotype_count"], "6")
        self.assertAlmostEqual(float(row["complete_haplotype_fraction"]), 0.6)
        self.assertEqual(
            row["pair_snp_selection_mode"],
            "pair_specific_fixed_complete_cohort_deferred_to_dprime",
        )

    def test_pair_complete_haplotype_fraction_is_enforced_in_screen_and_dprime(self):
        keys = [(f"S{index // 2}", index % 2) for index in range(20)]
        common = keys[:3]
        left_variant = Variant("1", 1, "L", "A", "G", 0.5, 1.0,
                               {key: index % 2 for index, key in enumerate(common)}, 1)
        right_variant = Variant("1", 11, "R", "A", "G", 0.5, 1.0,
                                {key: (index + 1) % 2 for index, key in enumerate(common)}, 2)
        data = VCFData(10, {}, {}, {}, 2, 0, 0, 0, 0, 0)
        args = existing_tests.BackgroundBlockLDTests.args(Path("."), "unused")
        args.min_pair_complete_haplotype_fraction = 0.20
        store = BlockStore(data, args)
        left = workflow.FixedBlock("left", 1, 0, 10, 1, 0, 0.0, True, 1,
                                   "deferred_pair_specific_fixed_cohort", (left_variant,))
        right = workflow.FixedBlock("right", 1, 10, 20, 1, 0, 0.0, True, 1,
                                    "deferred_pair_specific_fixed_cohort", (right_variant,))
        region = Region("1", "1", 10, 10, "pair", 1)
        self.assertFalse(store.pair_is_eligible_without_dprime(region, left, right))
        row = store.summarize(region, left, right)
        self.assertEqual(row["status"], "below_min_complete_haplotype_fraction")
        self.assertEqual(row["required_complete_haplotype_count"], "4")
        self.assertEqual(row["multiallelic_dprime"], "NA")

        args.min_pair_complete_haplotype_fraction = 0.10
        permissive = BlockStore(data, args)
        self.assertTrue(permissive.pair_is_eligible_without_dprime(region, left, right))
        self.assertEqual(permissive.summarize(region, left, right)["status"], "ok")

    def test_range_workflow_stratifies_pair_background_but_only_matches_available_region_blocks(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            (directory / "target.bed").write_text("1\t300\t400\ttarget\n", encoding="utf-8")
            (directory / "background.bed").write_text("1\t0\t1000\tbackground\n", encoding="utf-8")
            positions = [start + offset for start in range(0, 1000, 20) for offset in (4, 14, 17)]
            rows = [
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4",
            ]
            rng = random.Random(201)
            for pos in positions:
                mask = rng.randint(1, 7)
                alleles = [bin(haplotype & mask).count("1") % 2 for haplotype in range(8)]
                genotypes = "\t".join(f"{alleles[i]}|{alleles[i + 1]}" for i in range(0, 8, 2))
                rows.append(f"1\t{pos + 1}\tv{pos}\tA\tG\t.\tPASS\t.\tGT\t{genotypes}")
            (directory / "cohort.vcf").write_text("\n".join(rows) + "\n", encoding="utf-8")
            args = existing_tests.BackgroundBlockLDTests.args(directory, "range")
            args.snps_per_block = None
            args.min_snps_per_block, args.max_snps_per_block = 1, 2
            args.background_pairs_per_distance = 8
            args.background_regions = 5
            args.max_attempt_factor = 100
            result = run(args)
            self.assertGreaterEqual(result["background_pair_count"], 32)
            read_rows = existing_tests.BackgroundBlockLDTests.read_rows
            targets = read_rows(directory / "range.target_pairs.tsv")
            backgrounds = read_rows(directory / "range.background_pairs.tsv")
            expected_counts = {("left", 1): 2, ("left", 2): 1, ("right", 1): 1, ("right", 2): 2}
            self.assertEqual(region_snp_counts(targets, [(1, 1), (1, 2), (2, 1), (2, 2)]), expected_counts)
            target_keys = {pair_background_key(row, True) for row in targets}
            self.assertEqual(target_keys, {(110, 2, 1), (120, 2, 2), (120, 1, 1), (130, 1, 2)})
            for target in targets:
                key = pair_background_key(target, True)
                references = [float(row["multiallelic_dprime"]) for row in backgrounds if pair_background_key(row, True) == key]
                self.assertEqual(len(references), 8)
                expected_p = empirical_tail_p(float(target["multiallelic_dprime"]), sorted(references))
                self.assertAlmostEqual(float(target["pair_nominal_p"]), expected_p)
                self.assertEqual(target["pair_snp_count_stratum"], f"{key[1]},{key[2]}")
            for background in backgrounds:
                key = pair_background_key(background, True)
                self.assertEqual(int(background["distance_stratum_bp"]), key[0])
                self.assertEqual(int(background["snp_number_stratum_block1"]), key[1])
                self.assertEqual(int(background["snp_number_stratum_block2"]), key[2])
            regions = read_rows(directory / "range.background_regions.tsv")
            self.assertEqual(len(regions), 5)
            for region in regions:
                start, end = int(region["region_start"]), int(region["region_end"])
                self.assertGreaterEqual(start - 20, 0)
                self.assertLessEqual(end + 20, 1000)
                for side, rank in expected_counts:
                    window_start = start - rank * 10 if side == "left" else end + (rank - 1) * 10
                    selected_count = min(2, sum(window_start <= pos < window_start + 10 for pos in positions))
                    self.assertIn(selected_count, (1, 2), (region["background_region_id"], side, rank))
                self.assertEqual(region["required_block_snp_counts"], "NA")
                self.assertEqual(region["target_pair_mask"], "NA")
                self.assertEqual(region["n_eligible_left_blocks"], "2")
                self.assertEqual(region["n_eligible_right_blocks"], "2")
            metadata = json.loads((directory / "range.run.json").read_text(encoding="utf-8"))
            self.assertEqual(metadata["resolved_snp_limits"], {"min_snps_per_block": 1, "max_snps_per_block": 2})
            self.assertEqual(metadata["region_snp_count_matching"], "shared_range_and_entropy_cap")
            self.assertTrue(any(
                info["selected_snp_count_rejections"] > 0
                for info in metadata["pair_sampling_counts"].values()
            ))


if __name__ == "__main__":
    unittest.main()

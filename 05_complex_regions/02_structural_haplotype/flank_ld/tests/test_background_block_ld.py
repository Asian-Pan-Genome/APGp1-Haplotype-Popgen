import argparse
import csv
import json
import random
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from flanking_block_ld import background_block_ld as workflow
from flanking_block_ld.background_block_ld import (
    StartDomain,
    analyze_region,
    empirical_tail_p,
    merge_background_intervals,
    pair_sampling_domains,
    region_sampling_domains,
    random_starts,
    run,
    subtract_exclusions,
)
from flanking_block_ld.flanking_block_ld import Region


class BackgroundBlockLDTests(unittest.TestCase):
    @staticmethod
    def region(start, end, chrom="1", name="background"):
        return Region(chrom, f"chr{chrom}", start, end, name, 1)

    @staticmethod
    def bounds(regions):
        return [(region.chrom, region.start, region.end) for region in regions]

    @staticmethod
    def domain_starts(domains):
        return {
            (domain.chrom, start)
            for domain in domains
            for start in range(domain.start, domain.end + 1)
        }

    def test_background_intervals_merge_overlaps_and_touching_edges(self):
        merged = merge_background_intervals([
            self.region(20, 30),
            self.region(5, 15),
            self.region(0, 10),
            self.region(15, 20),
            self.region(0, 10, chrom="2"),
        ])
        self.assertEqual(self.bounds(merged), [("1", 0, 30), ("2", 0, 10)])

    def test_exclusion_subtraction_clips_and_preserves_other_chromosomes(self):
        retained = subtract_exclusions(
            [self.region(10, 100), self.region(0, 20, chrom="2")],
            [
                self.region(0, 20),
                self.region(30, 40),
                self.region(35, 50),
                self.region(90, 110),
                self.region(0, 100, chrom="3"),
            ],
        )
        self.assertEqual(
            self.bounds(retained),
            [("1", 20, 30), ("1", 50, 90), ("2", 0, 20)],
        )

    def test_pair_sampling_allows_intervening_gap_but_contains_both_windows(self):
        background = [self.region(0, 20), self.region(40, 60)]
        domains = pair_sampling_domains(background, window_bp=10, distance_bp=40)
        self.assertEqual(self.domain_starts(domains), {("1", i) for i in range(11)})
        self.assertTrue(all(domain.display_chrom == "chr1" for domain in domains))
        self.assertEqual(
            pair_sampling_domains(background, window_bp=10, distance_bp=20), []
        )

    def test_pair_sampling_does_not_bridge_chromosomes(self):
        domains = pair_sampling_domains(
            [self.region(0, 20), self.region(40, 60, chrom="2")],
            window_bp=10,
            distance_bp=40,
        )
        self.assertEqual(domains, [])

    def test_pair_start_domains_match_exhaustive_integer_containment(self):
        background = [self.region(0, 23), self.region(35, 61), self.region(80, 94)]
        width, distance = 7, 35
        observed = self.domain_starts(
            pair_sampling_domains(background, width, distance)
        )
        expected = set()
        for start in range(95):
            left_contained = any(r.start <= start and start + width <= r.end for r in background)
            right_contained = any(
                r.start <= start + distance and start + distance + width <= r.end
                for r in background
            )
            if left_contained and right_contained:
                expected.add(("1", start))
        self.assertEqual(observed, expected)

    def test_region_domains_contain_entire_region_and_both_flanks(self):
        domains = region_sampling_domains(
            [self.region(0, 40), self.region(60, 80)], central_bp=10, flank_span=10
        )
        self.assertEqual(self.domain_starts(domains), {("1", i) for i in range(10, 21)})
        self.assertEqual(
            region_sampling_domains(
                [self.region(0, 20), self.region(40, 60)],
                central_bp=10,
                flank_span=10,
            ),
            [],
        )
        exact = region_sampling_domains([self.region(100, 130)], 10, 10)
        self.assertEqual(self.domain_starts(exact), {("1", 110)})

    def test_empirical_tail_counts_ties_and_applies_plus_one(self):
        background = [0.1, 0.2, 0.2, 0.9]
        self.assertAlmostEqual(empirical_tail_p(0.2, background), 4 / 5)
        self.assertAlmostEqual(empirical_tail_p(1.0, background), 1 / 5)
        self.assertAlmostEqual(empirical_tail_p(0.0, background), 1.0)

    def test_empirical_tail_rejects_empty_or_nonfinite_inputs(self):
        with self.assertRaises(ValueError):
            empirical_tail_p(0.5, [])
        with self.assertRaises(ValueError):
            empirical_tail_p(float("nan"), [0.1])

    def test_empirical_tail_treats_roundoff_as_ties_not_real_differences(self):
        self.assertEqual(empirical_tail_p(1.0, [1.0 - 2e-16, 1.0]), 1.0)
        self.assertEqual(empirical_tail_p(1.0, [0.99999]), 0.5)

    def test_random_start_sampling_visits_all_legal_positions_without_replacement(self):
        domains = [StartDomain("1", "chr1", 5, 5), StartDomain("2", "chr2", 10, 14)]
        first = [(domain.chrom, start) for domain, start in random_starts(domains, 100, random.Random(17))]
        second = [(domain.chrom, start) for domain, start in random_starts(domains, 100, random.Random(17))]
        self.assertEqual(first, second)
        self.assertEqual(len(first), 6)
        self.assertEqual(set(first), {("1", 5)} | {("2", i) for i in range(10, 15)})

    @staticmethod
    def pair_row(left_rank, right_rank, dprime, nominal_p, status="ok"):
        left_start = (left_rank - 1) * 5
        right_start = 30 + (right_rank - 1) * 5
        return {
            "left_block_rank": str(left_rank),
            "right_block_rank": str(right_rank),
            "block1_rank": str(left_rank),
            "block2_rank": str(right_rank),
            "block1_side": "left",
            "block2_side": "right",
            "block1_label": f"L{left_rank}",
            "block2_label": f"R{right_rank}",
            "left_window_start": str(left_start),
            "left_window_end": str(left_start + 10),
            "right_window_start": str(right_start),
            "right_window_end": str(right_start + 10),
            "block_center_distance_bp": str(right_start - left_start),
            "multiallelic_dprime": str(dprime),
            "pair_nominal_p": str(nominal_p),
            "pair_type": "left_right",
            "status": status,
        }

    def test_region_summary_uses_target_mask_and_union_window_coverage(self):
        rows = [
            self.pair_row(1, 1, 0.9, 0.001),
            self.pair_row(1, 2, 0.8, 0.1),
            self.pair_row(2, 1, 0.7, 0.5),
            self.pair_row(2, 2, 1.0, 0.001, status="insufficient_left_snps"),
        ]
        summary = analyze_region(rows, [(1, 1), (1, 2), (2, 1)], 0.01)
        self.assertAlmostEqual(summary["max_dprime"], 0.9)
        self.assertEqual(summary["n_valid_pairs"], 3)
        self.assertEqual(summary["n_significant_pairs"], 1)
        self.assertAlmostEqual(summary["significant_pair_fraction"], 1 / 3)
        self.assertEqual(summary["n_supported_left_blocks"], 1)
        self.assertEqual(summary["n_supported_right_blocks"], 1)
        self.assertAlmostEqual(summary["supported_left_block_fraction"], 0.5)
        self.assertAlmostEqual(summary["supported_right_block_fraction"], 0.5)
        self.assertAlmostEqual(summary["supported_left_window_coverage_fraction"], 10 / 15)
        self.assertAlmostEqual(summary["supported_right_window_coverage_fraction"], 10 / 15)
        self.assertEqual(summary["peak_block1_label"], "L1")
        self.assertEqual(summary["peak_block2_label"], "R1")
        self.assertAlmostEqual(summary["peak_pair_nominal_p"], 0.001)

    def test_region_summary_rejects_background_missing_a_target_valid_edge(self):
        with self.assertRaises(ValueError):
            analyze_region(
                [self.pair_row(1, 1, 0.9, 0.001, "insufficient_left_snps")],
                [(1, 1)],
                0.01,
            )

    def test_pair_threshold_is_strict_less_than(self):
        summary = analyze_region([self.pair_row(1, 1, 0.9, 0.01)], [(1, 1)], 0.01)
        self.assertEqual(summary["n_significant_pairs"], 0)

    @staticmethod
    def read_rows(path):
        with path.open(encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    @staticmethod
    def args(directory, prefix):
        return argparse.Namespace(
            regions_bed=directory / "target.bed",
            background_bed=directory / "background.bed",
            vcf=directory / "cohort.vcf",
            vcf_access="preload",  # Plain-VCF legacy fixtures explicitly exercise full-scan compatibility.
            output_prefix=directory / prefix,
            window_bp=10,
            window_step_bp=None,
            block_count=2,
            max_flank_bp=20,
            snps_per_block=1,
            background_pairs_per_distance=8,
            pair_distance_tolerance=0.0,
            pair_match_snp_counts=True,  # Preserve legacy exact-pair regression fixtures.
            background_regions=5,
            region_length_tolerance_bp=0,  # These legacy fixtures test exact-length geometry.
            region_block_count_tolerance=0,
            region_block_count_allow_swap=False,  # Retain strict-count regression coverage.
            region_entropy_max_wasserstein=0.05,
            skip_region_entropy_matching=True,  # Legacy fixtures isolate the earlier controls.
            pair_significance_threshold=0.01,
            max_attempt_factor=20,
            seed=13,
            min_maf=0.1,
            max_site_missing_rate=0.05,
            max_window_missing_filtered_fraction=0.10,
            max_exact_entropy_combinations=100_000,
            min_complete_haplotypes=2,
            include_filtered=False,
            allow_unphased_missing=False,
            include_target_in_background=False,
        )

    def test_dense_phased_vcf_workflow_is_deterministic_and_contains_regions(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            (directory / "target.bed").write_text("chr1\t300\t400\ttarget\n", encoding="utf-8")
            (directory / "background.bed").write_text("chr1\t0\t1000\tbackground\n", encoding="utf-8")
            rows = [
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4",
            ]
            for pos in range(5, 1000, 5):
                rows.append(
                    f"1\t{pos}\tv{pos}\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|1\t0|1\t0|1"
                )
            (directory / "cohort.vcf").write_text("\n".join(rows) + "\n", encoding="utf-8")
            run(self.args(directory, "first"))
            run(self.args(directory, "second"))
            suffixes = (
                ".target_pairs.tsv", ".summary.tsv", ".background_pairs.tsv",
                ".background_regions.tsv", ".background_regions.bed", ".run.json",
            )
            for suffix in suffixes:
                self.assertTrue((directory / f"first{suffix}").is_file(), suffix)
                if suffix != ".run.json":
                    self.assertEqual(
                        (directory / f"first{suffix}").read_text(encoding="utf-8"),
                        (directory / f"second{suffix}").read_text(encoding="utf-8"),
                        suffix,
                    )
            target_pairs = self.read_rows(directory / "first.target_pairs.tsv")
            self.assertEqual(len(target_pairs), 4)
            self.assertTrue(all(row["status"] == "ok" for row in target_pairs))
            self.assertTrue(all(float(row["multiallelic_dprime"]) == 1.0 for row in target_pairs))
            self.assertTrue(all(float(row["pair_nominal_p"]) == 1.0 for row in target_pairs))
            with (directory / "first.target_pairs.tsv").open(encoding="utf-8") as handle:
                header = handle.readline().rstrip("\n").split("\t")
            self.assertEqual(header[7:11], [
                "multiallelic_dprime", "block1_snp_number_pass_threshold",
                "block2_snp_number_pass_threshold", "both_blocks_pass_snp_number_threshold",
            ])
            self.assertEqual(header[11:14], [
                "block1_missed_snp_number", "block2_missed_snp_number", "total_missed_snp_number",
            ])
            self.assertEqual(header[14], "pair_nominal_p")
            self.assertTrue(all(row["block1_snp_number_pass_threshold"] == "1" for row in target_pairs))
            metadata = json.loads((directory / "first.run.json").read_text(encoding="utf-8"))
            self.assertEqual(metadata["status"], "complete")
            self.assertEqual(metadata["region_qc_matching"],
                             "available_block_counts_with_tolerance_and_optional_swap; each region uses all its own QC-valid cross-flank pairs")
            summary = self.read_rows(directory / "first.summary.tsv")[0]
            self.assertEqual(summary["status"], "ok")
            self.assertEqual(int(summary["n_valid_pairs"]), 4)
            self.assertEqual(int(summary["n_significant_pairs"]), 0)
            self.assertEqual(float(summary["significant_pair_fraction"]), 0.0)
            self.assertEqual(float(summary["region_max_empirical_p"]), 1.0)
            self.assertEqual(float(summary["region_breadth_empirical_p"]), 1.0)
            background_pairs = self.read_rows(directory / "first.background_pairs.tsv")
            self.assertEqual(len(background_pairs), 24)
            for distance in (110, 120, 130):
                stratum = [row for row in background_pairs if int(row["distance_stratum_bp"]) == distance]
                self.assertEqual(len(stratum), 8)
                self.assertTrue(all(float(row["block_center_distance_bp"]) == distance for row in stratum))
            background_regions = self.read_rows(directory / "first.background_regions.tsv")
            self.assertEqual(len(background_regions), 5)
            for row in background_regions:
                start, end = int(row["region_start"]), int(row["region_end"])
                self.assertEqual(end - start, 100)
                self.assertGreaterEqual(start - 20, 0)
                self.assertLessEqual(end + 20, 1000)
                self.assertTrue(end + 20 <= 280 or start - 20 >= 420)

            pair_path = directory / "first.background_pairs.tsv"
            region_path = directory / "first.background_regions.tsv"
            original_pair_background = pair_path.read_text()
            original_region_background = region_path.read_text()
            with mock.patch.object(
                    workflow, "random_pair_windows",
                    side_effect=AssertionError("pair backgrounds must be reused")), \
                    mock.patch.object(
                        workflow, "random_region_intervals",
                        side_effect=AssertionError("region backgrounds must be reused")):
                reused = run(self.args(directory, "first"))
            self.assertEqual(reused["background_pair_count"], 24)
            self.assertEqual(reused["background_region_count"], 5)
            self.assertEqual(pair_path.read_text(), original_pair_background)
            self.assertEqual(region_path.read_text(), original_region_background)
            log = (directory / "first.log").read_text()
            self.assertIn("PAIR BACKGROUND REUSE", log)
            self.assertIn("REGION BACKGROUND REUSE", log)

            incompatible = self.args(directory, "first")
            incompatible.pair_distance_tolerance = 0.1
            with self.assertRaisesRegex(ValueError, "incompatible with current controls"):
                run(incompatible)

    def test_heterogeneous_workflow_p_values_use_the_correct_background_statistic(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
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
                genotypes = "\t".join(f"{alleles[index]}|{alleles[index + 1]}" for index in range(0, 8, 2))
                rows.append(f"1\t{pos}\tv{pos}\tA\tG\t.\tPASS\t.\tGT\t{genotypes}")
            (directory / "cohort.vcf").write_text("\n".join(rows) + "\n", encoding="utf-8")
            args = self.args(directory, "heterogeneous")
            args.background_pairs_per_distance = 40
            args.background_regions = 20
            args.pair_significance_threshold = 0.3
            run(args)
            target_pairs = self.read_rows(directory / "heterogeneous.target_pairs.tsv")
            background_pairs = self.read_rows(directory / "heterogeneous.background_pairs.tsv")
            for target_pair in target_pairs:
                distance = float(target_pair["block_center_distance_bp"])
                background_values = [
                    float(pair["multiallelic_dprime"])
                    for pair in background_pairs
                    if float(pair["distance_stratum_bp"]) == distance
                ]
                self.assertEqual(len(background_values), 40)
                expected_p = (1 + sum(value >= float(target_pair["multiallelic_dprime"]) - 1e-12 for value in background_values)) / 41
                self.assertAlmostEqual(float(target_pair["pair_nominal_p"]), expected_p)
            observed = self.read_rows(directory / "heterogeneous.summary.tsv")[0]
            backgrounds = self.read_rows(directory / "heterogeneous.background_regions.tsv")
            for statistic, p_field in (
                ("max_dprime", "region_max_empirical_p"),
                ("significant_pair_fraction", "region_breadth_empirical_p"),
            ):
                expected_p = (1 + sum(float(row[statistic]) >= float(observed[statistic]) - 1e-12 for row in backgrounds)) / 21
                self.assertAlmostEqual(float(observed[p_field]), expected_p)
            self.assertNotEqual(float(observed["region_max_empirical_p"]), float(observed["region_breadth_empirical_p"]))

    def test_no_valid_target_preserves_qc_and_does_not_mark_run_complete(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            (directory / "target.bed").write_text("1\t300\t400\ttarget\n", encoding="utf-8")
            (directory / "background.bed").write_text("1\t0\t1000\tbackground\n", encoding="utf-8")
            (directory / "cohort.vcf").write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n"
                "1\t295\tleft\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|1\t0|1\t0|1\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "no target cross-flank pair passes QC"):
                run(self.args(directory, "failed"))
            self.assertFalse((directory / "failed.run.json").exists())
            self.assertFalse((directory / "failed.summary.tsv").exists())
            qc_rows = self.read_rows(directory / "failed.target_pairs.tsv")
            self.assertEqual(len(qc_rows), 4)
            self.assertTrue(all("insufficient_right_snps" in row["status"] for row in qc_rows))
            self.assertTrue(all(row["pair_nominal_p"] == "NA" for row in qc_rows))


if __name__ == "__main__":
    unittest.main()

import contextlib
import io
import json
import math
import random
import tempfile
import threading
import unittest
from pathlib import Path

from flanking_block_ld import background_block_ld as workflow
from flanking_block_ld.flanking_block_ld import Region, read_vcf
from flanking_block_ld.tests import test_background_block_ld as existing_tests


class BackgroundPairSoftControlTests(unittest.TestCase):
    def test_cli_defaults_and_legacy_opt_in(self):
        argv = ["target.bed", "background.bed", "cohort.vcf.gz", "-o", "result"]
        args = workflow.build_parser().parse_args(argv)
        self.assertEqual(args.pair_distance_tolerance, 0.10)
        self.assertFalse(args.pair_match_snp_counts)
        self.assertEqual(args.region_entropy_max_wasserstein, 0.05)
        self.assertEqual(args.threads, 4)
        args = workflow.build_parser().parse_args(
            argv + ["--pair-distance-tolerance", "0", "--pair-match-snp-counts"])
        self.assertEqual(args.pair_distance_tolerance, 0)
        self.assertTrue(args.pair_match_snp_counts)
        # Rendering help exercises argparse's percent-sign interpolation.
        self.assertIn("+/-10%", workflow.build_parser().format_help())

    def test_ordered_thread_map_runs_concurrently_but_preserves_order(self):
        barrier = threading.Barrier(4)
        lock = threading.Lock()
        active = peak = 0

        def work(value):
            nonlocal active, peak
            with lock:
                active += 1
                peak = max(peak, active)
            if value < 4:
                barrier.wait(timeout=2)
            with lock:
                active -= 1
            return value * value

        self.assertEqual(list(workflow.ordered_thread_map(work, range(12), 4)),
                         [value * value for value in range(12)])
        self.assertEqual(peak, 4)

    def test_invalid_thread_count_fails_cli(self):
        argv = ["target.bed", "background.bed", "cohort.vcf.gz", "-o", "result"]
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as error:
            workflow.main(argv + ["--threads", "0"])
        self.assertEqual(error.exception.code, 2)

    def test_distance_band_rounding_and_nonoverlap(self):
        self.assertEqual(workflow.pair_distance_bounds(100_000, 0.10, 10_000), (90_000, 110_000))
        self.assertEqual(workflow.pair_distance_bounds(110, 0.10, 10), (99, 121))
        self.assertEqual(workflow.pair_distance_bounds(103, 0.10, 10), (93, 113))
        self.assertEqual(workflow.pair_distance_bounds(10, 0.10, 10), (10, 11))
        self.assertEqual(workflow.pair_distance_bounds(103, 0, 10), (103, 103))
        for tolerance in (-0.1, 1, math.nan, math.inf):
            with self.subTest(tolerance=tolerance), self.assertRaises(ValueError):
                workflow.pair_distance_bounds(100, tolerance, 10)
        with self.assertRaises(ValueError):
            workflow.pair_distance_bounds(9, 0.10, 10)

    def test_invalid_tolerance_fails_cli_before_loading_inputs(self):
        for value in ("-0.1", "1", "nan", "inf"):
            with self.subTest(value=value), contextlib.redirect_stderr(io.StringIO()), \
                    self.assertRaises(SystemExit) as error:
                workflow.main(["t.bed", "b.bed", "c.vcf", "-o", "result",
                               "--pair-distance-tolerance", value])
            self.assertEqual(error.exception.code, 2)

    def test_domain_count_and_decode_match_exhaustive_integer_pairs(self):
        rng = random.Random(107)
        for _ in range(200):
            first_start = rng.randrange(20)
            first_end = first_start + rng.randrange(10)
            second_start = rng.randrange(30)
            second_end = second_start + rng.randrange(10)
            lower, upper = sorted((rng.randrange(1, 20), rng.randrange(1, 20)))
            expected = [(a, b) for a in range(first_start, first_end + 1)
                        for b in range(second_start, second_end + 1) if lower <= b - a <= upper]
            domain = workflow.PairDistanceDomain("1", "chr1", first_start, first_end,
                                                  second_start, second_end, lower, upper)
            self.assertEqual(domain.size, len(expected))
            self.assertEqual([domain.decode(i) for i in range(domain.size)], expected)
            with self.assertRaises(ValueError):
                domain.decode(domain.size)

    def test_band_domains_containment_chromosomes_gaps_and_no_duplicate_pairs(self):
        region = existing_tests.BackgroundBlockLDTests.region
        background = [region(0, 13), region(8, 23), region(35, 61), region(80, 94),
                      region(0, 23, chrom="2"), region(35, 61, chrom="2")]
        width, distance = 7, 35
        domains = workflow.pair_distance_sampling_domains(background, width, distance, 0.10)
        observed = [(domain.chrom, *domain.decode(i)) for domain in domains for i in range(domain.size)]
        merged = workflow.merge_background_intervals(background)
        starts = [(r.chrom, a) for r in merged for a in range(r.start, r.end - width + 1)]
        lower, upper = workflow.pair_distance_bounds(distance, 0.10, width)
        expected = {(chrom, a, b) for chrom, a in starts for other_chrom, b in starts
                    if chrom == other_chrom and lower <= b - a <= upper}
        self.assertEqual(set(observed), expected)
        self.assertEqual(len(observed), len(expected))
        # The intervening gap can cross a BED gap, but windows cannot.
        self.assertTrue(any(a < 23 and b >= 35 for _, a, b in observed))
        sampled = [(domain.chrom, a, b) for domain, a, b in
                   workflow.random_pair_windows(domains, 10_000, random.Random(13))]
        self.assertEqual(set(sampled), expected)
        self.assertEqual(len(sampled), len(expected))

    def test_zero_tolerance_matches_legacy_exact_geometry_and_order(self):
        background = [existing_tests.BackgroundBlockLDTests.region(0, 23),
                      existing_tests.BackgroundBlockLDTests.region(35, 61)]
        old = [(domain.chrom, start, start + 35) for domain, start in
               workflow.random_starts(workflow.pair_sampling_domains(background, 7, 35),
                                      1000, random.Random(13))]
        new = [(domain.chrom, a, b) for domain, a, b in workflow.random_pair_windows(
            workflow.pair_distance_sampling_domains(background, 7, 35, 0), 1000, random.Random(13))]
        self.assertEqual(new, old)

    def test_large_distance_band_does_not_enumerate_bp_configurations(self):
        domains = workflow.pair_distance_sampling_domains(
            [existing_tests.BackgroundBlockLDTests.region(0, 100_000_000)], 10_000, 1_000_000, 0.10)
        self.assertEqual(len(domains), 1)
        expected = sum(100_000_000 - 10_000 - distance + 1
                       for distance in range(900_000, 1_100_001))
        self.assertEqual(domains[0].size, expected)
        first = [(a, b) for _, a, b in workflow.random_pair_windows(domains, 100, random.Random(1))]
        second = [(a, b) for _, a, b in workflow.random_pair_windows(domains, 100, random.Random(1))]
        self.assertEqual(first, second)
        self.assertEqual(len(set(first)), 100)
        self.assertTrue(all(900_000 <= b - a <= 1_100_000 for a, b in first))

    def test_pair_p_ignores_snp_counts_and_uses_reference_bank_not_sampled_distance(self):
        rows = [existing_tests.BackgroundBlockLDTests.pair_row(1, 1, 0.8, 0.005) for _ in range(2)]
        for row, counts in zip(rows, ((1, 2), (2, 1))):
            row.update(block1_snp_number_pass_threshold=str(counts[0]),
                       block2_snp_number_pass_threshold=str(counts[1]))
        distance = int(float(rows[0]["block_center_distance_bp"]))
        self.assertEqual(workflow.pair_background_key(rows[0]), (distance, 0, 0))
        workflow.attach_pair_p(rows, {(distance, 0, 0): [0.1, 0.2]}, 0.5, window_bp=10)
        for row in rows:
            self.assertAlmostEqual(float(row["pair_nominal_p"]), 1 / 3)
            self.assertEqual(row["pair_snp_count_stratum"], "NA")
            self.assertEqual(row["pair_snp_count_matching"], "none")
            self.assertEqual(row["pair_distance_background_n"], "2")
        background = dict(rows[0], distance_stratum_bp=str(distance),
                          block_center_distance_bp=str(distance + 1))
        self.assertEqual(workflow.pair_background_key(background), (distance, 0, 0))

    @staticmethod
    def write_inputs(directory):
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
                genotypes = "\t".join(f"{alleles[i]}|{alleles[i + 1]}" for i in range(0, 8, 2))
                rows.append(f"1\t{pos + 1}\tv{pos}\tA\tG\t.\tPASS\t.\tGT\t{genotypes}")
        (directory / "cohort.vcf").write_text("\n".join(rows) + "\n", encoding="utf-8")

    @staticmethod
    def args(directory, prefix, skip=False):
        argv = [str(directory / "target.bed"), str(directory / "background.bed"),
                str(directory / "cohort.vcf"), "-o", str(directory / prefix),
                "--vcf-access", "preload", "--window-bp", "10", "--block-count", "2",
                "--max-flank-bp", "20", "--min-snps-per-block", "1", "--max-snps-per-block", "2",
                "--background-pairs-per-distance", "12", "--background-regions", "6",
                "--region-length-tolerance-bp", "3", "--pair-significance-threshold", "0.3",
                "--max-attempt-factor", "100", "--seed", "13"]
        if skip:
            argv.append("--skip-background-region-pair-p")
        return workflow.build_parser().parse_args(argv)

    def test_default_workflow_soft_pairs_region_entropy_and_skip_mode(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            self.write_inputs(directory)
            for prefix, skip in (("full", False), ("skip", True)):
                args = self.args(directory, prefix, skip)
                with contextlib.redirect_stderr(io.StringIO()):
                    workflow.run(args)
            read_rows = existing_tests.BackgroundBlockLDTests.read_rows
            targets = read_rows(directory / "full.target_pairs.tsv")
            backgrounds = read_rows(directory / "full.background_pairs.tsv")
            banks = {}
            for row in backgrounds:
                key = workflow.pair_background_key(row)
                banks.setdefault(key, []).append(float(row["multiallelic_dprime"]))
                actual = int(float(row["block_center_distance_bp"]))
                self.assertLessEqual(int(row["pair_distance_min_bp"]), actual)
                self.assertLessEqual(actual, int(row["pair_distance_max_bp"]))
                self.assertEqual(row["snp_number_stratum_block1"], "NA")
                self.assertEqual(row["snp_number_stratum_block2"], "NA")
                self.assertEqual(row["pair_snp_count_matching"], "none")
            self.assertTrue(any(int(row["distance_stratum_bp"]) != int(float(row["block_center_distance_bp"]))
                                for row in backgrounds))
            self.assertEqual(len({workflow.pair_background_key(row) for row in targets}), 3)
            for row in targets:
                reference = banks[workflow.pair_background_key(row)]
                self.assertEqual(len(reference), 12)
                self.assertAlmostEqual(float(row["pair_nominal_p"]), workflow.empirical_tail_p(
                    float(row["multiallelic_dprime"]), sorted(reference)))
                self.assertEqual(row["pair_distance_background_n"], "12")
            target = next(row for row in targets if workflow.pair_background_key(row)[0] == 110)
            target_counts = (target["block1_snp_number_pass_threshold"], target["block2_snp_number_pass_threshold"])
            self.assertTrue(any((row["block1_snp_number_pass_threshold"], row["block2_snp_number_pass_threshold"])
                                != target_counts for row in backgrounds if int(row["distance_stratum_bp"]) == 110))
            data = read_vcf(directory / "cohort.vcf", [Region("1", "1", 0, 1000, "all", 1)],
                            0, 0.1, 0.95, False, False, include_region_interiors=True)
            store = workflow.BlockStore(data, self.args(directory, "check"))
            regions = read_rows(directory / "full.background_regions.tsv")
            for bg in regions:
                self.assertLessEqual(float(bg["region_entropy_wasserstein"]), 0.05 + 1e-12)
                rows = [row for row in store.region_pairs(Region(
                    "1", "1", int(bg["region_start"]), int(bg["region_end"]), "check", 0))
                        if workflow.valid_row(row)]
                ps = [workflow.empirical_tail_p(float(row["multiallelic_dprime"]),
                                               sorted(banks[workflow.pair_background_key(row)])) for row in rows]
                self.assertAlmostEqual(float(bg["significant_pair_fraction"]),
                                       sum(p < 0.3 for p in ps) / len(ps))
            self.assertEqual((directory / "full.target_pairs.tsv").read_text(),
                             (directory / "skip.target_pairs.tsv").read_text())
            self.assertEqual((directory / "full.background_regions.bed").read_text(),
                             (directory / "skip.background_regions.bed").read_text())
            full_summary = read_rows(directory / "full.summary.tsv")[0]
            skip_summary = read_rows(directory / "skip.summary.tsv")[0]
            self.assertEqual(full_summary["region_max_empirical_p"], skip_summary["region_max_empirical_p"])
            self.assertEqual(skip_summary["region_breadth_empirical_p"], "NA")
            metadata = json.loads((directory / "full.run.json").read_text())
            self.assertEqual(metadata["pair_snp_count_matching"], "none")
            self.assertEqual(metadata["pair_distance_tolerance"], 0.10)
            self.assertTrue(any(bank["source"] == "background_region_pairs"
                                for bank in metadata["pair_sampling_counts"].values()))
            self.assertTrue(all(bank["selected_snp_count_rejections"] == 0
                                for bank in metadata["pair_sampling_counts"].values()))
            self.assertTrue(all(":" not in key for key in metadata["pair_sampling_counts"]))
            self.assertIn("snps=unmatched", (directory / "full.log").read_text())


    def test_distance_and_snp_count_controls_can_be_selected_independently(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            self.write_inputs(directory)
            for prefix, tolerance, match_counts in (("exact_mixed", 0, False), ("soft_matched", 0.1, True)):
                args = self.args(directory, prefix, skip=True)
                args.pair_distance_tolerance = tolerance
                args.pair_match_snp_counts = match_counts
                with contextlib.redirect_stderr(io.StringIO()):
                    workflow.run(args)
                backgrounds = existing_tests.BackgroundBlockLDTests.read_rows(directory / f"{prefix}.background_pairs.tsv")
                targets = existing_tests.BackgroundBlockLDTests.read_rows(directory / f"{prefix}.target_pairs.tsv")
                banks = {}
                for row in backgrounds:
                    key = workflow.pair_background_key(row, match_counts)
                    banks.setdefault(key, []).append(float(row["multiallelic_dprime"]))
                    if match_counts:
                        self.assertEqual(row["snp_number_stratum_block1"], row["block1_snp_number_pass_threshold"])
                        self.assertEqual(row["snp_number_stratum_block2"], row["block2_snp_number_pass_threshold"])
                    if tolerance == 0:
                        self.assertEqual(int(row["distance_stratum_bp"]), int(float(row["block_center_distance_bp"])))
                self.assertEqual(len(banks), 4 if match_counts else 3)
                for row in targets:
                    reference = banks[workflow.pair_background_key(row, match_counts)]
                    self.assertEqual(len(reference), 12)
                    self.assertAlmostEqual(float(row["pair_nominal_p"]), workflow.empirical_tail_p(
                        float(row["multiallelic_dprime"]), sorted(reference)))

    def test_indexed_soft_pair_background_matches_preload_and_prefiltered_mode(self):
        from flanking_block_ld.tests.test_indexed_vcf import FakeBCFtools

        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            self.write_inputs(directory)
            Path(str(directory / "cohort.vcf") + ".csi").touch()
            for prefiltered in (False, True):
                fake = FakeBCFtools((directory / "cohort.vcf").read_text())
                for access in ("preload", "bcftools"):
                    args = self.args(directory, f"{access}_{prefiltered}")
                    args.vcf_access = access
                    args.vcf_prefiltered = prefiltered
                    with fake.patches(), contextlib.redirect_stderr(io.StringIO()):
                        workflow.run(args)
                for suffix in ("summary.tsv", "target_pairs.tsv", "background_pairs.tsv",
                               "background_regions.tsv", "background_regions.bed"):
                    self.assertEqual((directory / f"preload_{prefiltered}.{suffix}").read_text(),
                                     (directory / f"bcftools_{prefiltered}.{suffix}").read_text(), suffix)
                metadata = json.loads((directory / f"bcftools_{prefiltered}.run.json").read_text())
                self.assertEqual(metadata["vcf_access"]["mode"], "bcftools")
                self.assertEqual(metadata["pair_snp_count_matching"], "none")
                self.assertEqual(metadata["pair_distance_tolerance"], 0.10)
                self.assertGreater(metadata["vcf_access"]["query_batches"], 0)
                self.assertTrue(all(len(spans) <= 2 for spans in fake.queries))

    def test_one_and_four_threads_produce_identical_result_tables(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            self.write_inputs(directory)
            for threads in (1, 4):
                args = self.args(directory, f"threads_{threads}")
                args.threads = threads
                with contextlib.redirect_stderr(io.StringIO()):
                    workflow.run(args)
            for suffix in ("summary.tsv", "target_pairs.tsv", "background_pairs.tsv",
                           "background_regions.tsv", "background_regions.bed"):
                self.assertEqual((directory / f"threads_1.{suffix}").read_text(),
                                 (directory / f"threads_4.{suffix}").read_text(), suffix)
            one = json.loads((directory / "threads_1.run.json").read_text())
            four = json.loads((directory / "threads_4.run.json").read_text())
            self.assertEqual(one["parallelism"]["workers"], 1)
            self.assertEqual(four["parallelism"]["workers"], 4)
            self.assertIn("pair_reference_banks", four["parallelism"]["parallel_units"])
            self.assertIn("background_region_candidates", four["parallelism"]["parallel_units"])
            self.assertIn("BACKGROUND THREADS workers=4", (directory / "threads_4.log").read_text())


if __name__ == "__main__":
    unittest.main()

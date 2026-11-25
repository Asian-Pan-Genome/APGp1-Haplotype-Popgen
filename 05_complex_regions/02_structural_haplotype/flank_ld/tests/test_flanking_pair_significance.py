import csv
import gzip
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from flanking_block_ld import summarize_flanking_pair_significance as summary


class FlankingPairSignificanceTests(unittest.TestCase):
    @staticmethod
    def region(matrix, name="SCR1"):
        return summary.RegionPairs(
            dict(region_id=name, chrom="chr1", region_start="10000", region_end="20000"),
            {(i, j): (p is not None, p is not None and p < 0.01)
             for i, values in enumerate(matrix, 1) for j, p in enumerate(values, 1)},
        )

    @staticmethod
    def rows(matrix, name="SCR1"):
        return [dict(
            region_id=name, chrom="chr1", region_start="10000", region_end="20000",
            block1_side="left", block2_side="right", block1_rank=str(i), block2_rank=str(j),
            status="ok" if p is not None else "insufficient_left_snps",
            multiallelic_dprime="0.5" if p is not None else "NA",
            pair_nominal_p=str(p) if p is not None else "NA",
        ) for i, values in enumerate(matrix, 1) for j, p in enumerate(values, 1)]

    @staticmethod
    def write_input(path, rows):
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    def test_symmetric_prefixes_use_eligible_denominator_and_strict_threshold(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "pairs.tsv"
            self.write_input(path, self.rows([
                [0.005, 0.01, 0.2], [0.001, None, 0.2], [0.2, 0.2, 0.005],
            ]))
            region = summary.read_pairs(path, "pair_nominal_p", 0.01)["SCR1"]
        result, steps = summary.summarize_region(region)
        self.assertEqual(result["n_eligible_pairs"], 8)
        self.assertEqual(result["n_significant_pairs"], 3)
        self.assertEqual(result["overall_significant_pair_fraction"], 3 / 8)
        self.assertEqual([row["significant_pair_fraction"] for row in steps], [1, 2 / 3, 3 / 8])
        self.assertEqual(result["max_prefix_significant_pair_fraction"], 1)
        self.assertEqual((result["best_left_max_rank"], result["best_right_max_rank"]), (1, 1))

    def test_distal_significant_pair_cannot_be_selected_without_inner_pairs(self):
        region = self.region([[0.5, 0.5], [0.5, 0.001]])
        for mode in ("symmetric", "independent"):
            result, _ = summary.summarize_region(region, mode)
            self.assertEqual(result["max_prefix_significant_pair_fraction"], 0.25)

    def test_independent_extents_and_left_right_symmetry(self):
        matrix = [[0.5, 0.001, 0.001], [0.5, 0.5, 0.5], [0.5, 0.5, 0.5]]
        result, _ = summary.summarize_region(self.region(matrix))
        self.assertEqual(result["max_prefix_significant_pair_fraction"], 0.25)
        result, steps = summary.summarize_region(self.region(matrix), "independent")
        self.assertEqual(len(steps), 9)
        self.assertEqual(result["overall_significant_pair_fraction"], 2 / 9)
        self.assertEqual(result["max_prefix_significant_pair_fraction"], 2 / 3)
        self.assertEqual((result["best_left_max_rank"], result["best_right_max_rank"]), (1, 3))
        swapped, _ = summary.summarize_region(self.region(list(zip(*matrix))), "independent")
        self.assertEqual(swapped["max_prefix_significant_pair_fraction"], 2 / 3)
        self.assertEqual((swapped["best_left_max_rank"], swapped["best_right_max_rank"]), (3, 1))

    def test_minimum_pairs_only_restricts_maximum(self):
        region = self.region([[0.001, 0.3], [0.3, 0.3]])
        result, steps = summary.summarize_region(region, min_eligible_pairs=2)
        self.assertFalse(steps[0]["meets_min_eligible_pairs"])
        self.assertEqual(result["max_prefix_significant_pair_fraction"], 0.25)
        result, _ = summary.summarize_region(region, min_eligible_pairs=5)
        self.assertEqual(result["overall_significant_pair_fraction"], 0.25)
        self.assertIsNone(result["max_prefix_significant_pair_fraction"])
        self.assertEqual(result["status"], "no_prefix_meets_min_eligible_pairs")

    def test_zero_is_not_missing_and_no_eligible_is_na(self):
        result, _ = summary.summarize_region(self.region([[0.5, 0.5], [0.5, 0.5]]))
        self.assertEqual(result["max_prefix_significant_pair_fraction"], 0)
        self.assertEqual(result["best_n_eligible_pairs"], 4)
        self.assertEqual(result["n_maximizing_prefixes"], 2)
        result, steps = summary.summarize_region(self.region([[None, None], [None, None]]))
        self.assertIsNone(result["overall_significant_pair_fraction"])
        self.assertIsNone(result["max_prefix_significant_pair_fraction"])
        self.assertEqual(result["status"], "no_eligible_pairs")
        self.assertTrue(all(row["significant_pair_fraction"] is None for row in steps))

    def test_invalid_windows_do_not_renumber_ranks_and_ties_choose_closest(self):
        result, steps = summary.summarize_region(self.region([
            [None, None, None], [0.001, 0.5, None], [None, None, None],
        ]))
        self.assertIsNone(steps[0]["significant_pair_fraction"])
        self.assertEqual(result["best_left_max_rank"], 2)
        self.assertEqual(result["best_right_max_rank"], 2)
        self.assertEqual(result["n_maximizing_prefixes"], 2)

    def test_missing_eligible_p_is_error_not_nonsignificant(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "pairs.tsv"
            rows = self.rows([[0.001, 0.5]])
            rows[1]["pair_nominal_p"] = "NA"
            self.write_input(path, rows)
            with self.assertRaisesRegex(ValueError, "without --skip-pair-test"):
                summary.read_pairs(path, "pair_nominal_p", 0.01)

    def test_bad_p_and_rank_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "pairs.tsv"
            for field, value in (("pair_nominal_p", "1.5"), ("pair_nominal_p", "inf"),
                                 ("pair_nominal_p", "-0.1"), ("block1_rank", "0"),
                                 ("block1_rank", "1.5")):
                with self.subTest(field=field, value=value):
                    rows = self.rows([[0.001]])
                    rows[0][field] = value
                    self.write_input(path, rows)
                    with self.assertRaises(ValueError):
                        summary.read_pairs(path, "pair_nominal_p", 0.01)

    def test_reverse_orientation_normalized_and_duplicate_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "pairs.tsv"
            rows = self.rows([[0.001, 0.5]])
            reverse = dict(rows[1], block1_side="right", block2_side="left",
                           block1_rank="2", block2_rank="1")
            self.write_input(path, [rows[0], reverse])
            region = summary.read_pairs(path, "pair_nominal_p", 0.01)["SCR1"]
            self.assertEqual(set(region.pairs), {(1, 1), (1, 2)})
            self.write_input(path, rows + [reverse])
            with self.assertRaisesRegex(ValueError, "duplicate cross-flank pair"):
                summary.read_pairs(path, "pair_nominal_p", 0.01)

    def test_same_side_pairs_ignored_and_region_coordinates_checked(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "pairs.tsv"
            rows = self.rows([[0.001]])
            same_side = dict(rows[0], block2_side="left", block2_rank="2")
            self.write_input(path, rows + [same_side])
            region = summary.read_pairs(path, "pair_nominal_p", 0.01)["SCR1"]
            result, _ = summary.summarize_region(region)
            self.assertEqual(result["n_eligible_pairs"], 1)
            self.assertEqual(result["n_ignored_same_side_pairs"], 1)
            same_side["region_end"] = "30000"
            self.write_input(path, rows + [same_side])
            with self.assertRaisesRegex(ValueError, "inconsistent coordinates"):
                summary.read_pairs(path, "pair_nominal_p", 0.01)

    def test_matched_only_input_is_not_a_complete_grid(self):
        region = self.region([[0.001, 0.5], [0.5, 0.5]])
        del region.pairs[(1, 2)]
        del region.pairs[(2, 1)]
        with self.assertRaisesRegex(ValueError, "incomplete cross-flank grid"):
            summary.summarize_region(region)
        result, _ = summary.summarize_region(region, allow_incomplete_grid=True)
        self.assertEqual(result["n_unreported_pairs"], 2)

    def test_counts_agree_with_brute_force_for_random_rectangular_grids(self):
        rng = random.Random(193)
        for _ in range(30):
            left, right = rng.randint(1, 8), rng.randint(1, 8)
            region = self.region([[rng.choice([None, 0.005, 0.01, 0.5])
                                   for _ in range(right)] for _ in range(left)])
            for mode in ("symmetric", "independent"):
                result, steps = summary.summarize_region(region, mode)
                fractions = []
                for step in steps:
                    pairs = [value for (i, j), value in region.pairs.items()
                             if i <= step["left_max_rank"] and j <= step["right_max_rank"]]
                    n, s = sum(value[0] for value in pairs), sum(value[1] for value in pairs)
                    self.assertEqual((step["n_eligible_pairs"], step["n_significant_pairs"]), (n, s))
                    if n:
                        fractions.append(s / n)
                self.assertEqual(result["max_prefix_significant_pair_fraction"],
                                 max(fractions) if fractions else None)
                self.assertEqual(steps[-1]["n_eligible_pairs"], result["n_eligible_pairs"])

    def test_cli_multiregion_gzip_output_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            path, prefix = Path(tmp) / "pairs.tsv.gz", Path(tmp) / "out"
            self.write_input(path, self.rows([[0.001, 0.5], [None, 0.01]], "one")
                             + self.rows([[None]], "two"))
            command = [sys.executable, str(Path(summary.__file__).resolve()), str(path), "-o", str(prefix)]
            completed = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            output = Path(str(prefix) + ".summary.tsv")
            with output.open() as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                self.assertEqual(reader.fieldnames[:2], list(summary.SUMMARY_FIELDS[:2]))
                rows = list(reader)
            self.assertAlmostEqual(float(rows[0]["overall_significant_pair_fraction"]), 1 / 3)
            self.assertEqual(rows[0]["max_prefix_significant_pair_fraction"], "1")
            self.assertEqual(rows[1]["max_prefix_significant_pair_fraction"], "NA")
            previous = output.read_bytes()
            repeated = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(repeated.returncode, 2)
            self.assertEqual(previous, output.read_bytes())

    def test_actual_background_output_matches_existing_overall_fraction(self):
        from flanking_block_ld import background_block_ld as workflow
        from flanking_block_ld.tests.test_background_block_ld import BackgroundBlockLDTests
        from flanking_block_ld.tests.test_background_region_pair_p import BackgroundRegionPairPTests

        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            BackgroundRegionPairPTests.write_inputs(directory)
            args = BackgroundBlockLDTests.args(directory, "background")
            args.skip_region_test = True
            args.pair_significance_threshold = 0.3
            workflow.run(args)
            parsed = summary.build_parser().parse_args([
                str(directory / "background.target_pairs.tsv"), "-o", str(directory / "fractions"),
                "--p-threshold", "0.3",
            ])
            summary.run(parsed)
            with (directory / "background.summary.tsv").open() as handle:
                upstream = next(csv.DictReader(handle, delimiter="\t"))
            with (directory / "fractions.summary.tsv").open() as handle:
                observed = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(observed["n_eligible_pairs"], upstream["n_valid_pairs"])
            self.assertEqual(observed["n_significant_pairs"], upstream["n_significant_pairs"])
            self.assertEqual(float(observed["overall_significant_pair_fraction"]),
                             float(upstream["significant_pair_fraction"]))


if __name__ == "__main__":
    unittest.main()

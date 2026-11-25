"""Checks for independent cumulative-window maximum median r²."""

import csv
import statistics
import tempfile
import unittest
from pathlib import Path

from flanking_block_ld.flanking_block_ld import Region, VCFData, Variant, pairwise_r2
from flanking_block_ld.max_prefix_median_r2 import calculate, main, match_profile, Profile
from flanking_block_ld.median_flanking_r2 import position_index


class MaxPrefixMedianR2Test(unittest.TestCase):
    def test_short_flank_reports_descriptive_max_without_background_test(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "target.bed").write_text("chr1\t20\t40\ttarget\n")
            (root / "background.bed").write_text("chr1\t300\t350\n")
            header = ("##fileformat=VCFv4.2\n"
                      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\tD\n")
            row = "chr1\t{}\t.\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|1\t1|0\t1|1\n"
            (root / "variants.vcf").write_text(header + row.format(11) + row.format(41))
            prefix = root / "out"
            code = main([str(root / "target.bed"), str(root / "background.bed"),
                         str(root / "variants.vcf"), "--output-prefix", str(prefix),
                         "--flank-bp", "40", "--window-bp", "10", "--attempts", "20",
                         "--backgrounds", "5", "--length-tolerance-bp", "0",
                         "--vcf-access", "preload"])
            self.assertEqual(code, 0)
            with Path(str(prefix) + ".summary.tsv").open() as handle:
                summary = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(summary["status"], "target_short_flank")
            self.assertEqual(float(summary["max_prefix_median_r2"]), 1.0)
            self.assertEqual(summary["max_search_qualified"], "False")
            self.assertEqual(summary["empirical_p_upper"], "NA")

    def test_prefix_medians_match_direct_calculation(self):
        keys = [(str(i), 0) for i in range(8)]
        patterns = ([0, 0, 0, 0, 1, 1, 1, 1], [0, 0, 1, 1, 0, 0, 1, 1])
        region = Region("1", "chr1", 100, 120, "target", 0)
        positions = (65, 75, 85, 95, 120, 130, 140, 150)
        variants = tuple(Variant("1", pos + 1, str(pos), "A", "G", 0.5, 1,
                                 dict(zip(keys, patterns[index % 2])), index)
                         for index, pos in enumerate(positions))
        data = VCFData(8, {"1": variants}, {}, {}, 0, 0, 0, 0, 0, 0)
        result = calculate(region, data, position_index(data), 40, 10, 4, 1, 1)
        direct = []
        for left_rank in range(1, 5):
            for right_rank in range(1, 5):
                left = [v for v in variants if 100 - left_rank * 10 <= v.bed_pos < 100]
                right = [v for v in variants if 120 <= v.bed_pos < 120 + right_rank * 10]
                direct.append(statistics.median(pairwise_r2(a, b, keys) for a in left for b in right))
        self.assertEqual([step.median_r2 for step in result.steps], direct)
        self.assertEqual(result.maximum, max(direct))
        self.assertEqual(len(result.steps), 16)
        self.assertIn((result.best_left_rank, result.best_right_rank),
                      [(step.left_rank, step.right_rank) for step in result.steps
                       if step.median_r2 == result.maximum])
        too_few = calculate(region, data, position_index(data), 40, 10, 4, 1, 5)
        self.assertIsNotNone(too_few.maximum)
        self.assertFalse(too_few.threshold_met)

    def test_match_accepts_swapped_window_counts(self):
        region = Region("1", "chr1", 100, 120, "target", 0)
        target = Profile(region, (1, 2, 0, 3), (2, 0, 1, 0))
        candidate = Profile(region, target.right_counts, target.left_counts)
        match = match_profile(candidate, target, 0)
        self.assertTrue(match.strict)
        self.assertEqual(match.orientation, "swapped")

    def test_end_to_end_background_uses_same_max_search(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "target.bed").write_text("chr1\t100\t120\ttarget\n")
            (root / "background.bed").write_text("chr1\t300\t1300\n")
            header = ("##fileformat=VCFv4.2\n"
                      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\tD\n")
            row = "chr1\t{}\t.\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|1\t1|0\t1|1\n"
            (root / "variants.vcf").write_text(header + "".join(row.format(p) for p in range(1, 1301, 10)))
            prefix = root / "out"
            code = main([str(root / "target.bed"), str(root / "background.bed"),
                         str(root / "variants.vcf"), "--output-prefix", str(prefix),
                         "--flank-bp", "40", "--window-bp", "10", "--attempts", "50",
                         "--backgrounds", "5", "--length-tolerance-bp", "10",
                         "--min-snps-per-side", "2",
                         "--vcf-access", "preload"])
            self.assertEqual(code, 0)
            with Path(str(prefix) + ".profile.tsv").open() as handle:
                profile = list(csv.DictReader(handle, delimiter="\t"))
            with Path(str(prefix) + ".summary.tsv").open() as handle:
                summary = next(csv.DictReader(handle, delimiter="\t"))
            with Path(str(prefix) + ".backgrounds.tsv").open() as handle:
                backgrounds = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(profile), 16)
            self.assertEqual(int(summary["valid_backgrounds"]), 5)
            self.assertEqual(len(backgrounds), 5)
            self.assertEqual(float(summary["max_prefix_median_r2"]),
                             max(float(row["median_r2"]) for row in profile if row["median_r2"] != "NA"))
            self.assertEqual(float(profile[0]["median_r2"]), 1.0)
            self.assertEqual(profile[0]["eligible_for_max"], "False")
            self.assertIn("best_left_rank", summary)
            self.assertIn("best_right_rank", summary)
            self.assertNotEqual(summary["best_left_rank"], summary["best_right_rank"])


if __name__ == "__main__":
    unittest.main()

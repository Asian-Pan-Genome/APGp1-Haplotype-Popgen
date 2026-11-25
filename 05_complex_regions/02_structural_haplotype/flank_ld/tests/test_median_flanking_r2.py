"""Small end-to-end checks for median cross-flank r² sampling."""

import csv
import tempfile
import unittest
from pathlib import Path

from flanking_block_ld.flanking_block_ld import Region
from flanking_block_ld.max_prefix_median_r2 import main as max_prefix_main
from flanking_block_ld.median_flanking_r2 import Candidate, main, match_details


class MedianFlankingR2Test(unittest.TestCase):
    def test_background_flanks_extend_beyond_bed_when_requested(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "target.bed").write_text("chr1\t100\t120\ttarget\n")
            (root / "background.bed").write_text("chr1\t300\t350\n")
            header = ("##fileformat=VCFv4.2\n"
                      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\tD\n")
            row = "chr1\t{}\t.\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|1\t1|0\t1|1\n"
            (root / "variants.vcf").write_text(
                header + "".join(row.format(position) for position in range(1, 501, 5)))
            shared = [str(root / "target.bed"), str(root / "background.bed"),
                      str(root / "variants.vcf"), "--flank-bp", "40",
                      "--attempts", "100", "--backgrounds", "1",
                      "--length-tolerance-bp", "0", "--vcf-access", "preload",
                      "--allow-background-flank-extension"]

            median_prefix = root / "median"
            self.assertEqual(main([*shared, "--output-prefix", str(median_prefix)]), 0)
            with Path(str(median_prefix) + ".summary.tsv").open() as handle:
                median_summary = next(csv.DictReader(handle, delimiter="\t"))
            with Path(str(median_prefix) + ".backgrounds.tsv").open() as handle:
                median_background = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(median_summary["status"], "ok")
            self.assertNotEqual(median_summary["empirical_p_upper"], "NA")
            self.assertEqual(median_background["extended_beyond_background_bed"], "True")
            self.assertEqual(median_background["full_flanks"], "True")
            self.assertEqual(median_background["left_snps"], "8")
            self.assertEqual(median_background["right_snps"], "8")

            max_prefix = root / "max"
            self.assertEqual(max_prefix_main([
                *shared, "--output-prefix", str(max_prefix), "--window-bp", "10",
                "--min-pairs", "1", "--min-snps-per-side", "5"]), 0)
            with Path(str(max_prefix) + ".summary.tsv").open() as handle:
                max_summary = next(csv.DictReader(handle, delimiter="\t"))
            with Path(str(max_prefix) + ".backgrounds.tsv").open() as handle:
                max_background = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(max_summary["status"], "ok")
            self.assertNotEqual(max_summary["empirical_p_upper"], "NA")
            self.assertEqual(max_background["extended_beyond_background_bed"], "True")
            self.assertEqual(max_background["full_flanks"], "True")
            self.assertEqual(max_background["left_snps"], "8")
            self.assertEqual(max_background["right_snps"], "8")

    def test_extension_allows_p_with_background_short_at_chromosome_start(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "target.bed").write_text("chr2\t100\t120\ttarget\n")
            (root / "background.bed").write_text("chr1\t0\t50\n")
            header = ("##fileformat=VCFv4.2\n"
                      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\tD\n")
            row = "{}\t{}\t.\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|1\t1|0\t1|1\n"
            (root / "variants.vcf").write_text(
                header +
                "".join(row.format("chr1", position) for position in range(1, 96, 5)) +
                "".join(row.format("chr2", position) for position in range(61, 161, 5)))
            shared = [str(root / "target.bed"), str(root / "background.bed"),
                      str(root / "variants.vcf"), "--flank-bp", "40",
                      "--attempts", "100", "--backgrounds", "1",
                      "--length-tolerance-bp", "0", "--vcf-access", "preload",
                      "--allow-background-flank-extension"]
            for entry, extra in ((main, []),
                                 (max_prefix_main, ["--window-bp", "10",
                                                    "--min-pairs", "1",
                                                    "--min-snps-per-side", "1"])):
                prefix = root / ("median" if entry is main else "max")
                self.assertEqual(entry([*shared, "--output-prefix", str(prefix), *extra]), 0)
                with Path(str(prefix) + ".summary.tsv").open() as handle:
                    summary = next(csv.DictReader(handle, delimiter="\t"))
                with Path(str(prefix) + ".backgrounds.tsv").open() as handle:
                    background = next(csv.DictReader(handle, delimiter="\t"))
                self.assertEqual(summary["status"], "ok")
                self.assertNotEqual(summary["empirical_p_upper"], "NA")
                self.assertEqual(background["full_flanks"], "False")
                self.assertLess(int(background["left_flank_bp"]), 40)


    def test_partial_background_median_is_reported_without_p_value(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "target.bed").write_text("chr1\t100\t120\ttarget\n")
            (root / "background.bed").write_text("chr1\t300\t350\n")
            header = ("##fileformat=VCFv4.2\n"
                      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\tD\n")
            row = "chr1\t{}\t.\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|1\t1|0\t1|1\n"
            positions = [61, 71, 81, 91, 121, 131, 141, 151] + list(range(301, 351, 5))
            (root / "variants.vcf").write_text(header + "".join(row.format(p) for p in positions))
            prefix = root / "out"
            code = main([str(root / "target.bed"), str(root / "background.bed"),
                         str(root / "variants.vcf"), "--output-prefix", str(prefix),
                         "--flank-bp", "40", "--attempts", "100", "--backgrounds", "1",
                         "--length-tolerance-bp", "0", "--vcf-access", "preload"])
            self.assertEqual(code, 0)
            with Path(str(prefix) + ".summary.tsv").open() as handle:
                summary = next(csv.DictReader(handle, delimiter="\t"))
            with Path(str(prefix) + ".backgrounds.tsv").open() as handle:
                background = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(summary["status"], "background_short_flanks")
            self.assertEqual(summary["empirical_p_upper"], "NA")
            self.assertEqual(background["full_flanks"], "False")
            self.assertNotEqual(background["median_r2"], "NA")
            max_prefix = root / "max_out"
            max_code = max_prefix_main([
                str(root / "target.bed"), str(root / "background.bed"),
                str(root / "variants.vcf"), "--output-prefix", str(max_prefix),
                "--flank-bp", "40", "--window-bp", "10", "--attempts", "100",
                "--backgrounds", "1", "--length-tolerance-bp", "0",
                "--min-pairs", "1", "--min-snps-per-side", "1",
                "--vcf-access", "preload"])
            self.assertEqual(max_code, 0)
            with Path(str(max_prefix) + ".summary.tsv").open() as handle:
                max_summary = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(max_summary["status"], "background_short_flanks")
            self.assertNotEqual(max_summary["max_prefix_median_r2"], "NA")
            self.assertEqual(max_summary["empirical_p_upper"], "NA")

    def test_short_flank_and_no_full_length_background_still_reports_target(self):
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
                         "--flank-bp", "40", "--attempts", "20", "--backgrounds", "5",
                         "--length-tolerance-bp", "0", "--vcf-access", "preload"])
            self.assertEqual(code, 0)
            with Path(str(prefix) + ".summary.tsv").open() as handle:
                summary = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(summary["status"], "target_short_flank")
            self.assertEqual(float(summary["median_r2"]), 1.0)
            self.assertEqual(summary["empirical_p_upper"], "NA")

    def test_snp_counts_match_after_left_right_swap(self):
        target = Candidate(Region("1", "chr1", 100, 120, "target", 0), 4, 2)
        background = Candidate(Region("1", "chr1", 300, 320, "background", 0), 2, 4)
        strict, score, orientation = match_details(background, target, 0, 10_000)
        self.assertTrue(strict)
        self.assertEqual(score, 0)
        self.assertEqual(orientation, "swapped")

    def test_preload_selects_valid_backgrounds(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "target.bed").write_text("chr1\t100\t120\ttarget\n")
            (root / "background.bed").write_text("chr1\t300\t1300\n")
            header = ("##fileformat=VCFv4.2\n"
                      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                      "A\tB\tC\tD\n")
            rows = []
            for position in range(1, 1301, 10):
                pattern = ("0|0", "0|1", "1|0", "1|1") if position % 20 == 1 else (
                    "0|1", "1|0", "0|0", "1|1")
                rows.append("\t".join(("chr1", str(position), ".", "A", "G", ".", "PASS", ".",
                                       "GT", *pattern)) + "\n")
            (root / "variants.vcf").write_text(header + "".join(rows))
            prefix = root / "out"
            result = main([str(root / "target.bed"), str(root / "background.bed"),
                           str(root / "variants.vcf"), "--output-prefix", str(prefix),
                           "--flank-bp", "20", "--attempts", "50", "--backgrounds", "10",
                           "--length-tolerance-bp", "10", "--vcf-access", "preload"])
            self.assertEqual(result, 0)
            with Path(str(prefix) + ".summary.tsv").open() as handle:
                summary = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(summary["valid_backgrounds"], "10")
            self.assertEqual(summary["left_snps"], "2")
            self.assertEqual(summary["right_snps"], "2")
            self.assertEqual(float(summary["median_r2"]), 0.5)
            self.assertEqual(float(summary["max_r2"]), 1.0)
            self.assertIn("empirical_p_upper_max_r2", summary)
            with Path(str(prefix) + ".backgrounds.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 10)
            self.assertTrue(all(float(row["max_r2"]) >= float(row["median_r2"]) for row in rows))
            self.assertTrue(all(300 <= int(row["start"]) - 20 and
                                int(row["end"]) + 20 <= 1300 for row in rows))
            self.assertTrue(all(abs(int(row["length_bp"]) - 20) <= 10 for row in rows))

    def test_fallback_fills_when_snp_counts_differ(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "target.bed").write_text("chr1\t100\t120\ttarget\n")
            (root / "background.bed").write_text("chr1\t300\t1100\n")
            header = ("##fileformat=VCFv4.2\n"
                      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                      "A\tB\tC\tD\n")
            target_positions = [61, 71, 81, 91, 121, 131, 141, 151]
            background_positions = list(range(301, 1101, 20))
            row = "chr1\t{}\t.\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|1\t1|0\t1|1\n"
            (root / "variants.vcf").write_text(
                header + "".join(row.format(pos) for pos in target_positions + background_positions))
            prefix = root / "out"
            result = main([str(root / "target.bed"), str(root / "background.bed"),
                           str(root / "variants.vcf"), "--output-prefix", str(prefix),
                           "--flank-bp", "40", "--attempts", "100", "--backgrounds", "5",
                           "--length-tolerance-bp", "0", "--snp-count-tolerance", "0",
                           "--vcf-access", "preload"])
            self.assertEqual(result, 0)
            with Path(str(prefix) + ".summary.tsv").open() as handle:
                summary = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(summary["strict_available"], "0")
            self.assertEqual(summary["selected_fallback"], "5")
            self.assertEqual(summary["valid_backgrounds"], "5")


if __name__ == "__main__":
    unittest.main()

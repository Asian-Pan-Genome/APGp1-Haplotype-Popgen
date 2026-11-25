import argparse
import tempfile
import unittest
from pathlib import Path

from flanking_block_ld.fixed_window_block_ld import run


class FixedWindowBlockLDTests(unittest.TestCase):
    def write_inputs(self, directory: Path):
        bed = directory / "regions.bed"
        bed.write_text("chr1\t1000\t2000\ttarget\n", encoding="utf-8")
        vcf = directory / "cohort.vcf"
        rows = [
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4",
        ]
        for pos, variant_id in (
            (850, "L2"),
            (950, "L1"),
            (2050, "R1"),
            (2150, "R2"),
        ):
            rows.append(
                f"1\t{pos}\t{variant_id}\tA\tG\t.\tPASS\t.\tGT\t"
                "0|1\t0|1\t0|1\t0|1"
            )
        vcf.write_text("\n".join(rows) + "\n", encoding="utf-8")
        return bed, vcf

    def args(self, bed: Path, vcf: Path, output: Path, **overrides):
        values = dict(
            regions_bed=bed,
            vcf=vcf,
            output=output,
            window_bp=100,
            window_step_bp=None,
            block_count=2,
            max_flank_bp=200,
            pairing="matched",
            min_snps_per_block=1,
            max_snps_per_block=10,
            max_exact_entropy_combinations=100_000,
            min_maf=0.1,
            max_site_missing_rate=0.05,
            max_window_missing_filtered_fraction=0.10,
            permutations=0,
            seed=1,
            include_filtered=False,
            allow_unphased_missing=False,
        )
        values.update(overrides)
        return argparse.Namespace(**values)

    def read_rows(self, path: Path):
        lines = path.read_text(encoding="utf-8").splitlines()
        header = lines[0].split("\t")
        return [dict(zip(header, line.split("\t"))) for line in lines[1:]]

    def test_matched_pairs_are_emitted_near_to_far(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            bed, vcf = self.write_inputs(directory)
            output = directory / "matched.tsv"
            summary = run(self.args(bed, vcf, output))
            rows = self.read_rows(output)
            header = output.read_text(encoding="utf-8").splitlines()[0].split("\t")

            self.assertEqual(summary["block_pair_count"], 2)
            self.assertEqual(header[7], "multiallelic_dprime")
            self.assertEqual(
                header[8:14],
                [
                    "block1_snp_number_pass_threshold",
                    "block2_snp_number_pass_threshold",
                    "both_blocks_pass_snp_number_threshold",
                    "block1_missed_snp_number",
                    "block2_missed_snp_number",
                    "total_missed_snp_number",
                ],
            )
            self.assertEqual(header[14:17], [
                "dprime_excess",
                "dprime_permutation_mean",
                "dprime_empirical_p",
            ])
            self.assertEqual(rows[0]["block1_snp_number_pass_threshold"], "1")
            self.assertEqual(rows[0]["block2_snp_number_pass_threshold"], "1")
            self.assertEqual(rows[0]["both_blocks_pass_snp_number_threshold"], "true")
            self.assertEqual(rows[0]["total_missed_snp_number"], "0")
            self.assertEqual([row["pair_rank_near_to_far"] for row in rows], ["1", "2"])
            self.assertEqual(rows[0]["left_snps"], "L1@950")
            self.assertEqual(rows[0]["right_snps"], "R1@2050")
            self.assertEqual(rows[0]["left_window_start"], "900")
            self.assertEqual(rows[0]["left_window_end"], "1000")
            self.assertEqual(rows[1]["left_snps"], "L2@850")
            self.assertEqual(rows[1]["right_snps"], "R2@2150")
            self.assertAlmostEqual(float(rows[0]["multiallelic_dprime"]), 1.0)
            self.assertAlmostEqual(float(rows[1]["multiallelic_dprime"]), 1.0)
            self.assertTrue(all(row["status"] == "ok" for row in rows))

    def test_all_pairing_emits_full_cross_window_matrix(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            bed, vcf = self.write_inputs(directory)
            output = directory / "all.tsv"
            summary = run(self.args(bed, vcf, output, pairing="all"))
            rows = self.read_rows(output)

            self.assertEqual(summary["block_pair_count"], 4)
            pairs = [(row["left_block_rank"], row["right_block_rank"]) for row in rows]
            self.assertEqual(pairs, [("1", "1"), ("1", "2"), ("2", "1"), ("2", "2")])
            self.assertEqual(rows[1]["pair_rank_near_to_far"], ".")

    def test_all_blocks_pairing_emits_full_genomic_triangle(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            bed, vcf = self.write_inputs(directory)
            output = directory / "all_blocks.tsv"
            summary = run(self.args(bed, vcf, output, pairing="all-blocks"))
            rows = self.read_rows(output)

            self.assertEqual(summary["block_pair_count"], 6)
            self.assertEqual(
                [(row["block1_label"], row["block2_label"]) for row in rows],
                [
                    ("L2", "L1"),
                    ("L2", "R1"),
                    ("L2", "R2"),
                    ("L1", "R1"),
                    ("L1", "R2"),
                    ("R1", "R2"),
                ],
            )
            self.assertEqual(
                [
                    (row["block1_genomic_index"], row["block2_genomic_index"])
                    for row in rows
                ],
                [("1", "2"), ("1", "3"), ("1", "4"), ("2", "3"), ("2", "4"), ("3", "4")],
            )
            self.assertEqual(
                [row["pair_type"] for row in rows],
                [
                    "left_left",
                    "left_right",
                    "left_right",
                    "left_right",
                    "left_right",
                    "right_right",
                ],
            )

    def test_window_step_creates_overlapping_sliding_windows(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            bed, vcf = self.write_inputs(directory)
            output = directory / "sliding.tsv"
            run(
                self.args(
                    bed,
                    vcf,
                    output,
                    window_step_bp=50,
                    pairing="matched",
                )
            )
            rows = self.read_rows(output)

            self.assertEqual(
                [
                    (row["left_window_start"], row["left_window_end"])
                    for row in rows
                ],
                [("900", "1000"), ("850", "950")],
            )
            self.assertEqual(
                [
                    (row["right_window_start"], row["right_window_end"])
                    for row in rows
                ],
                [("2000", "2100"), ("2050", "2150")],
            )
            self.assertEqual(rows[0]["left_snps"], "L1@950")
            self.assertEqual(rows[1]["left_snps"], "L1@950")

    def test_window_is_filtered_when_too_many_sites_fail_missing_qc(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            bed = directory / "regions.bed"
            bed.write_text("1\t1000\t2000\ttarget\n", encoding="utf-8")
            vcf = directory / "missing.vcf"
            vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n"
                "1\t950\tLgood\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|1\t0|1\t0|1\n"
                "1\t960\tLmissing\tA\tG\t.\tPASS\t.\tGT\t.|.\t0|1\t0|1\t0|1\n"
                "1\t2050\tRgood\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|1\t0|1\t0|1\n",
                encoding="utf-8",
            )
            output = directory / "missing.tsv"
            run(
                self.args(
                    bed,
                    vcf,
                    output,
                    block_count=1,
                    max_flank_bp=100,
                )
            )
            row = self.read_rows(output)[0]

            self.assertEqual(row["left_pre_missing_snp_count"], "2")
            self.assertEqual(row["left_missing_filtered_snp_count"], "1")
            self.assertEqual(row["block1_missed_snp_number"], "1")
            self.assertEqual(row["block2_missed_snp_number"], "0")
            self.assertEqual(row["total_missed_snp_number"], "1")
            self.assertEqual(row["block1_snp_number_pass_threshold"], "0")
            self.assertEqual(row["left_missing_filtered_fraction"], "0.5")
            self.assertEqual(row["left_snp_count"], "0")
            self.assertEqual(row["multiallelic_dprime"], "NA")
            self.assertIn("left_window_excess_missing_filtered_sites", row["status"])

    def test_excess_snps_are_capped_by_exact_maximum_entropy(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            bed = directory / "regions.bed"
            bed.write_text("1\t1000\t2000\ttarget\n", encoding="utf-8")
            vcf = directory / "entropy.vcf"
            vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n"
                "1\t970\tB\tA\tG\t.\tPASS\t.\tGT\t0|0\t1|1\t0|0\t1|1\n"
                "1\t980\tAdup\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|0\t1|1\t1|1\n"
                "1\t990\tA\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|0\t1|1\t1|1\n"
                "1\t2050\tR\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|1\t0|1\t0|1\n",
                encoding="utf-8",
            )
            output = directory / "entropy.tsv"
            run(
                self.args(
                    bed,
                    vcf,
                    output,
                    block_count=1,
                    max_flank_bp=100,
                    max_snps_per_block=2,
                )
            )
            row = self.read_rows(output)[0]

            self.assertEqual(row["left_snp_count_before_cap"], "3")
            self.assertEqual(row["left_snp_count"], "2")
            self.assertEqual(row["left_snp_selection_method"], "exact_max_entropy")
            self.assertEqual(row["left_snps"], "A@990,B@970")

    def test_window_below_minimum_snp_count_is_filtered(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            bed, vcf = self.write_inputs(directory)
            output = directory / "minimum.tsv"
            run(self.args(bed, vcf, output, min_snps_per_block=2))
            rows = self.read_rows(output)
            self.assertTrue(all(row["multiallelic_dprime"] == "NA" for row in rows))
            self.assertTrue(all("insufficient_left_snps" in row["status"] for row in rows))
            self.assertTrue(all("insufficient_right_snps" in row["status"] for row in rows))


if __name__ == "__main__":
    unittest.main()

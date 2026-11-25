import argparse
import gzip
import tempfile
import unittest
from pathlib import Path
from typing import Optional

from flanking_block_ld.flanking_block_ld import multiallelic_dprime, read_regions, run


class FlankingBlockLDTests(unittest.TestCase):
    def write_bed(self, directory: Path) -> Path:
        path = directory / "regions.bed"
        path.write_text("chr1\t1000\t2000\ttarget\n", encoding="utf-8")
        return path

    def write_vcf(
        self,
        directory: Path,
        gzip_output: bool = False,
        unphased_variant: Optional[str] = None,
    ) -> Path:
        positions = [
            (801, "far_left"),
            (901, "L2"),
            (1000, "L1"),
            (1500, "inside"),
            (2001, "R1"),
            (2101, "R2"),
            (2201, "far_right"),
        ]
        lines = [
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4",
        ]
        for pos, variant_id in positions:
            gt = "0/1" if variant_id == unphased_variant else "0|1"
            lines.append(
                f"1\t{pos}\t{variant_id}\tA\tG\t.\tPASS\t.\tGT\t"
                + "\t".join([gt] * 4)
            )
        content = "\n".join(lines) + "\n"
        if gzip_output:
            path = directory / "cohort.vcf.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write(content)
        else:
            path = directory / "cohort.vcf"
            path.write_text(content, encoding="utf-8")
        return path

    def args(self, bed: Path, vcf: Path, output: Path, **overrides):
        values = dict(
            regions_bed=bed,
            vcf=vcf,
            output=output,
            block_snps=2,
            flank_bp=150,
            min_maf=0.1,
            min_call_rate=0.9,
            selection="nearest",
            permutations=0,
            seed=1,
            include_filtered=False,
            allow_unphased_missing=False,
            allow_partial_blocks=False,
        )
        values.update(overrides)
        return argparse.Namespace(**values)

    def read_output_row(self, path: Path):
        lines = path.read_text(encoding="utf-8").splitlines()
        return dict(zip(lines[0].split("\t"), lines[1].split("\t")))

    def test_multiallelic_dprime_known_cases(self):
        left = [(0,), (0,), (1,), (1,)]
        perfectly_associated = [(0,), (0,), (1,), (1,)]
        independent = [(0,), (1,), (0,), (1,)]
        self.assertAlmostEqual(multiallelic_dprime(left, perfectly_associated), 1.0)
        self.assertAlmostEqual(multiallelic_dprime(left, independent), 0.0)

    def test_phased_vcf_calculates_paper_dprime_and_rv(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            bed = self.write_bed(directory)
            vcf = self.write_vcf(directory, gzip_output=True)
            output = directory / "summary.tsv"
            summary = run(self.args(bed, vcf, output, permutations=100))
            row = self.read_output_row(output)

            self.assertEqual(summary["region_count"], 1)
            self.assertEqual(summary["ok_region_count"], 1)
            self.assertEqual(summary["sample_count"], 4)
            self.assertEqual(row["left_snps_nearest_or_selected_order"], "L1@1000,L2@901")
            self.assertEqual(row["right_snps_nearest_or_selected_order"], "R1@2001,R2@2101")
            self.assertEqual(row["complete_haplotype_count"], "8")
            self.assertEqual(row["left_haplotype_count"], "2")
            self.assertAlmostEqual(float(row["multiallelic_dprime"]), 1.0)
            self.assertAlmostEqual(float(row["cross_mean_r2"]), 1.0)
            self.assertAlmostEqual(float(row["block_rv_r2"]), 1.0)
            self.assertNotEqual(row["dprime_empirical_p"], "NA")
            self.assertEqual(row["status"], "ok")

    def test_unphased_heterozygous_variant_is_filtered_by_default(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            bed = self.write_bed(directory)
            vcf = self.write_vcf(directory, unphased_variant="L1")
            output = directory / "summary.tsv"
            summary = run(self.args(bed, vcf, output))
            row = self.read_output_row(output)

            self.assertEqual(summary["filtered_unphased_count"], 1)
            self.assertEqual(summary["ok_region_count"], 0)
            self.assertEqual(row["left_snp_count"], "1")
            self.assertEqual(row["multiallelic_dprime"], "NA")
            self.assertIn("insufficient_left_snps", row["status"])

    def test_bed_coordinates_are_preserved(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            bed = directory / "regions.bed"
            bed.write_text("chr1\t1000\t2000\n", encoding="utf-8")
            region = read_regions(bed)[0]
            self.assertEqual(region.chrom, "1")
            self.assertEqual(region.display_chrom, "chr1")
            self.assertEqual(region.region_id, "chr1:1000-2000")


if __name__ == "__main__":
    unittest.main()

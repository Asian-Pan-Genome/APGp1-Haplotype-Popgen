import argparse
import tempfile
import unittest
from pathlib import Path

from flanking_block_ld.entropy_capture_pilot import run


class EntropyCapturePilotTests(unittest.TestCase):
    def read_rows(self, path: Path):
        lines = path.read_text(encoding="utf-8").splitlines()
        header = lines[0].split("\t")
        return [dict(zip(header, line.split("\t"))) for line in lines[1:]]

    def test_reports_entropy_capture_and_sparse_window_fraction(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            bed = directory / "universe.bed"
            bed.write_text("1\t0\t400\tautosome_callable\n", encoding="utf-8")
            vcf = directory / "cohort.vcf"
            vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n"
                "1\t10\tA\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|0\t1|1\t1|1\n"
                "1\t20\tB\tA\tG\t.\tPASS\t.\tGT\t0|0\t1|1\t0|0\t1|1\n"
                "1\t30\tBdup\tA\tG\t.\tPASS\t.\tGT\t0|0\t1|1\t0|0\t1|1\n"
                "1\t110\tC\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|1\t0|1\t0|1\n"
                "1\t210\tD\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|0\t1|1\t1|1\n"
                "1\t220\tE\tA\tG\t.\tPASS\t.\tGT\t0|0\t1|1\t0|0\t1|1\n",
                encoding="utf-8",
            )
            prefix = directory / "pilot"
            args = argparse.Namespace(
                sampling_bed=bed,
                vcf=vcf,
                output_prefix=prefix,
                sample_windows=0,
                window_bp=100,
                sampling_step_bp=None,
                snp_counts=(2,),
                max_exact_entropy_combinations=100_000,
                min_complete_haplotypes=2,
                min_maf=0.1,
                max_site_missing_rate=0.05,
                max_window_missing_filtered_fraction=0.10,
                seed=1,
                include_filtered=False,
                allow_unphased_missing=False,
            )

            result = run(args)
            details = self.read_rows(result["detail_path"])
            summary = self.read_rows(result["summary_path"])[0]

            self.assertEqual(result["candidate_window_count"], 4)
            self.assertEqual(result["sampled_window_count"], 4)
            self.assertEqual(
                [row["eligible_snp_number"] for row in details],
                ["3", "1", "2", "0"],
            )
            self.assertEqual(summary["sparse_window_number_all"], "2")
            self.assertEqual(summary["sparse_window_fraction_all"], "0.5")
            self.assertEqual(summary["entropy_evaluable_window_number"], "2")
            self.assertAlmostEqual(float(summary["entropy_capture_median"]), 1.0)
            self.assertEqual(details[0]["selection_method"], "exact_max_entropy")
            self.assertAlmostEqual(float(details[0]["entropy_capture"]), 1.0)


if __name__ == "__main__":
    unittest.main()

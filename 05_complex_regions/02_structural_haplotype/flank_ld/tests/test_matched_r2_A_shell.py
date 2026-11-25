"""Integration check for equal-target-weight A aggregation."""

import csv
import subprocess
import tempfile
import unittest
from pathlib import Path


class MatchedR2AShellTest(unittest.TestCase):
    def test_ties_and_equal_target_weights(self):
        script = Path(__file__).resolve().parents[1] / "summarize_matched_r2_A.sh"
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for identifier, target, backgrounds in (
                ("a", 0.8, (0.2, 0.8, 1.0)),
                ("b", 0.9, (0.1, 0.2, 0.3)),
            ):
                prefix = root / (identifier + "_flanking_r2")
                Path(str(prefix) + ".summary.tsv").write_text(
                    "\r\n\ufeffmedian_r2\tstatus\r\n\r\n" + f"{target}\tok\r\n")
                if identifier == "a":
                    background_text = ("\r\n\ufefffull_flanks\tmedian_r2\r\n" +
                                       "".join(f"True\t{value}\r\n\r\n" for value in backgrounds))
                else:
                    background_text = ("\r\n\ufeffmedian_r2\tfull_flanks\r\n" +
                                       "".join(f"{value}\tTrue\r\n" for value in backgrounds))
                Path(str(prefix) + ".backgrounds.tsv").write_text(background_text)
            output = root / "result"
            result = subprocess.run(["bash", str(script), str(root), str(output)],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            with Path(str(output) + ".per_region.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            with Path(str(output) + ".overall.tsv").open() as handle:
                overall = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 2)
            self.assertAlmostEqual(float(rows[0]["A_i"]), 0.5)
            self.assertAlmostEqual(float(rows[1]["A_i"]), 1.0)
            self.assertAlmostEqual(float(overall["A"]), 0.75)


if __name__ == "__main__":
    unittest.main()

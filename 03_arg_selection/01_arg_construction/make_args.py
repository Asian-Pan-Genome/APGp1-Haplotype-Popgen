import argparse
import glob
import os
import sys

from util import load_config, run, remove_ext, which
from igd_run_tsinfer import tsinfer_from_igd, load_ratemap

CONFIG = load_config()


def mkoutdir(name: str):
    name = os.path.join(CONFIG.output_dir, name)
    if not os.path.exists(name):
        os.mkdir(name)
    return name


def make_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("prefix", help="input filename prefix")
    parser.add_argument("kind", help="singer or tsinfer")
    parser.add_argument("ref", help="chm13 or grch38")
    parser.add_argument(
        "-j",
        "--jobs",
        default=1,
        type=int,
        help="Number of jobs/threads to use. Defaults to 1.",
    )
    parser.add_argument(
        "-c",
        "--chromosome",
        default=None,
        help="Restrict to a specific chromosome number (1-22).",
    )
    args = parser.parse_args()
    assert args.kind in ("singer", "tsinfer")
    assert args.ref in ("chm13", "grch38")

    if args.chromosome is None:
        chromosomes = list(range(1, 23))
    else:
        chromosomes = [int(args.chromosome)]

    polarized = False
    if "06_polarized" in args.prefix or "07_intersect" in args.prefix:
        print("DATA ALREADY POLARIZED", file=sys.stderr)
        polarized = True

    def get_files_for(chrom_list, suffix):
        candidates = glob.glob(f"{args.prefix}*.{suffix}")
        if len(candidates) < len(chrom_list):
            print(
                f"Not enough {suffix} files for the chromosomes being analyzed. {len(candidates)} {suffix} files, and {len(chrom_list)} chromosomes.",
                file=sys.stderr,
            )
            exit(2)
        chrom2file = {}
        for c in candidates:
            parts = c.split("/")[-1].split(".")
            chrom = None
            for p in parts:
                if p.startswith("chr"):
                    try:
                        chrom = int(p[3:])
                        break
                    except ValueError:
                        pass
            if chrom is not None:
                assert chrom >= 1 and chrom <= 22, f"Invalid chromosome number: {chrom}"
                assert (
                    chrom not in chrom2file
                ), f"Ambiguous file prefix (multiple matches for chromosome {chrom})"
                chrom2file[chrom] = c
        return [chrom2file[chrom] for chrom in chrom_list]

    mkoutdir("ARGS")
    if args.kind == "tsinfer":
        file_suffix = "igd"
        print("Importing tsdate... (slow)", file=sys.stderr)
        import tsdate

        print("... done", file=sys.stderr)
        parallel_singer = None
        outdir = mkoutdir(os.path.join("ARGS", "tsinfer"))
    else:
        file_suffix = "vcf"
        tsdate = None
        parallel_singer = which("parallel_singer", required=True)
        outdir = mkoutdir(os.path.join("ARGS", "singer"))

    for filename, chrom in zip(get_files_for(chromosomes, file_suffix), chromosomes):
        if args.ref == "chm13":
            ratemap_file = os.path.join(
                CONFIG.chm13_ratemaps, f"ratemap.avgmask.chr{chrom}.txt"
            )
        else:
            ratemap_file = os.path.join(
                CONFIG.grch38_ratemaps, f"ratemap_Hg38_chr{chrom}.txt"
            )

        if args.kind == "tsinfer":
            rec_rate = load_ratemap(ratemap_file)

            if polarized:
                ancestral = "REF"
            elif args.ref == "chm13":
                ancestral = os.path.join(
                    CONFIG.chm13_ancestral, f"chr{chrom}.mleRecon.human_anc.fa"
                )
            else:
                ancestral = os.path.join(
                    CONFIG.grch38_ancestral, f"homo_sapiens_ancestor_{chrom}.fa"
                )

            print(f"Planning to run tsinfer on {filename}", file=sys.stderr)

            # Do the inference.
            ts = tsinfer_from_igd(
                filename,
                ancestral,
                [],
                [],
                rec_rate=rec_rate,
                jobs=args.jobs,
            )

            simplified_ts = tsdate.preprocess_ts(ts)
            dated_ts = tsdate.date(simplified_ts, mutation_rate=CONFIG.mut_rate)

            out_filename = os.path.join(
                outdir, f"{os.path.basename(filename)}.tsdate.trees"
            )
            print(f"Writing tree-seq to {out_filename}", file=sys.stderr)
            dated_ts.dump(out_filename)
        else:
            assert polarized, "SINGER can only be run on already-polarized data"
            jobs = 2 if args.jobs < 2 else args.jobs
            vcf_prefix = remove_ext(filename, ext="vcf")
            out_prefix = os.path.join(outdir, os.path.basename(vcf_prefix))
            cmd = [
                parallel_singer,
                "-Ne",
                str(CONFIG.ne),
                "-m",
                str(CONFIG.mut_rate),
                "-recomb_map",
                ratemap_file,
                "-mut_map",
                CONFIG.fake_mutmap,
                "-n",
                str(CONFIG.mcmc_samples),
                "-thin",
                str(CONFIG.mcmc_thin),
                "-num_cores",
                str(jobs),
                "-vcf",
                vcf_prefix,
                "-output",
                f"{out_prefix}.singer",
            ]
            print(f"Running command: {cmd}", file=sys.stderr)
            run(cmd)

if __name__ == "__main__":
    make_args()
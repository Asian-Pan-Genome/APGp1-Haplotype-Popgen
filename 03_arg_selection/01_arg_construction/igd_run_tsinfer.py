# Run tsinfer on a single IGD file
#
# 1. Any multi-allelic sites will be mapped but will not contribute to the tree
#    generation. These sites _MAY_ impact the later use of tsdate, which cares
#    about segregating sites (or mutations on branches).
# 2. Any missing alleles will be imputed by tsinfer. It is not clear how this affects
#    the resulting tsdate behavior, but it should not affect the tree topology.
# 3. Non-SNVs will cause your IGD file to be rejected, you must filter them out (e.g.
#    with `igdtools --drop-non-snvs`)
# 4. There are two options for polarization:
#    (a) Polarize the data such that the REF is the ancestral allele for all variants.
#    (b) Pass in another IGD file (which does not need to have any samples, just
#    variants) which has the REF as the ancestral allele.
#
# The population labels input file should be formatted as:
#  <INDIV ID> <POP NAME>
# where <INDIV ID> matches an individual ID from the IGD file. All individuals should
# be listed.
#
# The rate map file is the same as what SINGER takes as input for recombination/mutation
# rates, and is equivalent to an msprime.RateMap, with each line being:
#  <start position> <end position> <rate>
# The map will be extended to position 0.0 and position +500Mbp using the first and last
# rate entry.
#
# If you have a VCF(.GZ) file, use `igdtools` to generate an IGD file.
from tskit import RateMap
from typing import List, Tuple, Union
import argparse
import os
import pyfaidx
import pyigd
import sys
import tsinfer
import tskit

MIN_REC_RATE = 1e-12


def tsinfer_from_igd(
    igd_file: str,
    ancestral: str,
    popnames: List[str] = [],
    popmap: List[Tuple[str, int]] = [],
    rec_rate: Union[RateMap, float] = 1e-8,
    jobs: int = 1,
    skip_noaa: bool = False,
    skip_multi: bool = False,
    skip_missing: bool = False,
    verbose: bool = True,
):
    with open(igd_file, "rb") as f_igd:
        igd_file = pyigd.IGDReader(f_igd)
        # The length of the sequence is the last variant's position + 1
        position, flags = igd_file.get_position_and_flags(igd_file.num_variants - 1)
        seq_len = position + 1

        if ancestral == "REF":

            def get_ancestral(position, allele_list):
                return 0

        else:
            assert os.path.isfile(ancestral)
            ancestral_map = {}
            if ancestral.endswith(".igd"):
                with open(ancestral, "rb") as af:
                    ancestral_igd = pyigd.IGDReader(af)
                    for i in range(ancestral_igd.num_variants):
                        position, flags = ancestral_igd.get_position_and_flags(i)
                        ref = ancestral_igd.get_ref_allele(i)
                        if position not in ancestral_map:
                            ancestral_map[position] = ref
                        else:
                            assert (
                                ancestral_map[position] == ref
                            ), f"Multiple REF alleles at position {position}: {ref}, {ancestral_map[position]}"
            else:
                fasta_reader = pyfaidx.Fasta(ancestral)
                assert (
                    len(fasta_reader.values()) == 1
                ), "Your FASTA file has more than one sequence (should be only a single ancestral seqeuence)"
                for pos, allele in enumerate("X" + str(list(fasta_reader.values())[0])):
                    allele = allele.upper()
                    if allele in ("A", "C", "T", "G"):
                        ancestral_map[pos] = allele

            def get_ancestral(position, allele_list):
                ancestral_allele = ancestral_map.get(position)
                if ancestral_allele is not None:
                    try:
                        return allele_list.index(ancestral_allele)
                    except ValueError:
                        return tskit.MISSING_DATA
                return tskit.MISSING_DATA

        with tsinfer.SampleData(sequence_length=seq_len) as sample_data:
            # This order is gauranteed to match the order of first occurrence in the pop labels
            for name in popnames:
                sample_data.add_population(metadata={"name": name})
            # Define individuals mapped to populations
            for indiv, pop_id in popmap:
                sample_data.add_individual(
                    ploidy=igd_file.ploidy, population=pop_id, metadata={"name": indiv}
                )

            # Add the sites. IGD has only bi-allelic variants, so we construct each site from one or more variants
            sites = 0
            used = 0
            multi_sites = 0
            miss_sites = 0
            anc_sites = 0
            current_position = -1
            current_samples = None
            current_alleles = None
            current_has_missing = False
            for variant_index in range(igd_file.num_variants):
                ref = igd_file.get_ref_allele(variant_index)
                alt = igd_file.get_alt_allele(variant_index)
                position, is_missing, sample_list = igd_file.get_samples(variant_index)
                assert (
                    position <= seq_len
                ), "Variant {variant_index} has position {position} > {seq_len}"
                assert (
                    len(ref) == 1
                ), f"Non-SNV (allele is not a single nucleotide) at position {position}"
                assert (
                    len(alt) <= 1
                ), f"Non-SNV (allele is not a single nucleotide) at position {position}"

                if current_position != position:
                    if current_position >= 0:
                        ancestral_index = get_ancestral(
                            current_position, current_alleles
                        )
                        sites += 1
                        if ancestral_index != tskit.MISSING_DATA:
                            anc_sites += 1
                        add_site = True
                        if skip_noaa and (ancestral_index == tskit.MISSING_DATA):
                            add_site = False
                        if skip_multi and len(current_alleles) > 2:
                            add_site = False
                        if skip_missing and current_has_missing:
                            add_site = False
                        if add_site:
                            used += 1
                            sample_data.add_site(
                                current_position,
                                current_samples,
                                current_alleles,
                                ancestral_allele=ancestral_index,
                            )
                    current_position = position
                    current_samples = [0 for _ in range(igd_file.num_samples)]
                    current_alleles = [ref, alt]
                    current_has_missing = False
                elif not is_missing:
                    assert (
                        current_alleles[0] == ref
                    ), f"Multiple REF alleles at position {position}"
                    if len(current_alleles) == 2:
                        multi_sites += 1
                    current_alleles.append(alt)
                if is_missing:
                    current_has_missing = True
                    miss_sites += 1
                    for sample in sample_list:
                        current_samples[sample] = tskit.MISSING_DATA
                else:
                    allele_value = len(current_alleles) - 1
                    for sample in sample_list:
                        current_samples[sample] = allele_value
            if current_position >= 0:
                ancestral_index = get_ancestral(current_position, current_alleles)
                sites += 1
                if ancestral_index != tskit.MISSING_DATA:
                    anc_sites += 1
                add_site = True
                if skip_noaa and (ancestral_index == tskit.MISSING_DATA):
                    add_site = False
                if skip_multi and len(current_alleles) > 2:
                    add_site = False
                if skip_missing and current_has_missing:
                    add_site = False
                if add_site:
                    used += 1
                    sample_data.add_site(
                        current_position,
                        current_samples,
                        current_alleles,
                        ancestral_allele=ancestral_index,
                    )

            print("Finished constructing SampleData...", file=sys.stderr)
            print(f"  Total sites: {sites}", file=sys.stderr)
            print(f"  Multi-allelic: {multi_sites}", file=sys.stderr)
            print(f"  Sites with missing: {miss_sites}", file=sys.stderr)
            print(f"  Sites with ancestral: {anc_sites}", file=sys.stderr)
            print(f"  Used sites: {used}", file=sys.stderr)
            print("Running tsinfer", file=sys.stderr)
            sys.stderr.flush()

            sample_data.finalise()
            inferred_ts = tsinfer.infer(
                sample_data, recombination_rate=rec_rate, num_threads=jobs
            )
    return inferred_ts


def load_poplabels(poplabel_file: str) -> Tuple[List[str], List[Tuple[str, int]]]:
    popname_map = {}
    popnames = []
    popmap = []
    with open(poplabel_file) as f:
        for line in f:
            line = line.strip()
            if line:
                iid, popname = line.split()
                if popname not in popname_map:
                    popname_map[popname] = len(popnames)
                    popnames.append(popname)
                popmap.append((iid, popname_map[popname]))
    return popnames, popmap


def load_ratemap(ratemap_file: str) -> RateMap:
    positions = []
    rates = []
    with open(ratemap_file) as f:
        for i, line in enumerate(f):
            line = line.strip()
            if line:
                start, end, rate = line.split()
                rate = float(rate)
                if rate < MIN_REC_RATE:
                    rate = MIN_REC_RATE
                if i == 0:
                    positions.append(0.0)
                    rates.append(float(rate))
                positions.append(float(start))
                rates.append(float(rate))
    positions.append(max(float(end), 500_000_000))
    return RateMap(position=positions, rate=rates)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("igd_file", help="The IGD file to use for inference")
    parser.add_argument(
        "ancestral",
        help="Either the string 'REF' or a filename. If an IGD file, then assume the "
        "REF allele is the ancestral allele. If a FASTA file, assume it captures the ancestral allele starting at position 1.",
    )
    parser.add_argument(
        "--pop-labels",
        help="Two-column text file: col1 is individual ID, col2 is popname",
    )
    parser.add_argument(
        "--jobs",
        "-j",
        default=1,
        type=int,
        help="Number of jobs/threads to use for tsinfer",
    )
    parser.add_argument(
        "--recomb-rate",
        "-r",
        default=1e-8,
        help="Recombination rate map or constant rate to use.",
    )
    parser.add_argument(
        "--skip-noaa",
        action="store_true",
        help="Skip sites with unknown ancestral allele.",
    )
    parser.add_argument(
        "--skip-multi", action="store_true", help="Skip sites with more than 2 alleles."
    )
    parser.add_argument(
        "--skip-missing", action="store_true", help="Skip sites with missing data."
    )
    args = parser.parse_args()

    # Population labels
    popnames, popmap = ([], [])
    if args.pop_labels is not None:
        print(f"Loading population labels from {args.pop_labels}", file=sys.stderr)
        popnames, popmap = load_poplabels(args.pop_labels)

    # Recombination rate
    try:
        rec_rate = float(args.recomb_rate)
    except ValueError:
        rec_rate = load_ratemap(args.recomb_rate)

    # Do the inference.
    ts = tsinfer_from_igd(
        args.igd_file,
        args.ancestral,
        popnames,
        popmap,
        rec_rate=rec_rate,
        jobs=args.jobs,
        skip_noaa=args.skip_noaa,
        skip_multi=args.skip_multi,
        skip_missing=args.skip_missing,
    )
    out_filename = f"{os.path.basename(args.igd_file)}.trees"
    print(f"Writing tree-seq to {out_filename}", file=sys.stderr)
    ts.dump(out_filename)


if __name__ == "__main__":
    main()

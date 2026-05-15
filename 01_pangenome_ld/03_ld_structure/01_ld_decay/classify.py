#!/usr/bin/env python3

###############################################################################
# Classify structural variants (SVs) into internal and boundary SVs
# based on haplotype block coordinates and SV membership.
#
# Usage:
#   python classify.py \
#       --blocks autosome.block \
#       --svs input.svs.txt \
#       --output-prefix classified
#
# Description:
#   This script classifies structural variants (SVs) into:
#
#       1. Internal SVs
#          - SVs located within haplotype blocks
#          - Not located at block boundaries
#
#       2. Boundary SVs
#          - First or last SVs within a haplotype block
#
# Input:
#   1. Haplotype block file
#   2. SV table file
#
# Output:
#   <prefix>.internal_svs.txt
#   <prefix>.boundary_svs.txt
#
# Requirements:
#   - Python >= 3.8
###############################################################################

import argparse
from pathlib import Path


# ---------------------------- Parse haploblocks ----------------------------- #

def parse_haploblock_file(haploblock_file):
    """
    Parse haplotype block file.

    Expected columns:
        chrom start end ... sv_list

    Example SV list:
        >29395>29398|>29398>29401
    """

    blocks = []

    with open(haploblock_file, "r") as f:

        for line in f:

            fields = line.strip().split()

            if len(fields) < 6:
                continue

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])

            sv_list = fields[5].split("|")

            blocks.append({
                "chrom": chrom,
                "start": start,
                "end": end,
                "sv_ids": sv_list
            })

    return blocks


# -------------------------------- Parse SVs -------------------------------- #

def parse_sv_file(sv_file):
    """
    Parse SV file.

    Expected columns:
        chrom sv_id ... pos ref alt
    """

    svs = []

    with open(sv_file, "r") as f:

        for line in f:

            fields = line.strip().split()

            if len(fields) < 6:
                continue

            svs.append({
                "chrom": fields[0],
                "sv_id": fields[1],
                "pos": int(fields[3]),
                "ref": fields[4],
                "alt": fields[5],
                "original_line": line.strip()
            })

    return svs


# ------------------------------- Classification ----------------------------- #

def classify_svs(blocks, svs):
    """
    Classify SVs into:
        - internal SVs
        - boundary SVs
    """

    internal_svs = []
    boundary_svs = []

    # Record boundary SVs from each block
    boundary_sv_ids = set()

    for block in blocks:

        if len(block["sv_ids"]) == 0:
            continue

        boundary_sv_ids.add(block["sv_ids"][0])
        boundary_sv_ids.add(block["sv_ids"][-1])

    # Classify SVs
    for sv in svs:

        classified = False

        for block in blocks:

            if (
                sv["chrom"] == block["chrom"]
                and block["start"] <= sv["pos"] <= block["end"]
            ):

                if sv["sv_id"] in boundary_sv_ids:

                    boundary_svs.append(sv)

                elif sv["sv_id"] in block["sv_ids"]:

                    internal_svs.append(sv)

                classified = True
                break

        # Boundary SV outside coordinate range
        if (
            not classified
            and sv["sv_id"] in boundary_sv_ids
        ):
            boundary_svs.append(sv)

    return internal_svs, boundary_svs


# ------------------------------- Write output ------------------------------- #

def write_output(output_file, sv_list):

    with open(output_file, "w") as f:

        for sv in sv_list:
            f.write(sv["original_line"] + "\n")


# ----------------------------------- Main ---------------------------------- #

def main():

    parser = argparse.ArgumentParser(
        description="Classify SVs into internal and boundary SVs"
    )

    parser.add_argument(
        "-b", "--blocks",
        required=True,
        help="Input haplotype block file"
    )

    parser.add_argument(
        "-s", "--svs",
        required=True,
        help="Input SV file"
    )

    parser.add_argument(
        "-o", "--output-prefix",
        default="classified",
        help="Output file prefix (default: classified)"
    )

    args = parser.parse_args()

    # ----------------------------- Input check ------------------------------ #

    if not Path(args.blocks).exists():
        raise FileNotFoundError(
            f"Haplotype block file not found: {args.blocks}"
        )

    if not Path(args.svs).exists():
        raise FileNotFoundError(
            f"SV file not found: {args.svs}"
        )

    # -------------------------------- Logging ------------------------------ #

    print("============================================================")
    print("SV classification based on haplotype blocks")
    print(f"Block file : {args.blocks}")
    print(f"SV file    : {args.svs}")
    print("============================================================")

    # ------------------------------- Load data ------------------------------ #

    print("\n[Step 1] Parsing haplotype blocks...")
    blocks = parse_haploblock_file(args.blocks)

    print(f"Detected {len(blocks)} haplotype blocks")

    print("\n[Step 2] Parsing SV file...")
    svs = parse_sv_file(args.svs)

    print(f"Detected {len(svs)} SVs")

    # ----------------------------- Classification -------------------------- #

    print("\n[Step 3] Classifying SVs...")

    internal_svs, boundary_svs = classify_svs(blocks, svs)

    # ------------------------------- Output -------------------------------- #

    internal_out = f"{args.output_prefix}.internal_svs.txt"
    boundary_out = f"{args.output_prefix}.boundary_svs.txt"

    write_output(internal_out, internal_svs)
    write_output(boundary_out, boundary_svs)

    # ------------------------------ Statistics ----------------------------- #

    classified_n = len(internal_svs) + len(boundary_svs)
    unclassified_n = len(svs) - classified_n

    print("\n============================================================")
    print("Classification completed successfully")
    print("============================================================")

    print(f"Total SVs        : {len(svs)}")
    print(f"Internal SVs     : {len(internal_svs)}")
    print(f"Boundary SVs     : {len(boundary_svs)}")
    print(f"Unclassified SVs : {unclassified_n}")

    print("\nOutput files:")
    print(f"  - {internal_out}")
    print(f"  - {boundary_out}")

    print("============================================================")


if __name__ == "__main__":
    main()

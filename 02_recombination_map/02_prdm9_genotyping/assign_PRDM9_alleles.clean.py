#!/usr/bin/env python3

import argparse
import re
import sys
from collections import OrderedDict

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


ZF_DNA_LEN = 84
ZF_AA_LEN = 28

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def norm_dna(seq):
    return re.sub(r"\s+", "", str(seq)).upper().replace("U", "T")


def revcomp(seq):
    return norm_dna(seq).translate(_COMP)[::-1]


def chunks(seq, size):
    for i in range(0, len(seq), size):
        yield seq[i:i + size]


def split_zf_code(zf_code):
    zf_code = str(zf_code).strip()

    tokens = re.findall(r"[:|!][^:|!]+", zf_code)

    if "".join(tokens) != zf_code:
        raise ValueError(
            f"Cannot parse zf_code cleanly: {zf_code} -> {tokens}"
        )

    return tokens


def translate_dna(seq):
    if len(seq) % 3 != 0:
        return "NA"
    return str(Seq(seq).translate())


def get_dna_contact_aas(aa_seq):
    if aa_seq == "NA":
        return "NA"

    contacts = []
    for i in range(0, len(aa_seq), ZF_AA_LEN):
        block = aa_seq[i:i + ZF_AA_LEN]
        if len(block) < 13:
            return "NA"
        contacts.append(block[6] + block[7] + block[9] + block[12])

    return "".join(contacts)


def safe_get(row, col, default=""):
    if row is None or col not in row.index:
        return default
    val = row[col]
    if pd.isna(val) or str(val) == "":
        return default
    return val


def build_zf_reference(lib):
    unit_to_code = {}

    for idx, row in lib.iterrows():
        dna = row["dna_sequence"]
        zf_code = row["zf_code"]

        if len(dna) % ZF_DNA_LEN != 0:
            sys.exit(
                f"ERROR: reference allele dna_sequence length is not a multiple of 84 bp "
                f"at row {idx + 2}: {len(dna)} bp"
            )

        n_units = len(dna) // ZF_DNA_LEN

        try:
            zf_tokens = split_zf_code(zf_code)
        except ValueError as e:
            sys.exit(f"ERROR at reference row {idx + 2}: {e}")
        for token, unit_seq in zip(zf_tokens, chunks(dna, ZF_DNA_LEN)):
            if unit_seq in unit_to_code and unit_to_code[unit_seq] != token:
                sys.exit(
                    f"ERROR: the same 84-bp ZF unit has conflicting codes: "
                    f"{unit_to_code[unit_seq]} vs {token}\n{unit_seq}"
                )

            unit_to_code[unit_seq] = token

    return unit_to_code


def decompose_known(seq, unit_to_code):
    items = []

    for unit_seq in chunks(seq, ZF_DNA_LEN):
        rc_seq = revcomp(unit_seq)

        if unit_seq in unit_to_code:
            items.append({
                "token": unit_to_code[unit_seq],
                "sign": "+",
                "known": True,
                "canonical_seq": unit_seq
            })
        elif rc_seq in unit_to_code:
            items.append({
                "token": unit_to_code[rc_seq],
                "sign": "-",
                "known": True,
                "canonical_seq": rc_seq
            })
        else:
            items.append({
                "token": None,
                "sign": "+",
                "known": False,
                "canonical_seq": unit_seq
            })

    return items


def score_orientation(seq, lib, unit_to_code):
    items = decompose_known(seq, unit_to_code)

    known_count = sum(x["known"] for x in items)
    plus_count = sum(x["known"] and x["sign"] == "+" for x in items)
    minus_count = sum(x["known"] and x["sign"] == "-" for x in items)

    has_unknown = any(not x["known"] for x in items)
    zf_code = "".join(x["token"] if x["known"] else "?" for x in items)

    exact_code_seq = False
    exact_seq = (lib["dna_sequence"] == seq).any()
    exact_code = False

    if not has_unknown:
        exact_code_seq = (
            (lib["zf_code"] == zf_code) &
            (lib["dna_sequence"] == seq)
        ).any()
        exact_code = (lib["zf_code"] == zf_code).any()

    return (
        int(exact_code_seq),
        int(exact_seq),
        int(exact_code),
        known_count,
        plus_count,
        -minus_count
    )


def choose_orientation(raw_seq, lib, unit_to_code):
    forward_seq = raw_seq
    reverse_seq = revcomp(raw_seq)

    forward_score = score_orientation(forward_seq, lib, unit_to_code)
    reverse_score = score_orientation(reverse_seq, lib, unit_to_code)

    if reverse_score > forward_score:
        return reverse_seq, "reverse_complement", reverse_score

    return forward_seq, "forward", forward_score


def assign_zf_units(seq, unit_to_code, new_seq_to_name, new_seq_order):
    clean_units = []
    signed_units = []
    n_new_units_in_sample = 0

    for unit_seq in chunks(seq, ZF_DNA_LEN):
        rc_seq = revcomp(unit_seq)

        if unit_seq in unit_to_code:
            token = unit_to_code[unit_seq]
            sign = "+"
        elif rc_seq in unit_to_code:
            token = unit_to_code[rc_seq]
            sign = "-"
        else:
            n_new_units_in_sample += 1

            if unit_seq in new_seq_to_name:
                token = new_seq_to_name[unit_seq]
                sign = "+"
            elif rc_seq in new_seq_to_name:
                token = new_seq_to_name[rc_seq]
                sign = "-"
            else:
                token = f"ZF_New_{len(new_seq_order) + 1}"
                new_seq_to_name[unit_seq] = token
                new_seq_order.append(unit_seq)
                sign = "+"

        clean_units.append(token)
        signed_units.append(f"{token}{sign}")

    zf_code = "".join(clean_units)

    return zf_code, "|".join(clean_units), "|".join(signed_units), n_new_units_in_sample


def match_reference_allele(lib, zf_code, dna_seq):
    exact = lib[
        (lib["zf_code"] == zf_code) &
        (lib["dna_sequence"] == dna_seq)
    ]

    if len(exact) > 0:
        if len(exact) > 1:
            print(
                f"WARNING: multiple reference rows match zf_code={zf_code}; using the first one.",
                file=sys.stderr
            )
        return "published", exact.iloc[0]

    code_only = (lib["zf_code"] == zf_code).any()
    seq_only = (lib["dna_sequence"] == dna_seq).any()

    if code_only and seq_only:
        return "code_and_sequence_present_but_not_same_row", None
    if code_only:
        return "zf_code_only", None
    if seq_only:
        return "sequence_only", None

    return "novel", None


def write_new_zfs(path, new_seq_to_name, new_seq_order):
    with open(path, "w") as fout:
        for seq in new_seq_order:
            name = new_seq_to_name[seq]
            fout.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fout.write(seq[i:i + 80] + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Assign PRDM9 alleles from PRDM9 ZF-array FASTA using PRDM9.alleles.lib.tsv."
    )
    parser.add_argument("samples_fa", help="Input FASTA containing sample PRDM9 ZF-array sequences.")
    parser.add_argument("allele_lib", help="Reference PRDM9 allele library TSV. Required columns: zf_code, dna_sequence.")
    parser.add_argument(
        "--allele-out",
        default="samples.PRDM9_allele.list",
        help="Output allele assignment table. Default: samples.PRDM9_allele.list"
    )
    parser.add_argument(
        "--new-zf-out",
        default="ZF.new.fa",
        help="Output FASTA of newly observed 84-bp ZF units. Default: ZF.new.fa"
    )
    parser.add_argument(
        "--strict-conflict",
        action="store_true",
        help="Exit if zf_code-only or sequence-only conflicts are observed."
    )

    args = parser.parse_args()

    lib = pd.read_csv(args.allele_lib, sep="\t", dtype=str)

    required_cols = {"zf_code", "dna_sequence"}
    missing = required_cols - set(lib.columns)
    if missing:
        sys.exit(f"ERROR: missing required columns in allele library: {', '.join(sorted(missing))}")

    lib["zf_code"] = lib["zf_code"].astype(str).str.strip()
    lib["dna_sequence"] = lib["dna_sequence"].map(norm_dna)

    if "dna_contact_aas" in lib.columns:
        lib["dna_contact_aas"] = lib["dna_contact_aas"].astype(str).str.upper()

    unit_to_code = build_zf_reference(lib)

    records = []
    unknown_allele_to_id = OrderedDict()
    new_seq_to_name = OrderedDict()
    new_seq_order = []

    for rec in SeqIO.parse(args.samples_fa, "fasta"):
        sample = rec.id
        raw_seq = norm_dna(str(rec.seq))

        if len(raw_seq) % ZF_DNA_LEN != 0:
            sys.exit(
                f"ERROR: sample sequence length is not a multiple of 84 bp: "
                f"{sample}, length={len(raw_seq)}"
            )

        dna_seq, orientation, orientation_score = choose_orientation(
            raw_seq, lib, unit_to_code
        )

        zf_code, zf_units, zf_code_with_strand, n_new_units = assign_zf_units(
            dna_seq, unit_to_code, new_seq_to_name, new_seq_order
        )

        match_status, ref_row = match_reference_allele(lib, zf_code, dna_seq)

        if args.strict_conflict and match_status not in {"published", "novel"}:
            sys.exit(
                f"ERROR: reference conflict for {sample}: "
                f"zf_code={zf_code}, match_status={match_status}"
            )

        if match_status == "published":
            short_id = safe_get(ref_row, "short_ID", safe_get(ref_row, "ID", zf_code))
            published_allele = True
            in_pop = safe_get(ref_row, "in_pop", "")
            aa_sequence = safe_get(ref_row, "aa_sequence", translate_dna(dna_seq))
            dna_contact_aas = safe_get(ref_row, "dna_contact_aas", get_dna_contact_aas(aa_sequence))
            actype = safe_get(ref_row, "ACtype", "")
        else:
            key = (zf_code, dna_seq)
            if key not in unknown_allele_to_id:
                unknown_allele_to_id[key] = f"PRDM9_New_{len(unknown_allele_to_id) + 1}"

            short_id = unknown_allele_to_id[key]
            published_allele = False
            in_pop = False
            aa_sequence = translate_dna(dna_seq)
            dna_contact_aas = get_dna_contact_aas(aa_sequence)
            actype = "unknown"

        records.append({
            "sample": sample,
            "short_ID": short_id,
            "published_allele": published_allele,
            "in_pop": in_pop,
            "zf_code": zf_code,
            "zf_units": zf_units,
            "zf_code_with_strand": zf_code_with_strand,
            "orientation": orientation,
            "match_status": match_status,
            "n_zf_units": len(dna_seq) // ZF_DNA_LEN,
            "n_new_zf_units_in_sample": n_new_units,
            "dna_sequence": dna_seq,
            "aa_sequence": aa_sequence,
            "dna_contact_aas": dna_contact_aas,
            "ACtype": actype
        })

    out_df = pd.DataFrame(records)
    out_df.to_csv(args.allele_out, sep="\t", index=False)

    write_new_zfs(args.new_zf_out, new_seq_to_name, new_seq_order)

    print(f"Wrote allele assignments: {args.allele_out}", file=sys.stderr)
    print(f"Wrote new ZF units: {args.new_zf_out}", file=sys.stderr)
    print(f"New ZF unit count: {len(new_seq_order)}", file=sys.stderr)
    print(f"New PRDM9 allele count: {len(unknown_allele_to_id)}", file=sys.stderr)


if __name__ == "__main__":
    main()

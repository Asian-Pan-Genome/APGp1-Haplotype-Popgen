import sys
import argparse
import gzip
from typing import Tuple, List, Optional
from pathlib import Path

parser = argparse.ArgumentParser(
        description="Generate reference and alternative FASTA sequences from BED file with variants",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
BED file format (tab-separated):
    chrom    start    end    ref    alt    [optional_columns...]
    
Example:
    chr1    100    101    A    G
    chr2    500    502    AT    G
    
Note: BED format uses 0-based coordinates (start is 0-based, end is 1-based).
        """
    )
    
parser.add_argument("-b", "--bed", required=True, help="Input BED file with variants")
parser.add_argument("-fa", "--genome", required=True, help="Reference genome FASTA file")
parser.add_argument("-ref", default="ref_ale.fa", help="Output FASTA file for real sequences (with REF)")
parser.add_argument("-alt", default="alt_ale.fa", help="Output FASTA file for fake sequences (with ALT)")
parser.add_argument("-flk", "--flanking", type=int, default=5000, 
                       help="Number of flanking base pairs to include on each side (default: 5000)")
    
args = parser.parse_args()


def parse_bed_file(bed_path: str) -> List[Tuple[str, int, int, str, str,str]]:
    """
    Parse BED file with variant information.
    
    Args:
        bed_path: Path to BED file
        
    Returns:
        List of tuples (chrom, start, end, ref, alt)
    """
    variants = []
    
    # Check if file is gzipped
    open_func = gzip.open if bed_path.endswith('.gz') else open
    mode = 'rt' if bed_path.endswith('.gz') else 'r'
    
    with open_func(bed_path, mode) as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
                
            fields = line.split('\t')
            if len(fields) < 5:
                print(f"Warning: Line {line_num} has fewer than 5 fields, skipping: {line}", 
                      file=sys.stderr)
                continue
            
            try:
                chrom = fields[0]
                start = int(fields[1])  # 0-based start in BED format
                end = int(fields[2])    # 0-based end (exclusive)
                faid = fields[3]
                ref = fields[4].upper()
                alt = fields[5].upper()
                
                variants.append((chrom, start, end, faid, ref, alt))
            except ValueError as e:
                print(f"Warning: Line {line_num} has invalid format, skipping: {line}", 
                      file=sys.stderr)
                continue
    
    return variants

def readFasta(genome_file):
    seq = []
    info = ""
    open_func = gzip.open if genome_file.endswith('.gz') else open
    mode = 'rt' if genome_file.endswith('.gz') else 'r'
    
    with open_func(genome_file, mode) as fh:
        for line in fh:
            line = line.rstrip()
            line = line.replace("\n","")
            if line.startswith('>'):
                if not seq == []:
                    yield seq, info[1:]
                info = line.replace("\n","")
                seq.clear()
                continue
            else:
                seq.append(line.replace("\n",""))
        if info != "":
            yield seq, info[1:]

def generate_variant_sequences(variants, genome_file, flanking_bp = 5000) -> Tuple[List[str], List[str]]:
    """
    Generate real (with ALT) and fake (with REF) sequences.
    
    Args:
        variants: List of variants from BED file
        genome_file: Path to reference genome FASTA
        flanking_bp: Number of flanking base pairs to include on each side
        
    Returns:
        Tuple of (real_sequences, fake_sequences)
    """
    real_sequences = []
    fake_sequences = []
    
    i = 0
    for s, c in readFasta(genome_file):
        while i < len(variants):
            (chrom, start, end, faid, ref, alt) = variants[i]
            if c != chrom:
                break
            seq = "".join(s).upper()
            # Calculate positions (BED is 0-based)
            variant_length = end - start
            
            # Validate REF matches the reference genome
            ref_in_genome = seq[start:end]
            if ref_in_genome != ref:
                print(f"Warning: Variant {i}: REF '{ref}' doesn't match genome '{ref_in_genome}' "
                    f"at {chrom}:{start}-{end}. Using genome reference.", 
                    file=sys.stderr)
                ref = ref_in_genome
            
            # Extract flanking regions
            left_flank_start = max(0, start - flanking_bp)
            right_flank_end = min(len(seq), end + flanking_bp)
            
            left_flank = seq[left_flank_start:start]
            right_flank = seq[end:right_flank_end]
            
            # Generate sequences
            real_seq = left_flank + ref + right_flank
            fake_seq = left_flank + alt + right_flank
            
            # Create headers
            real_header = f">{faid} ref:{left_flank_start}-{right_flank_end}"
            fake_header = f">{faid} alt:{left_flank_start}-{right_flank_end}"
            
            real_sequences.append(f"{real_header}\n{real_seq}")
            fake_sequences.append(f"{fake_header}\n{fake_seq}")
            
            print(f"Processed variant {i}: {chrom}:{start+1}-{end} {ref}->{alt}", 
                file=sys.stderr)
            i = i + 1
    return real_sequences, fake_sequences


def write_fasta(sequences: List[str], output_file: str):
    """
    Write sequences to FASTA file.
    
    Args:
        sequences: List of FASTA sequences with headers
        output_file: Output file path
    """
    with open(output_file, 'w') as f:
        for seq in sequences:
            f.write(seq + '\n')


def main():
    # Check input files exist
    for file_path in [args.bed, args.genome]:
        if not Path(file_path).exists():
            print(f"Error: File not found: {file_path}", file=sys.stderr)
            sys.exit(1)
    
    print(f"Reading variants from {args.bed}...", file=sys.stderr)
    variants = parse_bed_file(args.bed)
    
    if not variants:
        print("Error: No valid variants found in BED file", file=sys.stderr)
        sys.exit(1)
    
    print(f"Found {len(variants)} variants", file=sys.stderr)
    print(f"Generating sequences with {args.flanking}bp flanking regions...", file=sys.stderr)
    
    real_seqs, fake_seqs = generate_variant_sequences(variants, args.genome, args.flanking)
    
    print(f"Writing real sequences to {args.ref}...", file=sys.stderr)
    write_fasta(real_seqs, args.ref)
    
    print(f"Writing fake sequences to {args.alt}...", file=sys.stderr)
    write_fasta(fake_seqs, args.alt)
    
    print(f"Done! Generated {len(real_seqs)} sequences in each file.", file=sys.stderr)


if __name__ == "__main__":
    main()

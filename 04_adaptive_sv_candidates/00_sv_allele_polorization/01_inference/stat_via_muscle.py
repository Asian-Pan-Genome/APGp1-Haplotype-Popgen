import sys
import argparse
from typing import Tuple, List, Optional
import gzip

parser = argparse.ArgumentParser(
        description="get the flanking and inner identity of a variant from muscle alignment",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python extract_alignment.py alignment.needle 10 50
        """
    )
    
parser.add_argument("-n", "--needle", required=True, help="Needle format file")
parser.add_argument("-idx", required=True, help="variant index")
parser.add_argument("-v", required=True, help="variant position in reference genome")
parser.add_argument("-f", required=True, help="flank region position in reference genome")
parser.add_argument('-ref', metavar='<STR>', default="REF", help='add a specific prefix to ID')
parser.add_argument('-alt', metavar='<STR>', default="ALT", help='add a specific prefix to ID')

if len(sys.argv) == 1:
    parser.print_help()
    parser.exit()    
args = parser.parse_args()

v_st = int(args.v.split("-")[0])
v_ed = int(args.v.split("-")[1])
f_st = int(args.f.split("-")[0])
f_ed = int(args.f.split("-")[1])
v_real_st = v_st - f_st
v_real_ed = v_ed - f_st

def sv_pos(ref_aligned, start, end):
    ref_pos = 0
    pos_map = {}
    
    for i, char in enumerate(ref_aligned):
        if char != '-':
            pos_map[ref_pos] = i
            ref_pos += 1
    
    # Get alignment columns
    #print(pos_map)
    #print(max(pos_map.keys()), min(pos_map.keys()))
    start_col = pos_map.get(start)
    end_col = pos_map.get(end)
    return start_col, end_col
    

def extract_alignment_interval(ref_aligned, query_aligned, start_col, end_col):
    """
    Extract alignment interval based on reference positions.
    
    Args:
        ref_aligned: Reference sequence with gaps
        query_aligned: Query sequence with gaps
        start: Start position (1-based, without gaps)
        end: End position (1-based, without gaps)
    """
    # Map reference positions to alignment columns
    
    if start_col is None or end_col is None:
        raise ValueError(f"Positions {start}-{end} not found in alignment")
    
    # Extract interval
    ref_interval = ref_aligned[start_col:end_col+1]
    query_interval = query_aligned[start_col:end_col+1]
    align_str = []
    align_base = 0
    gap_base = 0
    total_base = 0
    for r, q in zip(ref_interval, query_interval):
        if r == q and r != "-":
            align_str.append('|')
            align_base += 1
        elif r == '-' and q == '-':
            align_str.append('*')
            continue
        elif r == '-' or q == '-':
            align_str.append(' ')
            gap_base += 1
        else:
            align_str.append('.')
        total_base += 1
    
    return ref_interval, query_interval, ''.join(align_str), align_base, gap_base, total_base

def parse_msa(genome_file: str):
    seq = []
    info = ""
    open_func = gzip.open if genome_file.endswith('.gz') else open
    mode = 'rt' if genome_file.endswith('.gz') else 'r'
    outseq = {}

    with open_func(genome_file, mode) as fh:
        for line in fh:
            line = line.rstrip()
            line = line.replace("\n","")
            if line.startswith('>'):
                if not seq == []:
                    outseq[info] = ''.join(seq)
                info = line.replace("\n","").split()[0].replace(">", "")
                seq.clear()
                continue
            else:
                seq.append(line.replace("\n",""))
        if info != "":
            outseq[info] = ''.join(seq)
    return outseq




def main_simple():
    """Simple command-line interface."""
    
    
    needle_file = args.needle
    outnumbers = {}
    try:
        seqs = parse_msa(needle_file)
        #print(seqs.keys())
        for info in seqs.keys():
            if info == args.ref:
                ref_aligned = seqs[info]
            elif info == args.alt:
                alt_aligned = seqs[info]
        
        for info in seqs.keys():
            if info == args.ref:
                continue
            query_aligned = seqs[info]
            if info == args.alt:
                st, ed = sv_pos(ref_aligned, v_real_st-5000, v_real_st)
                ref_interval, query_interval, align_str, align_base, gap_base, total_base = extract_alignment_interval(
                    ref_aligned, query_aligned, st, ed
                )
                #print(f"Reference: {ref_interval}")
                #print(f"           {align_str}")
                #print(f"Query:     {query_interval}")
                alt_l = align_base/ total_base * 100
                st, ed = sv_pos(ref_aligned, v_real_ed, v_real_ed + 4999)
                ref_interval, query_interval, align_str, align_base, gap_base, total_base = extract_alignment_interval(
                    ref_aligned, query_aligned, st, ed
                )                   
                #print(f"Reference: {ref_interval}")
                #print(f"           {align_str}")
                #print(f"Query:     {query_interval}")
                
                alt_r = align_base/ total_base * 100
                #ref_interval, query_interval, align_str, align_base, gap_base, total_base = extract_alignment_interval(
                #    ref_aligned, query_aligned, v_real_st - 19, v_real_ed + 19
                #) 
                # Print alignment
                #print(f"Reference: {ref_interval}")
                #print(f"           {align_str}")
                #print(f"Query:     {query_interval}")
                continue
            else:
                st, ed = sv_pos(ref_aligned, v_real_st - 9, v_real_ed + 9)
                ref_interval, query_interval, align_str, align_base, gap_base, total_base = extract_alignment_interval(
                    ref_aligned, query_aligned, st, ed
                )
                if total_base <= 4:
                    variant_id_ref = 0
                elif align_base <= 4:
                    variant_id_ref = 0
                else:
                    variant_id_ref = (align_base - 4) / (total_base - 4) * 100
                print(f"REF:    {ref_interval}")
                print(f"        {align_str}")
                print(f"{info}: {query_interval}")
                st, ed = sv_pos(ref_aligned, v_real_st - 9, v_real_ed + 9)
                ref_interval, query_interval, align_str, align_base, gap_base, total_base = extract_alignment_interval(
                    alt_aligned, query_aligned, st, ed
                )
                if total_base - 4 <= 0:
                    variant_id_alt = 0
                elif align_base - 4 <= 0:
                    variant_id_alt = 0
                else:
                    variant_id_alt = (align_base - 4) / (total_base - 4) * 100
                print(f"ALT:    {ref_interval}")
                print(f"        {align_str}")
                print(f"{info}: {query_interval}")
                st, ed = sv_pos(ref_aligned, v_real_st-500, v_real_st)
                ref_interval, query_interval, align_str, align_base, gap_base, total_base = extract_alignment_interval(
                    ref_aligned, query_aligned, st, ed
                )
                    
                l_500 = align_base/ total_base * 100
                st, ed = sv_pos(ref_aligned, v_real_ed, v_real_ed + 50)
                ref_interval, query_interval, align_str, align_base, gap_base, total_base = extract_alignment_interval(
                    ref_aligned, query_aligned, st, ed
                )
                r_500 = align_base/ total_base * 100
                outnumbers[info] = (variant_id_ref, variant_id_alt,l_500, r_500)
        outstr = []
        if outnumbers == {}:
            print("%s\t%.2f\t%.2f\tNA" % (args.idx, alt_l, alt_r))
        for info in outnumbers.keys():
            variant_id_ref, variant_id_alt, l_500, r_500 = outnumbers[info]
            outstr.append("%s;%.2f;%.2f;%.2f;%.2f" % (info, variant_id_ref, variant_id_alt, l_500, r_500))
        print("%s\t%.2f\t%.2f\t%s" % (args.idx, alt_l, alt_r, '\t'.join(outstr)))
        
        #print(f"Reference: {ref_interval}")
        #print(f"           {align_str}")
        #print(f"Query:     {query_interval}")
        #print("variant identity: %.2f" % (align_base/ total_base * 100))
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main_simple()

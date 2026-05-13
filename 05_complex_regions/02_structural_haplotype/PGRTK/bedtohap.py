import argparse
import sys
import textwrap
import math
parser = argparse.ArgumentParser(prog='bedtohap',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
Get the structural haplotype from PGR-TK result
'''))


parser.add_argument('-bed', metavar='file', help='sample list')
parser.add_argument('-chr', metavar='str', help='chromosome')
parser.add_argument('-st', metavar='str', default="0.8", help='')
parser.add_argument('-ed', metavar='str', default="0.6", help='merge total ratio')
if len(sys.argv) == 1:
    parser.print_help()
    parser.exit()

args = parser.parse_args()

def readbed(fh):
    sample_info = {}
    for line in fh:
        line = line.rstrip()
        line = line.replace("\n", "")
        lining = line.split("\t")
        if line.startswith("#"):
            continue
        sam = lining[0]
        cp = lining[3]
        cping = cp.split(":")
        if ":" in cp:
            comp = "%s:%s" % (cping[0], cping[2])
        else:
            comp = cping[0]
        if sam in sample_info.keys():
            sample_info[sam].append(comp)
        else:
            sample_info[sam] = [comp]
    return sample_info

with open(args.bed, "r") as fh:
    sample_info = readbed(fh)
hap_info = {}
for sam in sample_info.keys():
    hap_type = sample_info[sam]
    hap_str = "|".join(hap_type)
    if hap_str in hap_info:
        hap_info[hap_str].append(sam)
    else:
        hap_info[hap_str] = [sam]
for haps in hap_info.keys():
    sample_list = hap_info[haps]
    print("%s\t%s\t%s\t%s\t%i\t%s" % (args.chr, args.st, args.ed, haps, len(sample_list), "|".join(sample_list)))

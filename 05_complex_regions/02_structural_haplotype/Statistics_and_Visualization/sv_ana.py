import sqlite3
import argparse
import sys
import textwrap
from operator import itemgetter
import math

parser = argparse.ArgumentParser(prog='',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
print parameters list 
'''))


parser.add_argument('-db', metavar='file', help='db')
parser.add_argument('-info', metavar='file', help='info')
parser.add_argument('-samorg', metavar='file', help='sample info')
parser.add_argument('-org', metavar='str',default="all", help='origin')
parser.add_argument('-chr', metavar='str',default="all", help='chr')
parser.add_argument('-toM', default=False, action='store_true', help='to matrix file')
parser.add_argument('-ehh', default=False, action='store_true', help='cal ehh')
parser.add_argument('-split', default=False, action='store_true', help='split hap')
#-o <STR>      output file or outdir for -cutf/-cuts
#-z            output a gzip type file
if len(sys.argv) == 1:
    parser.print_help()
    parser.exit()

args = parser.parse_args()
orgs=args.org.split(",")
if args.db:
    con = sqlite3.connect(args.db)
    cur = con.cursor()

    res = cur.execute("SELECT * FROM sample_info;")
    sample_list = res.fetchall()
    total_sample = 0
    sample_num = 0
    #print(sample_list)
    sample_info = {}
    sample_list_f = []
    for info in sample_list:
        sample_info[info[0]] = [info[1],info[2]]
        total_sample = total_sample + 1
    if args.org != "all" and info[1] in orgs:
        sample_num += 1
        sample_list_f.append(info[0])
    if args.org == "all":
        sample_num = total_sample
    con.close()

def read_org(fh):
    sample_info = {}
    for line in fh:
        line = line.rstrip()
        line = line.replace("\n", "")
        lining = line.split("\t")
        sample_info[lining[0]] = [lining[1],lining[2]]
    return sample_info
if not args.db:
    with open(args.samorg, 'r') as fh:
        sample_info = read_org(fh)
   
        
        

def read_info(fh):
    hap_corr = {}
    for line in fh:
        line = line.rstrip()
        line = line.replace("\n", "")
        lining = line.split("\t")
        chr = lining[0]
        block = "%s:%s-%s" % (chr , lining[1], lining[2])
        hap = lining[3]
        samp = []
        for sample in lining[5].split("|"):
            if sample != "ALL":
                for sam in sample_info:
                    if (sample in sam and sample_info[sam][0] in orgs) or args.org == "all":
                        samp.append(sample)
        if lining[5] == "ALL":
            samp = sample_list_f
        al_c = len(samp)
        if block in hap_corr.keys():
            hap_corr[block].append([hap,al_c, "|".join(samp), lining[4], lining[7]])
        else:
            hap_corr[block] = [[hap,al_c, "|".join(samp), lining[4], lining[7]]]
        #gene = "%s%s" % (lining[5], lining[7])
    for block in hap_corr.keys():
        hap_corr[block].sort(reverse=True, key=itemgetter(1))
    return hap_corr
        ##NO direction

def read_rawinfo(fh):
    hap_corr = {}
    for line in fh:
        line = line.rstrip()
        line = line.replace("\n", "")
        lining = line.split("\t")
        chr = lining[0]
        block = "%s:%s-%s" % (chr , lining[1], lining[2])
        hap = lining[3]
        samp = []
        for sample in lining[6].split("|"):
            if sample != "ALL": #and (sample_info[sample][0] == args.org or args.org == "all"):
                for sam in sample_info:
                    if (sam in sample and sample_info[sam][0] in orgs) or args.org == "all":
                        sample_index = sample.replace(sam, "")
                        if (not sample_index[0].isdigit()) and ("#" not in sample_index[1:3]):
                            samp.append(sample)
        if lining[6] == "ALL":
            samp = sample_list_f
        al_c = len(samp)
        if block in hap_corr.keys():
            hap_corr[block].append([hap,al_c, "|".join(samp), lining[4], lining[7]])
        else:
            hap_corr[block] = [[hap, al_c, "|".join(samp), lining[4], lining[7]]]
    for block in hap_corr.keys():
        hap_corr[block].sort(reverse=True, key=itemgetter(1))
    return hap_corr
        #gene = "%s%s" % (lining[5], lining[7])

def intersect(info1, info2):
    
    samp1 = set(info1[2].split("|"))
    samp2 = set(info2[2].split("|"))
    inter = samp1 & samp2
    if inter == set():
        return None
    else:
        return inter

def hap_join(block1, block2):
    i = 0
    join_set = {}
    for info1 in block1:
        for info2 in block2:
            samp_j = intersect(info1, info2)
            if samp_j:
                join_set[i] = samp_j
                i = i + 1
    return join_set

def str2r(block):
    chrom = block.split(":")[0]
    region = block.split(":")[1]
    st = region.split("-")[0]
    ed = region.split("-")[1]
    return chrom,st,ed

def ehh_cal(hap_corr):
    record = []
    for b1 in hap_corr.keys():
        block1 = hap_corr[b1]
        for b2 in hap_corr.keys():
            block2 = hap_corr[b2]
            if b1 != b2 and b2 in record:
                continue
            join_set = hap_join(block1, block2)
            total = 0
            hit_comb = 0
            for j in join_set:
                joint = join_set[j]
                hit_comb = hit_comb + math.comb(len(joint), 2)
                total = total + len(joint)
            ehh = hit_comb / math.comb(total, 2)
            chrom1,st1,ed1 = str2r(b1)
            chrom2,st2,ed2 = str2r(b2)
            print("%s\t%s\t%s\t%.6f\t%s\t%s\t%s" % (chrom1,st1,ed1,ehh, chrom2,st2,ed2))
        record.append(b1)
    return None


def to_matrix(hap_corr, c = 15):
    out_set = {}
    for block in hap_corr.keys():
        haplist = hap_corr[block]
        total = 0
        
        if c <= len(haplist):
            AC_outlist_init = [i[1] for  i in haplist[0:c]]
        else:
            AC_outlist_init = [i[1] for i in haplist] + [0] * (c-len(haplist))
        #print(AC_outlist_init)
        for a in AC_outlist_init:
            total = total + a
        if total == 500:
            r = 500 / total
        else:
            r = 500/total
        #print(AC_outlist_init)
        AC_outlist = [round(i * r, 2) for i in AC_outlist_init]
        AF_outlist = [round(i / sample_num, 4) for i in AC_outlist]
        chrom = block.split(":")[0]
        region = block.split(":")[1]
        st = region.split("-")[0]
        ed = region.split("-")[1]
        middle = (int(st) + int(ed)) / 2
        out_set[block] = [(chrom, middle),AC_outlist, AF_outlist]
    return out_set

with open(args.info, "r") as fh:
    if args.split:
        hap_corr = read_rawinfo(fh)
    else:
        hap_corr = read_info(fh)
    
if args.split:
    
    for block in hap_corr.keys():
        haplist = hap_corr[block]
        for hap in haplist:
            if hap[1] != 0:
                print("%s\t%s\t%s\t%s\t%i\t%s\t%s\t%s" %  (block.split(":")[0], block.split(":")[1].split("-")[0], block.split(":")[1].split("-")[1],
                                                                    hap[0], hap[1], hap[3], hap[2], hap[4]))
            else:
                print("%s\t%s\t%s\t%s\t%i\t%s\t%s\t%s" %  (block.split(":")[0], block.split(":")[1].split("-")[0], block.split(":")[1].split("-")[1],
                                                                    hap[0], hap[1], hap[3], "NA", hap[4]))
                
if args.toM:
    out_set = to_matrix(hap_corr)
    regionname = "window_centers_test_%s_%s.txt" % (args.org, args.chr)
    AC_name = "%s_test_K%i_count_Khap_%s.txt" % (args.chr, 15,args.org)
    AF_name = "%s_test_K%i_spectrum_Khap_%s.txt" % (args.chr, 15,args.org)
    for block in out_set.keys():
        info = out_set[block]
        #print(info)
        with open(regionname, "a") as fh:
            fh.write("%.1f\n" % (info[0][1]))
        with open(AC_name, "a") as fh:
            fh.write("%s\n" % " ".join([str(i) for i in info[1]]))
        with open(AF_name, "a") as fh:
            fh.write("%s\n" % " ".join([str(i) for i in info[2]]))
elif args.ehh:
    ehh_cal(hap_corr)

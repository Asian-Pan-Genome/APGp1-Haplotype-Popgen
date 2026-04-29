import argparse
import sys
import textwrap
import math
import numpy as np
from scipy.special import comb

parser = argparse.ArgumentParser(prog='',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
calculate
'''))
parser.add_argument('-info1', metavar='file', help='info1 file')
parser.add_argument('-info2', metavar='file', help='info2 file')
parser.add_argument('-eu', default=False, action='store_true', help='eu distance')
parser.add_argument('-simp', default=False, action='store_true', help='simp distance')
parser.add_argument('-topdiff', default=False, action='store_true', help='simp distance')
parser.add_argument('-merge', default=False, action='store_true', help='merge info files')
parser.add_argument('-Fst', default=False, action='store_true', help='calculate Fst')
parser.add_argument('-st', default=False, action='store_true', help='calculate stra') 

if len(sys.argv) == 1:
    parser.print_help()
    parser.exit()

args = parser.parse_args()

def read_info(fh):
    hap_corr = {}
    for line in fh:
        line = line.rstrip()
        line = line.replace("\n", "")
        lining = line.split("\t")
        hap = lining[2]
        af = float(lining[6])
        gene = "%s%s" % (lining[5], lining[7])
        ##NO direction
        #gene = lining[5]
         
        gene_order_pair = (lining[5], lining[7])
        samples = lining[9]
        start = lining[3]
        end = lining[4]
        if hap in hap_corr.keys():
            hap_corr[hap][0].append(gene)
            hap_corr[hap][2].append((start,end))
            hap_corr[hap][4].append(gene_order_pair)
        else:
            hap_corr[hap] = [[gene], af, [(start,end)], samples, [gene_order_pair]]
    info = {}
    for hap in hap_corr.keys():
        gene_list = tuple(hap_corr[hap][0])
        info[gene_list] = (hap_corr[hap][1], hap_corr[hap][2], hap_corr[hap][3], hap_corr[hap][4])
    return info

def merge_info(info1, info2):
    merge_info = {}
    count = 1
    hap_list = info1.keys()
    gene_info = {}
    for hap in info1:
        if hap in info2.keys():
            merge_info["Hap%i" % count] = [hap, info1[hap][0], info2[hap][0], info1[hap][2], info2[hap][2]]
            gene_info["Hap%i" % count] = [info1[hap][3], info1[hap][1], info1[hap][2], info2[hap][2]]
        else:
            merge_info["Hap%i" % count] = [hap, info1[hap][0], 0, info1[hap][2], "NA"]
            gene_info["Hap%i" % count] = [info1[hap][3], info1[hap][1], info1[hap][2], "NA"] 
        count = count + 1
    for hap in info2.keys():
        if hap not in hap_list:
            merge_info["Hap%i" % count] = [hap, 0, info2[hap][0], "NA", info2[hap][2]]
            gene_info["Hap%i" % count] = [info2[hap][3], info2[hap][1], "NA",info2[hap][2]]
    return merge_info, gene_info

def cat_stra(merge_info):
    stra = []
    info1_AC = []
    info2_AC = []
    for hap in merge_info.keys():
        info = merge_info[hap]
        if info[3] != "NA":
            AC1 = len(info[3].split("|"))
        else:
            AC1 = 0
        if info[4] != "NA":
            AC2 = len(info[4].split("|"))
        else:
            AC2 = 0
        info1_AC.append(AC1)
        info2_AC.append(AC2)
    info1_af = np.array(info1_AC) / sum(info1_AC)
    info2_af = np.array(info2_AC) / sum(info2_AC)
    i = 0
    #print(info1_af)
    #print(info2_af)
    while i < len(info1_af):
        if info1_af[i] > info2_af[i] and round(info2_af[i], 2) <= 0.1:
            stra.append(info1_af[i])
        i = i + 1
    if len(stra) == 0:
        return 0
    return max(stra)
    
def cal_fst(merge_info):
    info1_AC = []
    info2_AC = []
    for hap in merge_info.keys():
        info = merge_info[hap]
        if info[3] != "NA":
            AC1 = len(info[3].split("|"))
        else:
            AC1 = 0
        if info[4] != "NA":
            AC2 = len(info[4].split("|"))
        else:
            AC2 = 0
        info1_AC.append(AC1)
        info2_AC.append(AC2)
    if sum(info1_AC) < 2:
        Pi1 = 0
    else:
        Pi1 = 1 - sum([comb(i, 2, exact=True) for i in info1_AC]) / comb(sum(info1_AC), 2, exact=True)
    if sum(info2_AC) < 2:
        Pi2 = 0
    else:
        Pi2 = 1 - sum([comb(i, 2, exact=True) for i in info2_AC]) / comb(sum(info2_AC), 2, exact=True)
    Pixy = (Pi1 + Pi2) / 2
    if sum(info1_AC) == 0 and sum(info2_AC) == 0:
        Dxy = 0.0
    elif sum(info1_AC) == 0 or sum(info2_AC) == 0:
        Dxy = 1
    else:
        info1_af = np.array(info1_AC) / sum(info1_AC)
        info2_af = np.array(info2_AC) / sum(info2_AC)
        Dxy = 1 - sum(info1_af * info2_af)
    if Dxy == 0:
        Fst = 0.0
    else:
        if Pixy / Dxy >= 1:
            Fst = 0.0
        else:
            Fst = round(1 - Pixy / Dxy, 8)
    return Fst


def cal_euc(merge_info):
    dsum = 0
    for hap in merge_info.keys():
        info = merge_info[hap]
        #print(info[1])
        #print(info[2])
        dsum = (info[1] - info[2]) * (info[1] - info[2]) + dsum
    return math.sqrt(dsum)

def cal_diff(merge_info):
    maxhap = 0
    diff = 0
    for hap in merge_info.keys():
        info = merge_info[hap]
        if maxhap >= info[1]:
            continue
        maxhap = info[1]
        diff = info[1] - info[2]
    return diff

def cal_simp(merge_info):
    simp_sum_1 = 0
    simp_sum_2 = 0
    for hap in merge_info.keys():
        info = merge_info[hap]
        simp_sum_1 = simp_sum_1 + info[1] * info[1]
        simp_sum_2 = simp_sum_2 + info[2] * info[2]
    return 1 - simp_sum_1, 1 - simp_sum_2

with open(args.info1, 'r') as fh:
    info1 = read_info(fh)

with open(args.info2, 'r') as fh:
    info2 = read_info(fh)
#print(info1)
#print(info2)
merged, gene_info = merge_info(info1, info2)
#print(merged)
if args.eu:
    dis = cal_euc(merged)
    print("%.4f" % (dis))
elif args.simp:
    dis1, dis2 = cal_simp(merged)
    print("info1\t%.4f" % (dis1))
    print("info2\t%.4f" % (dis2))
elif args.merge:
    print("molecule\tstart\tend\tgene\tAF1\tAF2\torientation\tsample1\tsample2")
    for hapname in gene_info.keys():
        #print(merged[hapname])
        AF1 = merged[hapname][1]
        AF2 = merged[hapname][2]
        samples1 = merged[hapname][3]
        samples2 = merged[hapname][4] 
        gene_list = gene_info[hapname][0]
        region_list = gene_info[hapname][1]
        i = 0
        while i < len(gene_list):
            gene = gene_list[i]
            region = region_list[i]
            if gene[1] == "+":
                print("%s(%.2f,%.2f)\t%s\t%s\t%s\t%.4f\t%.4f\t1\t%s\t%s" % (hapname,AF1,AF2 ,region[0], region[1],gene[0],AF1,AF2,samples1,samples2))
            if gene[1] == "-":
                print("%s(%.2f,%.2f)\t%s\t%s\t%s\t%.4f\t%.4f\t0\t%s\t%s" % (hapname,AF1,AF2, region[0], region[1],gene[0],AF1,AF2,samples1,samples2))
            i = i + 1

elif args.topdiff:
    diff = cal_diff(merged)
    print("%.4f" % diff)
elif args.Fst:
    Fst = cal_fst(merged)
    print("%.4f" % Fst)
elif args.st:
    stra = cat_stra(merged)
    print("%.4f" % stra)

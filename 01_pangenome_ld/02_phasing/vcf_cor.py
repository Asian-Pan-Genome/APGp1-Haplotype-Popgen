import argparse
import sys
import textwrap
import pandas as pd
import random

parser = argparse.ArgumentParser(prog='VCF_cor',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
Description

    Correct the pangenome vcf to the original phasing state after shapit5 or beagle imputation
'''))

parser.add_argument('-org', metavar='file', help='original vcf')
parser.add_argument('-ph', metavar='file', help='phased vcf')
parser.add_argument('-l', metavar='file', help='variant list')
parser.add_argument('-c', default=False, action='store_true', help='check concordance when ignoring missing')
parser.add_argument('-m', default=False, action='store_true', help='random get missing sites')

if len(sys.argv) == 1:
    parser.print_help()
    parser.exit()

args = parser.parse_args()

def read_vcf(fh):
    header = []
    info = []
    for line in fh:
        if line.startswith('##'):
            header.append(line.strip())
            continue
        elif line.startswith('#'):
            #header.append(line.strip())
            head = line.strip().split('\t')
            continue
        fields = line.strip().split('\t')
        fields[1] = int(fields[1])
        info.append(fields)
    df = pd.DataFrame(info, columns=head)
    return header, df

def check(org, ph, sample):
    org.columns = ['pos', 'org', 'id']
    ph.columns = ['POS', 'ph', 'id']
    merged = pd.merge(org, ph, on='id', how='left')[['pos','org','ph','id']]
    total = 0
    match = 0
    for row in merged.itertuples():
        if pd.isna(row.ph):
            continue
        if row.org == ".|.":
            continue
        orgal = row.org.split("|")
        phal = row.ph.split("|")
        if orgal[0] == "." and (phal[1] != orgal[1] and phal[0] != orgal[1]):
        #if orgal[0] == "." and phal[1] != orgal[1]:
            print("%s\t%i\t%s\t%s\t%s" % (sample, row.pos, row.id, row.org, row.ph))
            continue
        if orgal[1] == "." and (phal[1] != orgal[0] and phal[0] != orgal[0]):
        #if orgal[1] == "." and phal[0] != orgal[0]:
            print("%s\t%i\t%s\t%s\t%s" % (sample, row.pos, row.id, row.org, row.ph))
            continue
        if orgal[0] != "." and orgal[1] != ".":
            if not((phal[1] == orgal[0] and phal[0] == orgal[1]) or (phal[1] == orgal[1] and phal[0] == orgal[0])):
                if phal[1] != orgal[1] or phal[0] != orgal[0]:
                    print("%s\t%i\t%s\t%s\t%s" % (sample, row.pos, row.id, row.org, row.ph))
            if phal[1] != orgal[1] and phal[0] != orgal[0]:
                if phal[1] == orgal[0] and phal[0] == orgal[1]:
                    print("%s\t%i\t%s\t%s\t%s\tsw" % (sample, row.pos, row.id, row.org, row.ph))
        #if orgal[1] != phal[1] and orgal[1] != ".":
        #    print("%s\t%i\t%s\t%s\t%s" % (sample, row.pos, row.id, row.org, row.ph))
     
def maskvialist(fh, l):
    changelog = []
    num = 0
    for line in fh:
        if line.startswith('##'):
            print(line.rstrip())
            continue
        elif line.startswith('#'):
            #header.append(line.strip())
            head = line.strip().split('\t')
            print(line.rstrip())
            num = len(head)
            continue
        fields = line.strip().split('\t')
        if fields[2] in l:
            random_number = random.randint(10, num)
            al = random.randint(0, 1)
            #print(random_number)
            #print(fields)
            orgalt = fields[random_number - 1]
            alt = fields[random_number - 1].split("|")
            alt[al] = "."
            fields[random_number - 1] = "|".join(alt)
            changelog.append("%s\t%s\t%s\t%s\t%s" % (head[random_number - 1], fields[1],fields[2], fields[random_number - 1], orgalt))
        print("\t".join(fields))
    return changelog

def fix_phased(org, ph, sample):
    org.columns = ['pos', 'org', 'id']
    ph.columns = ['POS', 'ph', 'id']
    merged = pd.merge(org, ph, on='id', how='left')[['pos','org','ph','id']]
    final = []
    strand = []
    bd = -1
    prepos = -1
    pre_dic = 0
    for row in merged.itertuples():
        if pd.isna(row.ph) or row.org == "0|0" or row.org == "1|1":
            final.append((row.pos, row.org, 0, row.id))
        elif row.org == "0|1" or row.org == "1|0":
            if row.ph == row.org:
                if pre_dic == -1:
                    strand.append((bd, prepos, -1))
                    pre_dic = 1
                    bd = row.pos
                if pre_dic == 0:
                    pre_dic = 1
                    bd = row.pos
                final.append((row.pos, row.org, 0, row.id))
            else:
                if pre_dic == 1:
                    strand.append((bd, prepos, 1))
                    pre_dic = -1
                    bd = row.pos
                if pre_dic == 0:
                    pre_dic = -1
                    bd = row.pos
                final.append((row.pos, row.org, 0, row.id))
        else:
            final.append((row.pos, row.ph, row.org, row.id))
        prepos = row.pos
    j = 0
    i = 0
    #for sd in strand:
    #    print(sd)
    while i < len(final):
        info = final[i]
        if j == len(strand):
            if info[2] == 0:
                final[i] = (info[0], info[1], 0, info[3])
            elif info[2] == ".|.":
                gt = info[1]
                final[i] = (info[0], gt, 0, info[3])
            elif info[2].split("|")[0] == ".":
                if info[2].split("|")[1] == info[1].split("|")[1]:
                    gt = info[1]
                else:
                    gt = info[1].split("|")[1] + "|" + info[2].split("|")[1]
                final[i] = (info[0], gt, "0", info[3])
            elif info[2].split("|")[1] == ".":
                if info[2].split("|")[0] == info[1].split("|")[0]:
                    gt = info[1]
                else:
                    gt = info[2].split("|")[0] + "|" + info[1].split("|")[0]
                final[i] = (info[0], gt, 0, info[3])
            i = i + 1
            continue
        inv = strand[j]
        if info[0] > strand[j][1] and info[0] > strand[j][0]:
            j += 1
        elif info[2] == 0:
            final[i] = (info[0], info[1], 0, info[3])
            i = i + 1
        elif info[0] <= strand[j][1] and info[0] >= strand[j][0] and info[2] != 0 and strand[j][2] == -1:
            if info[2] == ".|.":
                gt = info[1].split("|")[1] + "|" + info[1].split("|")[0]
                final[i] = (info[0], gt, 0, info[3])
            elif info[2].split("|")[0] == ".":
                gt = info[1].split("|")[1] + "|" + info[2].split("|")[1]
                final[i] = (info[0], gt, 0, info[3])
            elif info[2].split("|")[1] == ".":
                gt = info[2].split("|")[0] + "|" + info[1].split("|")[0]
                final[i] = (info[0], gt, 0, info[3])
            i = i + 1
        elif info[0] <= strand[j][1] and info[0] >= strand[j][0] and info[2] != 0 and strand[j][2] == 1:
            if info[2] == ".|.":
                gt = info[1]
                final[i] = (info[0], gt, 0, info[3])
            elif info[2].split("|")[0] == ".":
                gt = info[1].split("|")[0] + "|" + info[2].split("|")[1]
                final[i] = (info[0], gt, 0, info[3])
            elif info[2].split("|")[1] == ".":
                gt = info[2].split("|")[0] + "|" + info[1].split("|")[1]
                final[i] = (info[0], gt, 0, info[3])
            i = i + 1
        elif info[0] < strand[j][1] and info[0] < strand[j][0] and info[2] != 0:
            if info[2] == ".|.":
                gt = info[1]
                final[i] = (info[0], gt, 0, info[3])
            elif info[2].split("|")[0] == ".":
                if info[2].split("|")[1] == info[1].split("|")[1]:
                    gt = info[1]
                else:
                    gt = info[1].split("|")[1] + "|" + info[2].split("|")[1]
                final[i] = (info[0], gt, 0, info[3])
            elif info[2].split("|")[1] == ".":
                if info[2].split("|")[0] == info[1].split("|")[0]:
                    gt = info[1]
                else:
                    gt = info[2].split("|")[0] + "|" + info[1].split("|")[0]
                final[i] = (info[0], gt, 0, info[3])
            i = i + 1
    out = pd.DataFrame(final, columns=['pos', sample, "strand", "ID"])
    return out

if args.m:
    with open(args.l, 'r') as lfh:
        l = []
        for line in lfh:
            l.append(line.strip())
    with open(args.org, 'r') as org_fh:
        changelog = maskvialist(org_fh, l)
    with open("masking.log", 'w') as logfh:
        for log in changelog:
            logfh.write("%s\n" % log)
    sys.exit(0)

with open(args.org, 'r') as org_fh, open(args.ph, 'r') as ph_fh:
    org_header, org_df = read_vcf(org_fh)
    ph_header, ph_df = read_vcf(ph_fh)


if args.c:
    columns_index = org_df.columns.tolist()
    for sam in columns_index:
        if sam not in ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]:
            org_sample_df = org_df[["POS",sam, "ID"]]
            ph_sample_df = ph_df[["POS", sam, "ID"]]
            check(org_sample_df, ph_sample_df, sam)
    sys.exit(0)

columns_index = org_df.columns.tolist()
dfinfo = org_df[["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]]
i = 0
for sam in columns_index:
    i = i + 1
    if sam not in ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]:
        org_sample_df = org_df[["POS",sam, "ID"]]
        ph_sample_df = ph_df[["POS", sam, "ID"]]
        fixed_df = fix_phased(org_sample_df, ph_sample_df, sam)
        dfinfo = pd.merge(dfinfo, fixed_df[['ID', sam]], on='ID', how='left')

print("\n".join(org_header))


columns_index = dfinfo.columns.tolist()
ilen = len(columns_index)
print("\t".join([str(a) for a in list(dfinfo.columns)]))
for row in dfinfo.itertuples():
    outline = []
    i = 1
    while i <= ilen: 
        outline.append(str(row[i]))
        i = i + 1
    print("\t".join(outline))

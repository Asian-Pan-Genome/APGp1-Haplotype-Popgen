import argparse
import sys
import textwrap
import re

parser = argparse.ArgumentParser(prog='',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
heat to bed
'''))

parser.add_argument('-vcf', metavar='str', help='get bed')
parser.add_argument('-fam', metavar='str', help='get fam info')
parser.add_argument('-c', default=False, action='store_true', help='count missing')
parser.add_argument('-f', default=False, action='store_true', help='filter high missing')
parser.add_argument('-dip', default=False, action='store_true', help='change to dip')
parser.add_argument('-sex', default=False, action='store_true', help='process sex')
parser.add_argument('-b',  metavar='<STR>', default="0.04" ,help='input que name')
parser.add_argument('-af',  metavar='<STR>', default="0" ,help='input que name')


if len(sys.argv) == 1:
    parser.print_help()
    parser.exit()

args = parser.parse_args()


def count_missing(fh):
    for line in fh:
        line = line.rstrip()
        line = line.replace("\n","")
        if line.startswith("#"):
            continue
        lining = line.split("\t")
        i = 9
        missing = 0
        while i < len(lining):
            c = lining[i]
            if c == ".":
                missing = missing + 1
            i = i + 1
        print("%s\t%s\t%s\t%.2f" % (lining[0], lining[1], lining[2], missing/(len(lining) - 9)))


def filter_missing(fh, bound, af):
    for line in fh:
        line = line.rstrip()
        line = line.replace("\n","")
        if line.startswith("#"):
            print(line)
            continue
        lining = line.split("\t")
        alt = lining[4].split(",")
        if len(alt) >= 2:
            continue
        #alt_info = {}
        #while j <= len(alt):
            #alt_info[j] = []
        i = 9
        missing = 0
        ref_count = 0
        while i < len(lining):
            c = lining[i]
            if c == ".":
                missing = missing + 1
                lining[i] = "0"
                ref_count = ref_count + 1
            elif c == "0":
                ref_count = ref_count + 1
                #alt_info[0].append(i)
            #else:
            #    alt_info[int(c)].append(i)
            i = i + 1
        if missing/(len(lining) - 9) <= bound and ref_count/(len(lining) - 9) <= 1 - af:
            print("\t".join(lining))

def readfam(fh):
    momdic = {}
    fadic = {}
    for line in fh:
        line = line.rstrip()
        line = line.replace("\n","")
        lining = line.split("\t")
        momdic[lining[2]] = lining[0]
        fadic[lining[1]] = lining[0]
    return momdic, fadic

def todip(fh, momdic, fadic):
    for line in fh:
        line = line.rstrip()
        line = line.replace("\n","")
        if line.startswith("##"):
            print(line)
            continue
        elif line.startswith("#"):
            lining = line.split("\t")
            sample_ro = lining[9:]
            #samples = set([re.sub(r'_Mat|_Pat|-Mat|-Pat|_hap1|_hap2', '', a) for a in lining[10:]])
            #samples = list(samples)
            #if "CN1v1" in samples:
            #    samples.remove("CN1v1")
            #if "GRCh38" in samples:
            #    samples.remove("GRCh38")
            #print("%s\t%s" % ("\t".join(lining[0:9]), "\t".join(samples)))
            i = 0
            sondic={}
            while i < len(sample_ro):
                sam = sample_ro[i]
                if sam in momdic:
                    samid = momdic[sam]
                    if samid in sondic:
                        sondic[samid][1] = i
                    else:
                        sondic[samid] = [0, i]
                elif sam in fadic:
                    samid = fadic[sam]
                    if samid in sondic:
                        sondic[samid][0] = i
                    else:
                        sondic[samid] = [i, 0]
                i = i + 1
            samples = list(sondic.keys())
            print("%s\t%s" % ("\t".join(lining[0:9]), "\t".join(samples)))
            continue
        record = {}
        lining = line.split("\t")
        allele_ro = lining[9:]
        j = 0
        outinfo = []
        while j < len(samples):
            s = samples[j]
            order = sondic[s]
            patal = allele_ro[order[0]]
            matal = allele_ro[order[1]]
            outinfo.append("%s|%s" % (patal, matal))
            j = j + 1
        print("%s\t%s" % ("\t".join(lining[0:9]), "\t".join(outinfo)))

def tosex(fh):
    for line in fh:
        line = line.rstrip()
        line = line.replace("\n","")
        if line.startswith("##"):
            print(line)
            continue
        elif line.startswith("#"):
            lining = line.split("\t")
            sample_ro = lining[9:]
            samples = set([re.sub(r'_Mat|_Pat|-Mat|-Pat|_hap1|_hap2', '', a) for a in lining[10:]])
            samples = list(samples)
            print("%s\t%s" % ("\t".join(lining[0:9]), "\t".join(samples)))
            continue
        record = {}
        i = 9
        lining = line.split("\t")
        while i < len(lining):
            s = re.sub(r'_Mat|_Pat|-Mat|-Pat|_hap1|_hap2', '', sample_ro[i - 9])
            c = lining[i]
            if s in record:
                record[s].append(c)
            else:
                record[s] = [c]
            i = i + 1
        j = 0
        outinfo = []
        while j < len(samples):
            s = samples[j]
            r = record[s]
            if len(r) == 1:
                outinfo.append(r[0])
            else:
                if r[0] == "." and r[1] == ".":
                    outinfo.append(r[0])
                elif r[0] != ".": 
                    outinfo.append(r[0])
                elif r[1] != ".":
                    outinfo.append(r[1])
            j = j + 1
        print("%s\t%s" % ("\t".join(lining[0:9]), "\t".join(outinfo)))

with open(args.vcf, 'r') as fh:
    if args.c:
        count_missing(fh)
    elif args.f:
        filter_missing(fh, float(args.b), float(args.af))
    elif args.dip:
        with open(args.fam, 'r') as f:
            momdic, fadic = readfam(f)
        todip(fh, momdic, fadic)
    elif args.sex:
        tosex(fh)

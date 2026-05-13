import sys

# Read all lines into a list

id_set = {}
lines = sys.stdin.readlines()
for line in lines:
    line = line.rstrip()
    line = line.replace("\n","")
    lining = line.split()
    if len(lining) != 3:
        continue
    id = lining[1]
    chrom = lining[0]
    dic = lining[2]
    if id not in id_set.keys():
        if dic == "left":
            seq = {chrom: [1,0]}
        else:
            seq = {chrom: [0,1]}
        id_set[id] = seq
    else:
        if chrom not in id_set[id].keys():
            if dic == "left":
                id_set[id][chrom] = [1,0]
            else:
                id_set[id][chrom] = [0,1]
        else:
            if dic == "left":
                id_set[id][chrom][0] += 1
            else:
                id_set[id][chrom][1] += 1
id_filter = {}
for id in id_set.keys():
    for chrom in id_set[id]:
        chrom_filter = []
        left_count = id_set[id][chrom][0]
        right_count = id_set[id][chrom][1]
        if left_count == 1 and right_count == 1:
            chrom_filter.append(chrom)
    id_filter[id] = chrom_filter
for id in id_filter.keys():
    if len(id_filter[id]) == 1:
        print(f"{id}\t{id_filter[id][0]}")






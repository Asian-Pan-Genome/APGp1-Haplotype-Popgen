import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import tskit
import json
import os
import seaborn as sns
import pickle


import sys
import gzip
CHROM = 6
chromosome_map = [
    "INVALID",
    "NC_060925.1", # chr1
    "NC_060926.1",
    "NC_060927.1",
    "NC_060928.1",
    "NC_060929.1",
    "NC_060930.1",
    "NC_060931.1",
    "NC_060932.1",
    "NC_060933.1",
    "NC_060934.1",
    "NC_060935.1",
    "NC_060936.1",
    "NC_060937.1",
    "NC_060938.1",
    "NC_060939.1",
    "NC_060940.1",
    "NC_060941.1",
    "NC_060942.1",
    "NC_060943.1",
    "NC_060944.1",
    "NC_060945.1",
    "NC_060946.1", # chr22
]
# Ignored:
#NC_060947.1 chrX
#NC_060948.1 chrY

def gene_map_for_chr(map_file, chrom_num):
    map_result = []
    only_seq = chromosome_map[chrom_num]
    with gzip.open(map_file, "rt") as f:
        for i, line in enumerate(f):
            line = line.strip()
            if line.startswith("#"):
                continue
            line = line.split("\t")
            assert len(line) == 9
            seqname, source, feature, start, end, score, strand, frame, attribute = line
            if feature == "gene" and seqname == only_seq:
                map_result.append([int(start), int(end), attribute])
    return map_result

gene_map_6 = sorted(gene_map_for_chr("/fs/cbsubscb18/storage/Jinmin/ncbi_chm13_anno/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz", CHROM))



gtf_path = "/fs/cbsubscb18/storage/Jinmin/ncbi_chm13_anno/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz"
chrom = "NC_060930.1" #6
start, end = int(2.9e7), int(3.4e7)

records = []

def parse_attributes(attr_str):
    attrs = {}
    for item in attr_str.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        key, value = item.split(" ", 1)
        attrs[key] = value.strip('"')
    return attrs


with gzip.open(gtf_path, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue

        fields = line.rstrip().split("\t")
        if len(fields) != 9:
            continue

        chrom_name, source, feature, f_start, f_end, score, strand, frame, attributes = fields
        f_start, f_end = int(f_start), int(f_end)

        if chrom_name == chrom and f_end >= start and f_start <= end:
            record = {
                "chrom": chrom_name,
                "source": source,
                "feature": feature,
                "start": f_start,
                "end": f_end,
                "score": score,
                "strand": strand,
                "frame": frame,
            }

            record.update(parse_attributes(attributes))
            records.append(record)

df = pd.DataFrame(records)


biotype_colors = {
    "protein_coding": "black",
    "lncRNA": "#fc8803",
    "pseudogene": "blue",
    "transcribed_pseudogene": "blue",
    # TODO: look up
    "miRNA": "purple",
    "other": "purple",
    "ncRNA": "purple",
    "snRNA": "purple",
    "snoRNA": "purple",
    "C_region": "purple",
    "misc_RNA": "purple",
    "C_region_pseudogene": "purple",
    "J_segment": "purple",
    "V_segment": "purple",
    "tRNA": "grey",
    "V_segment_pseudogene": "blue",
    "NA":"red",
    "ncRNA_pseudogene":"yellow"
}


selected_genes=[
    "HLA-V",
    "HLA-A",
    "HLA-C",
    "HLA-B",
    "HLA-DRA",
    "HLA-DRB6",
    "HLA-DRB1",
    "HLA-DQA1",
    "HLA-DQB1",
    "TAP1",
    "TAP2",
    "HLA-DPB1",
    "HLA-DPA1",
    "TAPBP"
]
selected_df = df[(df['feature'] == 'gene') & (df['gene_id'].isin(selected_genes))]

selected_df = selected_df.sort_values(['chrom','start'])

selected_df["x"]=(selected_df['start']+selected_df['end'])/2
selected_df["y"]=20

selected_df["color"] = selected_df["gene_biotype"].map(biotype_colors).fillna("gray")

print(f"Number of selected genes found: {len(selected_df)}")
print(selected_df[['gene_id','chrom','start','end']])

def get_tree_intervals(ts:tskit.TreeSequence,interval=None):
    if interval is None:
        interval=[0,ts.sequence_length]
    tree=ts.at(interval[0])
    interval_pos=[]
    interval_pos.append(tree.interval[0])
    interval_pos.append(tree.interval[1]-0.01)
    while (tree.interval[1]<interval[1] and tree.interval[1]<ts.sequence_length):
        tree.next()
        interval_pos.append(tree.interval[0])
        interval_pos.append(tree.interval[1]-0.01)
        
    return np.array(interval_pos)
def get_value_per_tree(ts:tskit.TreeSequence,type,interval=None):
    if interval is None:
        interval=[0,ts.sequence_length]
    tree=ts.at(interval[0])
    doubled_values=[]
    nds=ts.nodes_time
    if type=="diversity":
        bp=[i for i in ts.breakpoints()]
        all_div=ts.diversity(windows=bp,mode="branch")
        for i,v in enumerate(all_div):
            if bp[i+1]<interval[0] or bp[i]>interval[1]:
                pass
            else:
                doubled_values.append(v)
                doubled_values.append(v)
        return np.array(doubled_values)/2
    value=None
    if type=="root":
        value=nds[tree.root]
    doubled_values.append(value)
    doubled_values.append(value)
    while (tree.interval[1]<interval[1] and tree.interval[1]<ts.sequence_length):
        tree.next()
        value=None
        if type=="root":
            value=nds[tree.root]
        doubled_values.append(value)
        doubled_values.append(value)

        
    return np.array(doubled_values)
def to_interval(a):
    extended=[]
    assert(len(a)>1)
    for i in range(len(a)):
        extended.append(a[i])
        if i==0 or i==len(a)-1:
            pass
        else:
            extended.append(a[i]+0.1)
    return np.array(extended)
def doubled(a):
    extended=[]
    for i in a:
        extended.append(i)
        extended.append(i)
    return np.array(extended)
def get_genes_for_plot(df,positions,nearby=10000,must_be_gene=True):
    if must_be_gene:
        df = df[df["feature"] == "gene"]
    dfs = []

    for pos in positions:
        left = pos - nearby
        right = pos + nearby

        mask = (df["end"] >= left) & (df["start"] <= right)
        dfs.append(df.loc[mask])

    if dfs:
        return pd.concat(dfs).drop_duplicates().reset_index(drop=True)
    else:
        return df.iloc[0:0]

step=1000
cuts = np.array(list(range(0, int(5e6) + 1, step)))


df1 = pd.DataFrame({
    "Position (BP)": to_interval(2.9e7+cuts),
    "Pairwise TMRCA (Mya)": doubled(singer_sr_div),
    "Method": "SR (Singer)",
})
df2 = pd.DataFrame({
    "Position (BP)": to_interval(2.9e7+cuts),
    "Pairwise TMRCA (Mya)": doubled(singer_lr_div),
    "Method": "LR (Singer)",
})

df3 = pd.DataFrame({
    "Position (BP)": to_interval(2.9e7+cuts),
    "Pairwise TMRCA (Mya)": doubled(tsinfer_lr_div),
    "Method": "LR (tsinfer)",
})

df4 = pd.DataFrame({
    "Position (BP)": to_interval(2.9e7+cuts),
    "Pairwise TMRCA (Mya)": doubled(tsinfer_sr_div),
    "Method": "SR (tsinfer)",
})
df5= pd.DataFrame({
    "Position (BP)": to_interval(2.9e7+cuts),
    "Pairwise TMRCA (Mya)": doubled(relate_lr_div),
    "Method": "LR (relate)",
})
df6= pd.DataFrame({
    "Position (BP)": to_interval(2.9e7+cuts),
    "Pairwise TMRCA (Mya)": doubled(relate_sr_div),
    "Method": "SR (relate)",
})


singer_df = pd.concat([df1, df2])
tsinfer_df=pd.concat([df3, df4])
relate_df=pd.concat([df5, df6])





def build_interesting_genes_panel(
    df,
    tmrcas,
    cuts,
    threshold,
    selected_df,
):

    idx = np.where(tmrcas > threshold)[0]
    interesting_pos = 2.9e7+500 + cuts[idx]

    genes = get_genes_for_plot(df, interesting_pos)
    genes["color"] = genes["gene_biotype"].map(biotype_colors).fillna("gray")

    genes = genes.copy()
    genes["x"] = (genes["start"] + genes["end"]) / 2
    genes = genes.sort_values("x").reset_index(drop=True)

    #protein coding subset
    protein = genes[genes["gene_biotype"] == "protein_coding"].copy()

    #combine with selected genes
    combined = pd.concat([
        protein,
        selected_df[~selected_df["gene_id"].isin(protein["gene_id"])]
    ]).drop_duplicates("gene_id").reset_index(drop=True)

    #randomly assign y value
    combined = random_assign_y_value(combined)

    mask = combined["gene_id"].isin(selected_df["gene_id"])
    combined.loc[mask, "color"] = "grey"

    return combined


#creating singer combined
singer_plot_data=build_interesting_genes_panel(
    df,
    singer_df,
    cuts,
    12,
    selected_df,
)
#dump
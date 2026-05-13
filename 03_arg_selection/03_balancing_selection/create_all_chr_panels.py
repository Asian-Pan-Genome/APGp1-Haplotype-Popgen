import json
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import tskit
import os
import seaborn as sns


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

def read_results_json(json_file):
    with open(json_file, "r") as f:
        data = json.load(f)
    return pd.DataFrame(data["results"])

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
    "NC_060947.1",
    "NC_060948.1"
]
chrom_name_to_id = {name: idx for idx, name in enumerate(chromosome_map) if idx != 0}
biotype_colors = {
    "protein_coding": "green",
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
    "tRNA": "purple",
    "V_segment_pseudogene": "blue",
    "NA":"red"
}
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

gtf_path = f"{SCRIPT_DIR}/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz"

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

        # Save all annotations, no filtering
        record = {
            "chr": chrom_name_to_id[chrom_name],
            "source": source,
            "feature": feature,
            "start": f_start,
            "end": f_end,
            "score": score,
            "strand": strand,
            "frame": frame,
        }

        # parse attributes into the record
        record.update(parse_attributes(attributes))
        records.append(record)

# Convert to DataFrame
df = pd.DataFrame(records)

print(df.head())
print("Total annotations:", len(df))
print("Chromosomes:", df['chr'].unique())

annotation_df = df[(df["chr"] < 23)]

def get_gene_df_by_prefix(annotation_df, gene_names, exact=False):

    if isinstance(gene_names, str):
        gene_names = [gene_names]

    if "gene_id" not in annotation_df.columns:
        raise ValueError("gene_id not found in dataframe")

    mask = annotation_df["feature"] == "gene"

    if exact:
        name_mask = annotation_df["gene_id"].isin(gene_names)
    else:
        name_mask = False
        for name in gene_names:
            name_mask |= annotation_df["gene_id"].str.startswith(name)

    return annotation_df[mask & name_mask].copy()

pass_biotypes=["protein_coding"]
def get_genes_for_plot(df,hits_df,nearby=10000,must_be_gene=True,pass_biotypes=None):
    if must_be_gene:
        df = df[df["feature"] == "gene"]
    if pass_biotypes is None:
        pass
    else:
        df = df[df["gene_biotype"].isin(pass_biotypes)]
    dfs = []

    for _, record in hits_df.iterrows():
        mid_pos = int(record["start"]) + (int(record["step"])/2)
        left = mid_pos - nearby
        right = mid_pos + nearby

        mask = (df["end"] >= left) & (df["start"] <= right) &(df["chr"]==record["chr"])
        dfs.append(df.loc[mask])


    if dfs:
        return pd.concat(dfs).drop_duplicates().reset_index(drop=True)
    else:
        return df.iloc[0:0]
    

lr_all_chr_df=pd.read_csv(f"{SCRIPT_DIR }/mask_div/1e3_lr_mask.csv")
#sr_all_chr_df=pd.read_csv(f"{SCRIPT_DIR }/mask_div/1e3_sr_mask.csv")

all_chr_df=lr_all_chr_df

# all_chr_df["value"]=all_chr_df["value"]*25/2/1e6
manhattan_df=all_chr_df.copy()
manhattan_df["chr_parity"] = manhattan_df["chr"] % 2

chr_max = (
    manhattan_df
    .groupby("chr")["start"]
    .max()
    .sort_index()
)

chr_offsets = chr_max.cumsum().shift(fill_value=0)

manhattan_df["x"] = (
    manhattan_df["start"]
    + manhattan_df["chr"].map(chr_offsets)
)


#mask
mask_manhattan_df = manhattan_df[manhattan_df["good_bp"] > 950]
mask_all_chr_df=all_chr_df[all_chr_df["good_bp"]>950]

#excluding HLA
no_hla_df = mask_all_chr_df[
    (mask_all_chr_df['chr'] != 6) | (mask_all_chr_df['start'] < 2.9e7) | (mask_all_chr_df['start'] > 3.4e7)
]

percentile_9999 = no_hla_df['value'].quantile(0.9999)
hits_df=mask_all_chr_df[mask_all_chr_df['value'] > percentile_9999]


#assign plotting y values
gene_df=get_genes_for_plot(annotation_df,hits_df,nearby=10000,pass_biotypes=pass_biotypes)
gene_df["y"]=20
gene_df["x"] = (gene_df["start"] + gene_df["end"]) / 2 + gene_df["chr"].map(chr_offsets)

#filter HLA

gene_df = gene_df[
    ~(
        (gene_df["chr"] == 6) &
        ((gene_df["x"] - (3.2e7+chr_offsets[6])).abs() < 0.3e7)
    )
]

gene_df.loc[len(gene_df)] = {
    "chr": 6,
    "x": 3.2e7+chr_offsets[6],
    "gene_id": "HLA",
    "gene_biotype": "protein_coding",
}


print(len(gene_df))

ys=10

# 11 groups, total 55
partitions = [10,10,8]+[12 for i in range(5)]+[9]
start_ys = [10 for i in range(10)]+[10]
step = 2
ys = []

for n, y0 in zip(partitions, start_ys):
    ys.extend([y0 + i * step for i in range(n)])

# assert len(ys) == len(gene_df)
gene_df["y"] = ys[:len(gene_df)]

gene_df.loc[gene_df.index[-1], "y"] = 30



mask_manhattan_df.to_csv(f"{SCRIPT_DIR}/plot_data/lr_mask_manhattan_df.csv", index=False)
gene_df.to_csv(f"{SCRIPT_DIR}/plot_data/lr_mask_gene_df.csv", index=False)

# Save plot configuration
plot_config = {
    "percentile_9999": float(percentile_9999),
    "figsize": (14, 4),
    "marker_size": 5,
    "alpha": 0.5,
    "colors": {0: 0.5, 1: 0.7},
    "ymax": 35,
    "ytick_step": 5,
    "title": "SR, 0.01% threshold excluding MHC region; protein_coding genes only"
}

with open(f"{SCRIPT_DIR}/plot_data/lr_mask_plot_config.json", "w") as f:
    json.dump(plot_config, f, indent=4)
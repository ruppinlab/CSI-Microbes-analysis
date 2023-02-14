import pandas as pd
from os.path import join


f = join("raw", "SRPRISM", "{patient}", "{sample}", "{genome}-unpaired-count.gff")

col_names = ["seqname", "source", "type", "start", "end", "score", "strand", "frame", "attributes", "read_count"]

cells = pd.read_csv(snakemake.input[0], sep="\t")
samples = cells[["patient", "sample"]].drop_duplicates()
genome = snakemake.wildcards["genome"]

gene_level_output = []

# focus on "type" == "gene" or "pseudogene"

for _, sample in samples.iterrows():
    sample_file = f.format(patient=sample["patient"], sample=sample["sample"], genome=genome)
    df = pd.read_csv(sample_file, sep="\t", names=col_names, comment="#")
    df["GeneID"] = df["attributes"].str.extract("Dbxref=GeneID:([0-9]{8})")
    df = df.loc[df["GeneID"].notna()]
    df.GeneID = df.GeneID.astype("int")
    df = df[["GeneID", "read_count"]].drop_duplicates()
    df["sample"] = sample["sample"]
    gene_level_output.append(df)


# read in the feature_df and clean it up
feature_df = pd.read_csv("data/GCF_001457555.1_NCTC10562_feature_table.txt", sep="\t")
feature_df = feature_df.rename(columns={"# feature": "feature"})
feature_df = feature_df.loc[feature_df["class"].isin(["protein_coding", "tRNA", "pseudogene", "rRNA", "SRP_RNA", "RNase_P_RNA", "tmRNA"])]
# drop two duplicated genes
feature_df = feature_df.loc[~(feature_df.GeneID.isin([45636153, 45636130]) & feature_df.name.isna())]

gene_level_df = pd.concat(gene_level_output)
gene_level_df = gene_level_df.groupby(["sample", "GeneID"])["read_count"].sum()
gene_level_df = gene_level_df.reset_index()
gene_level_df = feature_df[["feature", "class", "seq_type", "name", "symbol", "GeneID"]].merge(gene_level_df, on="GeneID")
gene_level_df.groupby(["sample", "class"])["read_count"].sum().to_csv(snakemake.output[0], sep="\t")

# gene_level_df.to_csv(snakemake.output[0], sep="\t")
# type_level_df = pd.concat(type_level_output, axis=1)
# type_level_df.to_csv(snakemake.output[1], sep="\t")

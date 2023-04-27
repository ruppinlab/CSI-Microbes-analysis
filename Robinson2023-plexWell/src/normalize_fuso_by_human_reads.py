import pandas as pd


# get the number of Fn reads
meta_df = pd.read_csv("output/P1/genus_PathSeq_All_metadata.tsv", sep="\t")
read_df = pd.read_csv("output/P1/genus_PathSeq_All_reads.tsv", sep="\t", index_col=0)

meta_df = meta_df.loc[meta_df["Condition"] == "Infected"]
fn_df = read_df.loc[848]
fn_df[meta_df.cell].sum() # 481184

# get the number of human reads
human_df = pd.read_csv("output/star/P1/human_genes.tsv", sep="\t", index_col=0)
human_df[meta_df.cell].sum().sum() # 301507801

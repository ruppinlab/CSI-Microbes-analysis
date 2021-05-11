import pandas as pd

read_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0).transpose()

meta_df = pd.read_csv(snakemake.input[1], sep="\t")
df = read_df.merge(meta_df[["sample", "patient"]], right_on="sample", left_index=True).drop(columns="sample")

df.groupby("patient").sum().to_csv(snakemake.output[0], sep="\t")

import pandas as pd

read_df1 = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
read_df2 = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)

read_df = pd.concat([read_df1, read_df2], join="outer", axis=1).fillna(0)

meta_df1 = pd.read_csv(snakemake.input[2], sep="\t")
meta_df2 = pd.read_csv(snakemake.input[3], sep="\t")
meta_df = pd.concat([meta_df1, meta_df2])

tax_df1 = pd.read_csv(snakemake.input[4], sep="\t")
tax_df2 = pd.read_csv(snakemake.input[5], sep="\t")
tax_df = pd.concat([tax_df1, tax_df2]).drop_duplicates()

read_df.to_csv(snakemake.output[0], sep="\t")
meta_df.to_csv(snakemake.output[1], sep="\t", index=False)
tax_df.to_csv(snakemake.output[2], sep="\t", index=False)

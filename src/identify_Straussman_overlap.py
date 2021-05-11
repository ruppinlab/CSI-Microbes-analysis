import pandas as pd

df = pd.read_csv(snakemake.input[0], sep="\t")

read_df = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)

hits_df = pd.read_csv(snakemake.input[2], sep="\t")

df = df.loc[df.name.isin(read_df.index)]

df.merge(hits_df, on="tax_id").to_csv(snakemake.output[0], sep="\t", index=False)

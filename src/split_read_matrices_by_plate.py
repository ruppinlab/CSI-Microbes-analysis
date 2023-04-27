import pandas as pd

read_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
metadata_df = pd.read_csv(snakemake.input[1], sep="\t")
metadata_df = metadata_df.loc[metadata_df["plate"].astype("str") == snakemake.wildcards["plate"]]
read_df = read_df[metadata_df["cell"]]
read_df.to_csv(snakemake.output[0], sep="\t")
# metadata_df.to_csv(snakemake.output[1], sep="\t", index=False)

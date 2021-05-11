import pandas as pd

read_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
#print(read_df)
metadata_df = pd.read_csv(snakemake.input[1], sep="\t")
#print(metadata_df)
metadata_df = metadata_df.loc[metadata_df.patient == snakemake.wildcards["patient"]]
#print(metadata_df)
read_df = read_df[metadata_df["cell"]]
#print(read_df.sum(axis=1))
read_df.to_csv(snakemake.output[0], sep="\t")
metadata_df.to_csv(snakemake.output[1], sep="\t", index=False)

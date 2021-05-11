import pandas as pd

hits_df = pd.read_csv(snakemake.input[0], sep="\t")
read_df = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)
read_df = read_df.loc[read_df.index.str.replace("_", " ").isin(hits_df.species)]
read_df.to_csv(snakemake.output[0], sep="\t")

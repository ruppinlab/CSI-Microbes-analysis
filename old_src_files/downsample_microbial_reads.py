import pandas as pd
import numpy as np

meta_df = pd.read_csv(snakemake.input[0], sep="\t")
read_df = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)

celltype = snakemake.wildcards["celltype"]
celltype_of_interest = snakemake.wildcards["celltype_of_interest"]
celltype_comparison = snakemake.wildcards["celltype_comparison"]
# filter by cell-types of interest
meta_df = meta_df.loc[meta_df[celltype].isin([celltype_of_interest, celltype_comparison])]
read_df = read_df[meta_df["cell"]]
# instantiate numpy random number generator
g = np.random.default_rng(seed=int(snakemake.wildcards["seed"]))
n_reads = int(snakemake.wildcards["nreads"])
# subsample read_df
read_df.loc[snakemake.wildcards["microbe"]] = g.binomial(n=n_reads, p=(read_df.loc[snakemake.wildcards["microbe"]]/read_df.loc[snakemake.wildcards["microbe"]].sum()))

meta_df.to_csv(snakemake.output[0], sep="\t", index=False)
read_df.to_csv(snakemake.output[1], sep="\t")

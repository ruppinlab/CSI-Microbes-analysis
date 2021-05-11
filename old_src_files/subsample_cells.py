import pandas as pd

meta_df = pd.read_csv(snakemake.input[0], sep="\t")
read_df = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)

celltype = snakemake.wildcards["celltype"]
celltype_of_interest = snakemake.wildcards["celltype_of_interest"]
celltype_comparison = snakemake.wildcards["celltype_comparison"]
#meta_df = meta_df.loc[meta_df.plate == snakemake.wildcards.plate]
bystander_df = meta_df.loc[meta_df[celltype] == celltype_comparison]
n_samples = min(int(snakemake.wildcards["ncells"]), bystander_df.shape[0])
bystander_df = bystander_df.sample(n=n_samples, random_state=int(snakemake.wildcards["seed"]))
infected_df = meta_df.loc[meta_df[celltype] == celltype_of_interest]
n_samples = min(int(snakemake.wildcards["ncells"]), infected_df.shape[0])
infected_df = infected_df.sample(n=n_samples, random_state=int(snakemake.wildcards["seed"]))
meta_df = pd.concat([bystander_df, infected_df])
read_df = read_df[meta_df["cell"]]

meta_df.to_csv(snakemake.output[0], sep="\t", index=False)
read_df.to_csv(snakemake.output[1], sep="\t")

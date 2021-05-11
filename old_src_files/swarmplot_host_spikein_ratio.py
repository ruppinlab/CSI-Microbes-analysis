import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

read_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
meta_df = pd.read_csv(snakemake.input[1], sep="\t")
# split by patient
meta_df = meta_df.loc[meta_df.patient == snakemake.wildcards["patient"]]
read_df = read_df[meta_df.cell.values]
spike_ins = read_df.index[read_df.index.str.startswith("ERCC-00")]
human_genes = read_df.index[~read_df.index.str.startswith("ERCC-00")]
df = pd.DataFrame(data=[read_df.loc[spike_ins].sum(), read_df.loc[human_genes].sum()]).T
df.columns = ["spike_in_reads", "human_reads"]
df = meta_df.merge(df, left_on="cell", right_index=True)
df["human_spikein_ratio"] = df["human_reads"]/df["spike_in_reads"]

fig, ax = plt.subplots(figsize=(20,10))
ax = sns.swarmplot(x=snakemake.wildcards["celltype"], y="human_spikein_ratio", data=df, ax=ax)
ax.set(ylabel="ratio of human reads to ERCC spike-in reads")

plt.legend(bbox_to_anchor=(1.1,1), loc="upper right")

plt.savefig(snakemake.output[0])

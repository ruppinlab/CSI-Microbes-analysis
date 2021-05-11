import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import ranksums

star_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
read_df = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)
meta_df = pd.read_csv(snakemake.input[2], sep="\t")
# split by patient
meta_df = meta_df.loc[meta_df.patient == snakemake.wildcards["patient"]]
star_df = star_df[meta_df.cell.values]
spike_ins = star_df.index[star_df.index.str.startswith("ERCC-00")]
human_genes = star_df.index[~star_df.index.str.startswith("ERCC-00")]
df = pd.DataFrame(data=[star_df.loc[spike_ins].sum(), star_df.loc[human_genes].sum()]).T
print(df)
df.columns = ["spike_in_reads", "human_reads"]
df = meta_df.merge(df, left_on="cell", right_index=True)
print(df)
df["human_spikein_ratio"] = df["human_reads"]/df["spike_in_reads"]
df = df.merge(read_df.T, left_on="cell", right_index=True)
print(df)
df["norm_microbe"] = np.log2((df[snakemake.wildcards["microbe"]]*df["human_spikein_ratio"])/(df["human_reads"]+df["spike_in_reads"])*1000000+1)
fig, ax = plt.subplots(figsize=(20,10))
ax = sns.swarmplot(x=snakemake.wildcards["celltype"], y="norm_microbe", data=df, ax=ax)
ax.set(ylabel="reads normalized by spike-in and sequencing depth".format(snakemake.wildcards["microbe"]))
df = df.loc[df[snakemake.wildcards["microbe"]] > 0]
print(ranksums(df.loc[df.Tumor == "Tumor"]["norm_microbe"], df.loc[df.Tumor == "nonTumor"]["norm_microbe"],))
plt.legend(bbox_to_anchor=(1.1,1), loc="upper right")

plt.savefig(snakemake.output[0])

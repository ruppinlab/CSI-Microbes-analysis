import pandas as pd
from scipy.stats import ranksums
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


df = pd.read_csv("output/metadata_for_host_transcriptome_analysis.tsv", sep="\t", index_col=0)
df = df.loc[df.infection == "infected"]

df["Myeloid"] = df.apply(lambda x: "Myeloid" if x["celltype1"] == "Myeloid" else "nonMyeloid", axis=1)

# df = df.loc[(df[read_df.columns] >= 2).any(axis=1)]
# df.groupby("celltype1").apply(lambda x: x[read_df.columns].sum(axis=1).median())
#
# myeloid_df = df.loc[df["Myeloid"] == "Myeloid"]
# non_myeloid_df = df.loc[df["Myeloid"] == "nonMyeloid"]
# ranksums(myeloid_df[read_df.columns].sum(axis=1), non_myeloid_df[read_df.columns].sum(axis=1))
#
# stromal_df = df.loc[df["Stromal"] == "Stromal"]
# non_stromal_df = df.loc[df["Stromal"] == "nonStromal"]
# ranksums(stromal_df[read_df.columns].sum(axis=1), non_stromal_df[read_df.columns].sum(axis=1))
df = df.sort_values(by=["celltype1", "celltype2"])
my_pal = {"nonMyeloid": "#4DBBD5FF", "Myeloid": "#E64B35FF"}
font = {"size": 5}
plt.rc('font', **font)
fig, ax = plt.subplots(figsize=(1.5, 2))
# n_umis_df = df[read_df.columns].sum(axis=1).to_frame("n_UMIs").merge(df[meta_df.columns], left_index=True, right_index=True)
df["log2_UMIs"] = np.log2(df["n_microbial_UMIs"])

print(ranksums(df.loc[df.celltype1 == "Myeloid"]["n_microbial_UMIs"], df.loc[df.celltype1 != "Myeloid"]["n_microbial_UMIs"]))
p = sns.violinplot(x="Myeloid", y="log2_UMIs", data=df, scale="width", order=["nonMyeloid", "Myeloid"], palette=my_pal)
ax.set(xlabel=None, ylabel=None)
plt.tight_layout()
plt.savefig("output/plots/n_umis_celltype1.pdf")
my_pal = {"Mast":"#00a087", "DC": "#4DBBD5FF", "Macrophage": "#3c5488", "Monocyte": "#E64B35FF"}
# ("#00a087", "#4DBBD5FF","#E64B35FF")
fig, ax = plt.subplots(figsize=(2, 2))
myeloid_df = df.loc[df.celltype1 == "Myeloid"]
celltypes = ["Monocyte", "Macrophage", "Mast", "DC"]
output = []
for celltype1 in celltypes:
    for celltype2 in celltypes:
        if celltype1 != celltype2:
            out = ranksums(myeloid_df.loc[myeloid_df.celltype2 == celltype1, "n_microbial_UMIs"], myeloid_df.loc[myeloid_df.celltype2 == celltype2, "n_microbial_UMIs"])
            output.append(pd.Series({"celltype1": celltype1, "celltype2": celltype2, "pval": out[1], "statistic": out[0]}))

print(pd.concat(output, axis=1).T)

order = ["Mast", "DC", "Macrophage", "Monocyte"]
p = sns.violinplot(x="celltype2", y="log2_UMIs", data=myeloid_df, scale="width", order=order, palette=my_pal)
ax.set(xlabel=None, ylabel=None)
plt.tight_layout()
plt.savefig("output/plots/n_umis_myeloid_celltype2.pdf")

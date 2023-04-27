import pandas as pd
from scipy.stats import ranksums
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


# read_df = pd.read_csv("output/infection_phenotype.tsv", sep="\t", index_col=0)
df = pd.read_csv("output/metadata_for_host_transcriptome_analysis.tsv", sep="\t", index_col=0)
#
# meta_df = pd.read_csv("data/units.tsv", sep="\t")
# meta_df = meta_df.loc[(meta_df["10x_chemistry"] == "v3") & (meta_df["Is_Tumor"] == "Yes")]
# meta_df.index = meta_df.apply(lambda x: "{}-{}".format(x["sample"], x["barcode"]), axis=1)
# df = meta_df[["patient", "sample", "barcode", "celltype1", "celltype2"]].merge(read_df, left_index=True, right_index=True)
#
# nUMI_df = pd.read_csv("output/n_microbial_UMIs_per_cell.tsv", sep="\t", index_col=0)
#
# df = df.merge(nUMI_df, left_index=True, right_index=True)
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
plt.savefig("output/plots/n_umis_myeloid.pdf")

fig, ax = plt.subplots(figsize=(2, 2))
my_pal = {"Stromal":"#00a087", "B": "orange", "TNKILC": "#4DBBD5FF", "Epithelial": "#3c5488", "Myeloid": "#E64B35FF"}
p = sns.violinplot(x="celltype1", y="log2_UMIs", data=df, scale="width", palette=my_pal)
ax.set(xlabel=None, ylabel=None)
plt.tight_layout()
plt.savefig("output/plots/n_umis_celltype1.pdf")

myeloid_df = df.loc[df.celltype1 == "Myeloid"]
celltypes = ["Monocyte", "Macrophage", "Mast", "DC", "Granulocyte"]
output = []
for celltype1 in celltypes:
    for celltype2 in celltypes:
        if celltype1 != celltype2:
            out = ranksums(myeloid_df.loc[myeloid_df.celltype2 == celltype1, "n_microbial_UMIs"], myeloid_df.loc[myeloid_df.celltype2 == celltype2, "n_microbial_UMIs"])
            output.append(pd.Series({"celltype1": celltype1, "celltype2": celltype2, "pval": out[1], "statistic": out[0]}))

print(pd.concat(output, axis=1).T)
fig, ax = plt.subplots(figsize=(2, 2))
my_pal = {"Mast":"#00a087", "Granulocyte": "orange", "DC": "#4DBBD5FF", "Macrophage": "#3c5488", "Monocyte": "#E64B35FF"}
order = ["Mast", "Granulocyte", "DC", "Macrophage", "Monocyte"]
p = sns.violinplot(x="celltype2", y="log2_UMIs", data=myeloid_df, scale="width", order=order, palette=my_pal)
ax.set(xlabel=None, ylabel=None)
plt.tight_layout()
plt.savefig("output/plots/n_umis_myeloid_celltype2.pdf")

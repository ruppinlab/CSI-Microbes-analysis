import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import ranksums

sns.set(font_scale=.5)
sns.set_style("white")
df = pd.read_csv("data/units.tsv", sep="\t")
df = df.loc[df["sample"] == "SCAF2963_3_Live"]
read_df = pd.read_csv("output/P1/genus_PathSeq_All_reads.tsv", sep="\t", index_col=0)
read_df = read_df.loc[848].to_frame("Fusobacterium")
df.index = df.apply(lambda x: "{}-{}".format(x["sample"], x["barcode"]), axis=1)
df = read_df.merge(df, left_index=True, right_index=True)
df = df.loc[df.Fusobacterium > 0]
print(ranksums(df.loc[df.celltype1 == "HCT116"]["Fusobacterium"], df.loc[df.celltype1 == "THP1"]["Fusobacterium"]))
print(ranksums(df.loc[df.celltype1 == "HCT116"]["Fusobacterium"], df.loc[df.celltype1 == "Jurkat"]["Fusobacterium"]))
print(ranksums(df.loc[df.celltype1 == "THP1"]["Fusobacterium"], df.loc[df.celltype1 == "Jurkat"]["Fusobacterium"]))
df["log2_Fusobacterium"] = np.log2(df["Fusobacterium"])
my_pal = {"Jurkat": "#4DBBD5FF", "THP1": "#3c5488", "HCT116": "#E64B35FF"}
plt.subplots(figsize=(2,2))
ax = sns.violinplot(x="celltype1", y="log2_Fusobacterium", data=df, scale="width", palette=my_pal, order=["Jurkat", "THP1", "HCT116"], saturation=1)
ax.set(xlabel=None, ylabel=None)
# plt.tick_params(axis='both', which='major', labelsize=4)
# plt.tick_params(axis='both', which='major', labelsize=4)
plt.savefig(snakemake.output[0], bbox_inches="tight")

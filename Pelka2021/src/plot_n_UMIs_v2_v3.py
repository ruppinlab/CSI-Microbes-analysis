import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import ranksums


df = pd.read_csv("data/units.tsv", sep="\t")
df = df.loc[df["Is_Tumor"] == "Yes"]

tax_level = "genus"
read_file = "output/{}/{}/{}_PathSeq_All_reads.tsv"

sample_df = df[["patient", "sample", "10x_chemistry"]].drop_duplicates()

v2_sample_df = sample_df.loc[sample_df["10x_chemistry"] == "v2"]

output = []
for x, row in v2_sample_df.iterrows():
    read_df = pd.read_csv(read_file.format(row["patient"], row["sample"], tax_level), sep="\t", index_col=0)
    output.append(read_df)

v2_df = pd.concat(output, axis=1).fillna(0)
v2_df = v2_df.sum()


v3_sample_df = sample_df.loc[sample_df["10x_chemistry"] == "v3"]

output = []
for x, row in v3_sample_df.iterrows():
    read_df = pd.read_csv(read_file.format(row["patient"], row["sample"], tax_level), sep="\t", index_col=0)
    output.append(read_df)

v3_df = pd.concat(output, axis=1).fillna(0)
v3_df = v3_df.sum()

v2_df = v2_df.loc[v2_df > 0]
v3_df = v3_df.loc[v3_df > 0]
print(ranksums(v2_df.values, v3_df.values))

v2_df = v2_df.to_frame("n_UMIs")
v3_df = v3_df.to_frame("n_UMIs")

v2_df["chemistry"] = "3' v2"
v3_df["chemistry"] = "3' v3"
df = pd.concat([v2_df, v3_df])
df["log2_UMIs"] = np.log2(df["n_UMIs"])

sns.set(font_scale=.5)
sns.set_style("white")

my_pal = {"3' v2": "#4DBBD5FF", "3' v3": "#E64B35FF"}
plt.subplots(figsize=(2,2))
ax = sns.violinplot(x="chemistry", y="log2_UMIs", data=df, scale="width", palette=my_pal, order=["3' v2", "3' v3"], saturation=1)
ax.set(xlabel=None, ylabel=None)
# plt.tick_params(axis='both', which='major', labelsize=4)
# plt.tick_params(axis='both', which='major', labelsize=4)
plt.savefig("output/plots/v2_v3_nUMIs_comparison.svg", bbox_inches="tight")

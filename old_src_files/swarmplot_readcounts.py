import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

read_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
metadata_df = pd.read_csv(snakemake.input[1], sep="\t")
metadata_df = metadata_df.loc[metadata_df.patient == snakemake.wildcards["patient"]]
print(read_df)
print(metadata_df.cell.values)
read_df = read_df[metadata_df.cell.values]
readcount_df = read_df.sum().to_frame("total_reads")
readcount_df["log2(reads+1)"] = np.log2(readcount_df["total_reads"]+1)
print(readcount_df)
df = readcount_df.merge(metadata_df, left_index=True, right_on="cell")
fig, ax = plt.subplots(figsize=(20,10))
ax = sns.swarmplot(x=snakemake.wildcards["celltype"], y="log2(reads+1)", data=df, ax=ax)
plt.legend(bbox_to_anchor=(1.1,1), loc="upper right")
# ax.legend(loc='upper right', bbox_to_anchor=(1.03, 1))

#test_results
plt.savefig(snakemake.output[0])
plt.close(fig)

# fig, ax = plt.subplots(figsize=(20,10))
# ax = sns.boxplot(x="patient", y="contaminants", hue=snakemake.wildcards["celltype"], data=contam_df, ax=ax)
# test_results = add_stat_annotation(ax, data=contam_df, x="patient", box_pairs=box_pairs,
#                                    y="contaminants", hue=snakemake.wildcards["celltype"], test="Mann-Whitney", text_format="star", loc="outside")
# plt.legend(bbox_to_anchor=(1.1,1), loc="upper right")
#
# #test_results
# plt.savefig(snakemake.output[0])
# plt.close(fig)

# df.groupby("cancer").boxplot(column=snakemake.wildcards["microbe"], subplots=False, sharex=True)
# fig1, ax1 = plt.subplots(figsize=(20, 10))
# ax1.boxplot([read_df[tumor_samples],
#              read_df[nontumor_samples]],
#             labels=["Tumor", "Non-Tumor"])
# ax1.set_xlabel("cell type")
# ax1.set_ylabel("number of unambiguous reads per cell")
# ax1.set_title("Number of reads from {} across cell types".format(snakemake.wildcards["microbe"]))
#plt.show()

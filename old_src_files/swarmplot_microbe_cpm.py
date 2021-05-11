import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

read_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
metadata_df = pd.read_csv(snakemake.input[1], sep="\t")
metadata_df = metadata_df.loc[metadata_df.patient == snakemake.wildcards["patient"]]
read_df = read_df[metadata_df.cell.values]
read_df.loc["read_count"] = read_df.sum()
print(read_df)
print(metadata_df)
# log2 normalize the readcounts for plotting
microbe_log2_read_df = np.log2((read_df.loc[snakemake.wildcards["microbe"]]/read_df.loc["read_count"])*1000000+1)
# print(microbe_log2_read_df)
microbe_log2_read_df.name = snakemake.wildcards["microbe"]
print(microbe_log2_read_df)
# drop celltypes with unknown celltype
# metadata_df = metadata_df.loc[metadata_df[snakemake.wildcards["celltype"]] != "unknown"]
# drop patients with only one celltype
#celltypes_per_patient = metadata_df.groupby("patient")[snakemake.wildcards["celltype"]].nunique()
#patients_to_keep = celltypes_per_patient.loc[celltypes_per_patient.values > 1].index
#metadata_df = metadata_df.loc[metadata_df.patient.isin(patients_to_keep)]
# read_df = read_df.loc[read_df.index.str.startswith(snakemake.wildcards["microbe"])].sum().to_frame(name=snakemake.wildcards["microbe"])
df = metadata_df.merge(microbe_log2_read_df, right_index=True, left_on="cell")
fig, ax = plt.subplots(figsize=(20,10))
ax = sns.swarmplot(x=snakemake.wildcards["celltype"], y=snakemake.wildcards["microbe"], data=df, ax=ax)
ax.set(ylabel="log2({} cpm+1)".format(snakemake.wildcards["microbe"]))
plt.legend(bbox_to_anchor=(1.1,1), loc="upper right")

plt.savefig(snakemake.output[0])

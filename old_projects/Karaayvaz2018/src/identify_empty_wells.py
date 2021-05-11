import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

human_read_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
meta_df = pd.read_csv(snakemake.input[1], sep="\t")
microbe_read_df = pd.read_csv(snakemake.input[2], sep="\t", index_col=0)

num_human_reads = human_read_df.sum().to_frame(name="total_human_reads")
num_human_genes = (human_read_df > 0).astype(int).sum().to_frame(name="num_detected_human_genes")
human_df = pd.concat([num_human_reads, num_human_genes], axis=1)
human_df["log2(total_human_reads+1)"] = human_df["total_human_reads"].apply(lambda x: np.log2(x+1))
human_df["log2(num_detected_human_genes+1)"] = human_df["num_detected_human_genes"].apply(lambda x: np.log2(x+1))

num_microbe_reads = microbe_read_df.sum().to_frame(name="total_microbe_reads")
num_microbe_genes = (microbe_read_df > 0).astype(int).sum().to_frame(name="num_detected_microbial_genera")
microbe_df = pd.concat([num_microbe_reads, num_microbe_genes], axis=1)
microbe_df["log2(total_microbe_reads+1)"] = microbe_df["total_microbe_reads"].apply(lambda x: np.log2(x+1))
microbe_df["log2(num_detected_microbial_genera+1)"] = microbe_df["num_detected_microbial_genera"].apply(lambda x: np.log2(x+1))


# filter by the three plates of interest
meta_df = meta_df.loc[meta_df.plate.isin(["B1-P1", "B5-P1", "B6-P7"])]
# we have 252 rows but we should have 288 (96*3) - weird

# human_df["well_status"] = human_df["log2(num_detected_human_genes+1)"].apply(lambda x: "empty" if x < 6 else "not-empty")

df = human_df.merge(microbe_df, left_index=True, right_index=True)
df = df.merge(meta_df, left_index=True, right_index="sample")

fig, ax = plt.subplots(figsize=(20,10))
sns.scatterplot(x="log2(num_detected_microbial_genera+1)", y="log2(num_detected_human_genes+1)", hue="QC_status", data=df, ax=ax)
plt.savefig(snakemake.output[2])
plt.close(fig)

microbe_df.to_csv(snakemake.output[1], sep="\t")
meta_df.to_csv(snakemake.output[0], sep="\t", index=False)

# let's focus only on B1-P1 (PT039 - 24 empty wells), B5-P1 (PT084 - 12 empty wells) and B6-P7 (PT089 - 11 empty wells)

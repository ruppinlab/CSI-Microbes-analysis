import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv("data/units.tsv", sep="\t")
df = df.loc[df["Is_Tumor"] == "Yes"]

tax_level = "genus"
read_file = "output/{}/{}/{}_PathSeq_All_reads.tsv"

sample_df = df[["patient", "sample", "10x_chemistry"]].drop_duplicates()

v2_sample_df = sample_df.loc[sample_df["10x_chemistry"] == "v2"]

output = []
for x, row in v2_sample_df.iterrows():
    read_df = pd.read_csv(read_file.format(row["patient"], row["sample"], tax_level), sep="\t", index_col=0)
    output.append(read_df.sum(axis=1))

v2_summed_reads = pd.concat(output, axis=1).fillna(0).sum(axis=1)

v3_sample_df = sample_df.loc[sample_df["10x_chemistry"] == "v3"]

output = []
for x, row in v3_sample_df.iterrows():
    read_df = pd.read_csv(read_file.format(row["patient"], row["sample"], tax_level), sep="\t", index_col=0)
    output.append(read_df.sum(axis=1))

v3_summed_reads = pd.concat(output, axis=1).fillna(0).sum(axis=1)


tax_map_file = "output/{}/tax_id_map_All_PathSeq.tsv"
output = []
for p in df.patient.unique():
    df = pd.read_csv(tax_map_file.format(p), sep="\t", index_col=0)
    output.append(df)

tax_df = pd.concat(output).drop_duplicates()
tax_df = tax_df.loc[tax_df["taxa_level"] == tax_level]
d = dict(zip(tax_df["tax_id"], tax_df.index))
v2_summed_reads = v2_summed_reads.rename(index=d)
v3_summed_reads = v3_summed_reads.rename(index=d)

# define for piecharts
colors = sns.color_palette()

# generate piechart for v2
v2_df = v2_summed_reads.to_frame("nreads")
v2_df = v2_df.sort_values(by="nreads", ascending=False)
top_v2_df = v2_df.iloc[0:9].copy()
top_v2_df.loc["Other"] = v2_df.iloc[9:].sum()
top_v2_df["percentage"] = top_v2_df/top_v2_df.sum()

top_v2_df["name"] = top_v2_df.index
legend_labels = top_v2_df.apply(lambda x: "{} {:.0%}".format(x["name"], x["percentage"]), axis=1)

plt.subplots(figsize=(1.8, 1.8)) # , rotatelabels=True
plt.pie(top_v2_df["nreads"], colors=colors, startangle=90)

plt.legend(labels = legend_labels, bbox_to_anchor=(1,1), prop={'size': 4})
plt.savefig("output/plots/v2_genera_piechart.svg", bbox_inches="tight")

# generate piechart for v3
v3_df = v3_summed_reads.to_frame("nreads")
v3_df = v3_df.sort_values(by="nreads", ascending=False)
top_v3_df = v3_df.iloc[0:9].copy()
top_v3_df.loc["Other"] = v3_df.iloc[9:].sum()
top_v3_df["percentage"] = top_v3_df/top_v3_df.sum()

top_v3_df["name"] = top_v3_df.index
legend_labels = top_v3_df.apply(lambda x: "{} {:.0%}".format(x["name"], x["percentage"]), axis=1)

plt.subplots(figsize=(1.8, 1.8)) # , rotatelabels=True
plt.pie(top_v3_df["nreads"], colors=colors, startangle=90)

plt.legend(labels = legend_labels, bbox_to_anchor=(1,1), prop={'size': 4})
plt.savefig("output/plots/v3_genera_piechart.svg", bbox_inches="tight")

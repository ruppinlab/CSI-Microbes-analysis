import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# df = pd.read_csv("output/Pt0/S0/Fn-read-type-table.tsv")
df = pd.read_csv(snakemake.input[0], sep="\t")

live_5p_df = df.loc[df["sample"] == "SCAF2965_5_Live"].copy()

live_5p_df["percentage"] = live_5p_df["read_count"]/live_5p_df["read_count"].sum()
live_5p_df = live_5p_df.sort_values(by="percentage", ascending=False)
other_df = live_5p_df.loc[live_5p_df.percentage < .001]
top_df = live_5p_df.loc[live_5p_df.percentage > .001]

other_df = other_df.set_index("class").sum()
other_df["class"] = "Other"
top_df = pd.concat([top_df, other_df.to_frame().T])

# piechart_labels = ["Fusobacterium"] + [None]*9
legend_labels = top_df.apply(lambda x: "{} {:.0%}".format(x["class"], x["percentage"]), axis=1)
colors = sns.color_palette()

plt.subplots(figsize=(2, 2))
plt.pie(top_df["read_count"], colors=colors, startangle=90)
plt.legend(labels = legend_labels, bbox_to_anchor=(1,1), prop={'size': 4})
plt.savefig(snakemake.output[0], bbox_inches="tight")

live_3p_df = df.loc[df["sample"] == "SCAF2963_3_Live"].copy()

live_3p_df["percentage"] = live_3p_df["read_count"]/live_3p_df["read_count"].sum()
live_3p_df = live_3p_df.sort_values(by="percentage", ascending=False)
other_df = live_3p_df.loc[live_3p_df.percentage < .001]
top_df = live_3p_df.loc[live_3p_df.percentage > .001]

other_df = other_df.set_index("class").sum()
other_df["class"] = "Other"
top_df = pd.concat([top_df, other_df.to_frame().T])

# piechart_labels = ["Fusobacterium"] + [None]*9
legend_labels = top_df.apply(lambda x: "{} {:.0%}".format(x["class"], x["percentage"]), axis=1)
colors = sns.color_palette()

plt.subplots(figsize=(2, 2))
plt.pie(top_df["read_count"], colors=colors, startangle=90)
plt.legend(labels = legend_labels, bbox_to_anchor=(1,1), prop={'size': 4})
plt.savefig(snakemake.output[1], bbox_inches="tight")

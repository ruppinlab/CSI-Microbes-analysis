import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


read_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
tax_map = pd.read_csv(snakemake.input[2], sep="\t")

read_df = tax_map[["name", "tax_id"]].merge(read_df, left_on="tax_id", right_index=True)
read_df = read_df.set_index("name")
read_df = read_df.drop(columns="tax_id")
meta_df = pd.read_csv(snakemake.input[1], sep="\t")

meta_df = meta_df.loc[meta_df["infected"].isin(["infected", "bystander"])]

read_df = read_df[meta_df.cell]

summed_df = read_df.sum(axis=1).sort_values(ascending=False)
# take the top 9 rows
top_df = summed_df.iloc[0:9]
top_df.loc["Other"] = summed_df.iloc[9:].sum()


top_df = top_df.to_frame(name="reads")
top_df = top_df.reset_index()

top_df["percentage"] = top_df["reads"]/top_df["reads"].sum()

piechart_labels = ["Salmonella"] + [None]*9
legend_labels = top_df.apply(lambda x: "{} {:.0%}".format(x["name"], x["percentage"]), axis=1)
explode = [.01] + [0]*9
colors = sns.color_palette()
red = colors[3]
blue = colors[0]
orange = colors[1]
colors[0] = red
colors[1] = blue
colors[3] = orange


plt.subplots(figsize=(1.8,1.8))
plt.pie(top_df["reads"], explode=explode, colors=colors, startangle=90, labels=piechart_labels, labeldistance=0.1, textprops=dict(color="w", size=5))
plt.legend(labels = legend_labels, bbox_to_anchor=(1,1), prop={'size': 4})
# plt.title("Aulicino2018 (plate-based)")
# (1,0) - bottom middle of the graph
# (0,0) - bottom left of the graph
# (0,1) - top left of the graph - this is what I want
# plt.tight_layout()
# plt.subplots_adjust(left=.2)
# plt.show()
plt.savefig(snakemake.output[0], bbox_inches="tight")

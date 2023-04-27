import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


read_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
tax_map = pd.read_csv(snakemake.input[1], sep="\t")

read_df = tax_map[["name", "tax_id"]].merge(read_df, left_on="tax_id", right_index=True)
read_df = read_df.set_index("name")
read_df = read_df.drop(columns="tax_id")
meta_df = pd.read_csv(snakemake.input[2], sep="\t")

meta_df = meta_df.loc[meta_df["sample"] == "GSM3454529"]

read_df = read_df[meta_df.cell]

summed_df = read_df.sum(axis=1).sort_values(ascending=False)
# take the top 9 rows
top_df = summed_df.iloc[0:9]
top_df.loc["Other"] = summed_df.iloc[9:].sum()


top_df = top_df.to_frame(name="reads")
top_df = top_df.reset_index()
top_df.to_csv("output/Pt0/summed_microbe_reads.tsv", sep="\t")


top_df = pd.read_csv("output/Pt0/summed_microbe_reads.tsv", sep="\t")
top_df["percentage"] = top_df["reads"]/top_df["reads"].sum()

piechart_labels = [None] + ["Salmonella"] + [None]*8
legend_labels = top_df.apply(lambda x: "{} {:.0%}".format(x["name"], x["percentage"]), axis=1)
explode = [0] + [0.01] + [0]*8
colors = sns.color_palette()
red = colors[3]
blue = colors[0]
orange = colors[1]
colors[1] = red
colors[0] = blue
colors[3] = orange


plt.subplots(figsize=(1.8,1.8)) # , rotatelabels=True
plt.pie(top_df["reads"], explode=explode, colors=colors, rotatelabels=True,
                         startangle=90, labeldistance=0.1, labels=piechart_labels,
                         textprops=dict(color="w", size=5))
# for t in texts:
#     t.set_rotation(75)
    # t.set_verticalalignment('center')
plt.legend(labels = legend_labels, bbox_to_anchor=(1,1), prop={'size': 4})
# plt.title("Ben-Moshe2019 (10x 3' v2)")
# (1,0) - bottom middle of the graph
# (0,0) - bottom left of the graph
# (0,1) - top left of the graph - this is what I want
# plt.tight_layout()
# plt.subplots_adjust(left=.2)
#plt.show()
plt.savefig(snakemake.output[0], bbox_inches="tight")

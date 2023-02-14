import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

top_df = pd.read_csv("output/P1/summed_microbe_reads.tsv", sep="\t")
top_df["percentage"] = top_df["reads"]/top_df["reads"].sum()


piechart_labels = ["Fusobacterium"] + [None]*9
legend_labels = top_df.apply(lambda x: "{} {:.0%}".format(x["name"], x["percentage"]), axis=1)
explode = [.01] + [0]*9
colors = sns.color_palette()
red = colors[3]
blue = colors[0]
orange = colors[1]
colors[0] = red
colors[1] = blue
colors[3] = orange


plt.subplots(figsize=(1.8, 1.8))
plt.pie(top_df["reads"], explode=explode, colors=colors, startangle=90, labels=piechart_labels, labeldistance=0.05, textprops=dict(color="w", size=5))
plt.legend(labels = legend_labels, bbox_to_anchor=(1,1), prop={'size': 4})
# plt.title("Robinson2023-plexWell (plate-based)")
# plt.legend(labels = labels, bbox_to_anchor=(0,1))
# (1,0) - bottom middle of the graph
# (0,0) - bottom left of the graph
# (0,1) - top left of the graph - this is what I want
# plt.tight_layout()
# plt.subplots_adjust(left=.2)
plt.savefig("output/plots/piechart.svg", bbox_inches="tight")

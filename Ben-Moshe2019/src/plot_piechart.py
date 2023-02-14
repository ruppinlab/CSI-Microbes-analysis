import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

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
plt.savefig("output/plots/piechart.svg", bbox_inches="tight")

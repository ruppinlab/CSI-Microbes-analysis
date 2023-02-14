import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv(snakemake.input[0], sep="\t")

df["percentage"] = df["read_count"]/df["read_count"].sum()
df = df.sort_values(by="percentage", ascending=False)
other_df = df.loc[df.percentage < .001]
top_df = df.loc[df.percentage > .001]

other_df = other_df.set_index("class").sum()
other_df["class"] = "Other"
top_df = pd.concat([top_df, other_df.to_frame().T])

# piechart_labels = ["Fusobacterium"] + [None]*9
legend_labels = top_df.apply(lambda x: "{} {:.0%}".format(x["class"], x["percentage"]), axis=1)
# explode = [.01] + [0]*9
colors = sns.color_palette()
# red = colors[3]
# blue = colors[0]
# orange = colors[1]
# colors[0] = red
# colors[1] = blue
# colors[3] = orange


plt.subplots(figsize=(2, 2))
plt.pie(top_df["read_count"], colors=colors, startangle=90)
plt.legend(labels = legend_labels, bbox_to_anchor=(1,1), prop={'size': 4})
# plt.title("Robinson2023-10x (10x 3' v3 and 5' v2)")
# (1,0) - bottom middle of the graph
# (0,0) - bottom left of the graph
# (0,1) - top left of the graph - this is what I want
# plt.tight_layout()
# plt.subplots_adjust(left=.2)
# plt.show()
plt.savefig(snakemake.output[0], bbox_inches="tight")

# taken from https://stackoverflow.com/questions/51908395/plot-circles-where-size-depends-number-of-votes-python
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as mcolors
import matplotlib.cm as cm

pelka_df = pd.read_csv("../Pelka2021/output/expected_actual_genera_infected_greater_2umis_cells_per_celltype1.tsv", sep="\t")
zhang_df = pd.read_csv("../Zhang2021/output/expected_actual_genera_infected_greater_2umis_cells_per_celltype1.tsv", sep="\t")

pelka_num_infected_cells = pelka_df.groupby("microbe")["actual_infected_cells"].sum()

zhang_num_infected_cells = zhang_df.groupby("microbe")["actual_infected_cells"].sum()

zhang_top_genera = zhang_num_infected_cells.loc[zhang_num_infected_cells >= 10].index
pelka_top_genera = pelka_num_infected_cells.loc[pelka_num_infected_cells >= 10].index

top_genera = list(set(list(zhang_top_genera) + list(pelka_top_genera)))

pelka_myeloid_df = pelka_df.loc[pelka_df["celltype1"] == "Myeloid"]
zhang_myeloid_df = zhang_df.loc[zhang_df["celltype1"] == "Myeloid"]

pelka_top_myeloid_df = pelka_myeloid_df.loc[pelka_myeloid_df.microbe.isin(top_genera)].copy()
pelka_other_myeloid_df = pelka_myeloid_df.loc[~pelka_myeloid_df.microbe.isin(top_genera)].copy()
pelka_other_myeloid_df = pelka_other_myeloid_df[["expected_infected_cells", "actual_infected_cells"]].sum()
pelka_other_myeloid_df["log2FC"] = np.log2(pelka_other_myeloid_df["actual_infected_cells"]/pelka_other_myeloid_df["expected_infected_cells"])
pelka_other_myeloid_df["microbe"] = "Other"
pelka_other_myeloid_df["celltype1"] = "Myeloid"
pelka_myeloid_df = pd.concat([pelka_top_myeloid_df, pelka_other_myeloid_df.to_frame().T])
pelka_myeloid_df["cancer_type"] = "Pelka2021"
# repeat for zhang
zhang_top_myeloid_df = zhang_myeloid_df.loc[zhang_myeloid_df.microbe.isin(top_genera)].copy()
zhang_other_myeloid_df = zhang_myeloid_df.loc[~zhang_myeloid_df.microbe.isin(top_genera)].copy()
zhang_other_myeloid_df = zhang_other_myeloid_df[["expected_infected_cells", "actual_infected_cells"]].sum()
zhang_other_myeloid_df["log2FC"] = np.log2(zhang_other_myeloid_df["actual_infected_cells"]/zhang_other_myeloid_df["expected_infected_cells"])
zhang_other_myeloid_df["microbe"] = "Other"
zhang_other_myeloid_df["celltype1"] = "Myeloid"
zhang_myeloid_df = pd.concat([zhang_top_myeloid_df, zhang_other_myeloid_df.to_frame().T])
zhang_myeloid_df["cancer_type"] = "Zhang2021"

df = pd.concat([pelka_myeloid_df, zhang_myeloid_df])
df["log2FC"] = df.log2FC.replace(-np.Inf, np.nan)
df["log2FC"] = df.log2FC.fillna(0)
df = df[["microbe", "log2FC", "cancer_type"]]

# get the number of infected cells
pelka_num_infected_cells = pelka_num_infected_cells.reset_index()
pelka_num_infected_cells["cancer_type"] = "Pelka2021"
other_infected_cells = pelka_num_infected_cells.loc[~pelka_num_infected_cells.microbe.isin(top_genera)]["actual_infected_cells"].sum()
pelka_num_infected_cells = pelka_num_infected_cells.loc[pelka_num_infected_cells.microbe.isin(top_genera)]
other_df = pd.Series({"cancer_type": "Pelka2021", "microbe": "Other", "actual_infected_cells": other_infected_cells}).to_frame().T
pelka_num_infected_cells = pd.concat([pelka_num_infected_cells, other_df])

zhang_num_infected_cells = zhang_num_infected_cells.reset_index()
zhang_num_infected_cells["cancer_type"] = "Zhang2021"
other_infected_cells = zhang_num_infected_cells.loc[~zhang_num_infected_cells.microbe.isin(top_genera)]["actual_infected_cells"].sum()
zhang_num_infected_cells = zhang_num_infected_cells.loc[zhang_num_infected_cells.microbe.isin(top_genera)]
other_df = pd.Series({"cancer_type": "Zhang2021", "microbe": "Other", "actual_infected_cells": other_infected_cells}).to_frame().T
zhang_num_infected_cells = pd.concat([zhang_num_infected_cells, other_df])

# normalize for the number of cells sequenced - 171052 for Zhang and 90312 for Pelka
pelka_num_infected_cells["actual_infected_cells"] = pelka_num_infected_cells["actual_infected_cells"]*(171052/90312)
n_infected_cells = pd.concat([pelka_num_infected_cells, zhang_num_infected_cells])
df = df.merge(n_infected_cells, on=["microbe", "cancer_type"])
df["log2_infected_cells"] = np.log2(df["actual_infected_cells"].astype("float")+1)
df = df.sort_values(by="microbe", ascending=False)

# update so "Other" is last
other_df = df.loc[df.microbe == "Other"]
df = df.loc[df.microbe != "Other"]
df = pd.concat([other_df, df])

# mcolors.CenteredNorm() - keeps the center at zero
# cmap="RdBu" uses a divergent colormap, which is ideal for showing our results
# setup the colorbar
# can we use something other than ScalarMappable?
scalarmappaple = cm.ScalarMappable(norm=mcolors.CenteredNorm(), cmap="RdBu_r")
scalarmappaple.set_array(df["log2FC"])

#sc = plt.scatter(x, y, s=a2, alpha=0.5)
#plt.legend(*sc.legend_elements("sizes", num=6))

font = {"size": 6}
plt.rc('font', **font)
fig, ax = plt.subplots(figsize=(3, 4.1))
sc = ax.scatter(x=df["cancer_type"], y=df["microbe"], s=df["log2_infected_cells"]*15, c=df["log2FC"], norm=mcolors.CenteredNorm(), cmap="RdBu_r")
plt.legend(*sc.legend_elements("sizes", num=4), bbox_to_anchor=(0, -.05), loc='upper left', ncols=3)

plt.margins(x=.5)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="10%", pad=0.5)
# scatter_ax = divider.append_axes("left", size="90%", pad=0.25)
# first goal, let's tilt the genera to be on the x-axis and 45 degrees

# scatter_ax.xticks(rotation=45)

# plt.colorbar -
plt.colorbar(scalarmappaple, cax=cax)
plt.tight_layout() # this needs to come after plt.colorbar
# plt.show()
# plt.legend()
# plt.show()
plt.savefig("output/plots/genera_overview.pdf")

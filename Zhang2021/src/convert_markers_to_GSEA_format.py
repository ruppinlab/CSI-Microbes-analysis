import pandas as pd
import numpy as np
from math import log10


df = pd.read_csv("output/infected_uninfected_tumor_cell_subset_markers.tsv", sep="\t")

df.index.name = "Name" # this is required
df["GSEA_statistic"] = df.apply(lambda x: np.sign(x["avg_log2FC"])*-log10(x["p_val"]+1e-300), axis=1)
df[["GSEA_statistic", "avg_log2FC", "p_val_adj"]].sort_values(by="GSEA_statistic", ascending=False).to_csv("output/infected_uninfected_tumor_cell_subset_markers.rnk", sep="\t")

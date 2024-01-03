import pandas as pd
from statsmodels.stats.multitest import fdrcorrection

f = "output/microbial_load_corr/celltype1_{}_{}_load_corr.tsv"
celltype = "TNKILC"
# genera = ["Treponema", "Leptotrichia", "Bacteroides", "Capnocytophaga", "Fusobacterium", "Campylobacter", "Streptococcus"]
genera = ["Pseudomonas", "Sphingomonas", "Mycoplasma", "Treponema", "Leptotrichia", "Shigella", "Bacteroides", "Capnocytophaga", "Fusobacterium", "Alloprevotella", "Campylobacter", "Streptococcus"]

output = []
for g in genera:
    try:
        df = pd.read_csv(f.format(celltype, g), sep="\t")
        df = df.loc[~df.pval.isna()]
        df["p_val_adj"] = fdrcorrection(df["pval"])[1]
        df["genus"] = g
        output.append(df)
    except Exception as e:
        print("missing file for {} and {}".format(celltype, g))

df = pd.concat(output)
df = df.loc[df.p_val_adj < .05]
df["celltype"] = celltype

df.to_csv("output/microbial_load_corr/{}_load_corr.tsv".format(celltype), sep="\t")
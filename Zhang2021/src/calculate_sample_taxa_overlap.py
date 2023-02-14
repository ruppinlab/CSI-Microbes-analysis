import sys

import pandas as pd
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import hypergeom

min_umis = int(snakemake.wildcards["min_umis"])

# read_df = pd.read_csv(snakemake.input[0], sep="\t")
read_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)

# meta_df = pd.read_csv(snakemake.input[1], sep="\t")
meta_df = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)
#meta_df = meta_df.loc[meta_df["celltype1"] == "Myeloid"]
#print(meta_df)
# tax_df = pd.read_csv(snakemake.input[2], sep="\t")
tax_df = pd.read_csv(snakemake.input[2], sep="\t")

read_df = read_df.merge(tax_df[["tax_id", "name"]], left_index=True, right_on="tax_id").set_index("name").drop(columns="tax_id")
read_df = read_df.T
# read_df = read_df.loc[meta_df.index]
top_bacteria = read_df.columns[(read_df >= min_umis).sum() > 5]
if len(top_bacteria) <= 1:
    pd.DataFrame(data={"b1": [], "b2": [], "pval": [], "b1_n": [], "b2_n": [], "overlap_n": [], "fdr": []}).to_csv(snakemake.output[0], sep="\t", index=False)
    sys.exit()
output = []
for b1 in top_bacteria:
    for b2 in top_bacteria:
            b1_n = read_df.loc[read_df[b1] >= min_umis].shape[0]
            b2_n = read_df.loc[read_df[b2] >= min_umis].shape[0]
            overlap_n = read_df.loc[(read_df[b1] >= min_umis) & (read_df[b2] >= min_umis)].shape[0]
            #print("co-occurrence p-value for {} and {}".format(b1, b2))
            pval = hypergeom.sf(overlap_n-1, read_df.shape[0], b1_n, b2_n)
            output.append(pd.Series(data={"b1": b1, "b2": b2, "pval": pval, "b1_n": b1_n, "b2_n": b2_n, "overlap_n": overlap_n}))

overlap_df = pd.concat(output, axis=1).T
overlap_df = overlap_df.drop_duplicates(subset=["b1", "b2"])
overlap_df = overlap_df.loc[overlap_df["b2"] > overlap_df["b1"]]
overlap_df["fdr"] = fdrcorrection(overlap_df["pval"], method="i")[1]
overlap_df = overlap_df.sort_values(by="pval")
print(overlap_df)
overlap_df.to_csv(snakemake.output[0], sep="\t", index=False)

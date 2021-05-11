import pandas as pd
from statsmodels.stats.multitest import fdrcorrection


def combine_results(x):
    # figure out the direction
    up_x = x.loc[x.direction == "up"]
    down_x = x.loc[x.direction == "down"]
    #print(up_x["p.value"].mean())
    #print(down_x["p.value"].mean())
    if up_x["p.value"].mean() < down_df["p.value"].mean():
        #print("choosing up")
        x = up_x
    else:
        x = down_x
    return pd.Series({"p.value": min(x["p.value"].multiply(2).mean(), 1), "summary.AUC": x["summary.AUC"].mean()})
    # need to convert two-sided p-value to one-sided p-value

gene_list = []
up_output = []
down_output = []
for i in snakemake.input["up"]:
    df = pd.read_csv(i, sep="\t", index_col=0)
    df = df.rename_axis("gene").reset_index()[["gene", "p.value", "summary.AUC"]]
    up_output.append(df)
    gene_list.append(df.gene)

for i in snakemake.input["down"]:
    df = pd.read_csv(i, sep="\t", index_col=0)
    df = df.rename_axis("gene").reset_index()[["gene", "p.value", "summary.AUC"]]
    down_output.append(df)
    gene_list.append(df.gene)

overlapped_genes = set(gene_list[0]).intersection(*gene_list)
up_df = pd.concat(up_output)
up_df["direction"] = "up"
down_df = pd.concat(down_output)
down_df["direction"] = "down"
df = pd.concat([up_df, down_df])
df = df.loc[df.gene.isin(overlapped_genes)]

# now let's perform aggregation
new_df = df.groupby("gene").apply(combine_results)
print(fdrcorrection(new_df["p.value"]))
new_df["FDR"] = fdrcorrection(new_df["p.value"])[1]
new_df = new_df.sort_values(by="FDR")
new_df.to_csv(snakemake.output[0], sep="\t")
#print(new_df)

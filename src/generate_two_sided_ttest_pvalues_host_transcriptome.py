import pandas as pd
from statsmodels.stats.multitest import fdrcorrection


def combine_results(x):
    #print(x)
    # figure out the direction
    up_x = x.loc[x.direction == "up",]
    down_x = x.loc[x.direction == "down",]
    #print(up_x["p.value"].median())
    #print(down_x["p.value"].median())
    if up_x["p.value"].median() < down_df["p.value"].median():
        # use the results from the direction="up" test
        # minimum AUC should .5
        x = up_x
        auc = max(x["summary.AUC"].median(), .5)
    else:
        # use the results from direction="down" test
        # maximum AUC should be .5
        x = down_x
        auc = min(x["summary.AUC"].median(), .5)
    # multiply by 2 to convert from one-sided to two-sided (but p-values can't be greater than 1)
    return pd.Series({"p.value": min(x["p.value"].multiply(2).median(), 1), "summary.AUC": auc})
    # need to convert two-sided p-value to one-sided p-value

#ind = snakemake.params[["index"]]

#up_output = []
#down_output = []
#for i in snakemake.input["up"]:
up_df = pd.read_csv(snakemake.input["up"], sep="\t", index_col=0)
up_df = up_df.reset_index()[["index", "p.value", "summary.AUC"]]
# up_output.append(df)
    #gene_list.append(df.gene)

#for i in snakemake.input["down"]:
down_df = pd.read_csv(snakemake.input["down"], sep="\t", index_col=0)
down_df = down_df.reset_index()[["index", "p.value", "summary.AUC"]]
# down_output.append(df)
    #gene_list.append(df.gene)


#overlapped_genes = set(gene_list[0]).intersection(*gene_list)
#up_df = pd.concat(up_output)
up_df["direction"] = "up"
#down_df = pd.concat(down_output)
down_df["direction"] = "down"
df = pd.concat([up_df, down_df])
#df = df.loc[df.gene.isin(overlapped_genes)]

#auc_df = pd.concat(aucs_output).groupby("gene").median()
# now let's perform aggregation
new_df = df.groupby("index").apply(combine_results)
#print(fdrcorrection(new_df["p.value"]))
new_df["FDR"] = fdrcorrection(new_df["p.value"])[1]
new_df = new_df.sort_values(by="FDR")
new_df.index.name = None
new_df.index = new_df.index.astype("str")
# new_df = new_df.set_index("index")
#new_df = new_df.merge(auc_df, left_index=True, right_index=True)
new_df.to_csv(snakemake.output[0], sep="\t")

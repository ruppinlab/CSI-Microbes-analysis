import pandas as pd
from scipy.stats import spearmanr, kendalltau, combine_pvalues, pearsonr


spearman_rho=pd.NamedAgg(column="spearman_rho", aggfunc="mean"),
spearman_pval=pd.NamedAgg(column="spearman_pval", aggfunc=(lambda x: combine_pvalues(x)[1])),
pearson_rho=pd.NamedAgg(column="pearson_rho", aggfunc="mean"),
pearson_pval=pd.NamedAgg(column="pearson_pval", aggfunc=(lambda x: combine_pvalues(x)[1])),
kendall_rho=pd.NamedAgg(column="kendall_rho", aggfunc="mean"),
kendall_pval=pd.NamedAgg(column="kendall_pval", aggfunc=(lambda x: combine_pvalues(x)[1])),

def aggregate_results(df, microbe):
    summed_microbe_counts = df["summed_microbe_counts"].sum()
    percentage_zero_cells = df["percentage_zero_cells"].mean()
    spearman_rho = df["spearman_rho"].mean()
    pearson_rho = df["pearson_rho"].mean()
    kendall_rho = df["kendall_rho"].mean()
    spearman_pvals = []
    pearson_pvals = []
    kendall_pvals = []
    for index, row in df.iterrows():
        if row["spearman_rho"] > 0:
            spearman_pvals.append(row["spearman_pval"]/2)
        else:
            spearman_pvals.append(1-(row["spearman_pval"]/2))
        if row["pearson_rho"] > 0:
            pearson_pvals.append(row["pearson_pval"]/2)
        else:
            pearson_pvals.append(1-(row["pearson_pval"]/2))
        if row["kendall_rho"] > 0:
            kendall_pvals.append(row["kendall_pval"]/2)
        else:
            kendall_pvals.append(1-(row["kendall_pval"]/2))
    print(spearman_pvals)
    spearman_pval = combine_pvalues(spearman_pvals)[1]
    pearson_pval = combine_pvalues(pearson_pvals)[1]
    kendall_pval = combine_pvalues(kendall_pvals)[1]
    d = {
        "microbe": microbe,
        "summed_microbe_count": summed_microbe_counts,
        "spearman_pval": spearman_pval,
        "spearman_rho": spearman_rho,
        "pearson_pval": pearson_pval,
        "pearson_rho": pearson_rho,
        "kendall_pval": kendall_pval,
        "kendall_rho": kendall_rho,
        "percentage_zero_cells": percentage_zero_cells
    }
    return pd.Series(data=d)


microbe_read_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
spikein_read_df = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)
meta_df = pd.read_csv(snakemake.input[2], sep="\t")
taxid_df = pd.read_csv(snakemake.input[3], sep="\t")
print(taxid_df)
print(taxid_df["tax_id"])
# we only care about the total number of spike-ins
spikein_read_df = spikein_read_df.sum().to_frame("spikein_reads")

df = spikein_read_df.merge(microbe_read_df.T, left_index=True, right_index=True)

df = df.merge(meta_df, left_index=True, right_on="cell")

plates = df.plate.unique()
microbes = microbe_read_df.index

output = []
for microbe in microbes:
    for plate in plates:
        plate_df = df.loc[df.plate == plate]
        plate_microbe_df = plate_df[microbe]
        spearman_res = spearmanr(plate_df["spikein_reads"], plate_microbe_df)
        pearson_res = pearsonr(plate_df["spikein_reads"], plate_microbe_df)
        kendall_res = kendalltau(plate_df["spikein_reads"], plate_microbe_df)
        d = {
            "taxid": microbe,
            "spearman_pval": spearman_res[1],
            "spearman_rho": spearman_res[0],
            "pearson_pval": pearson_res[1],
            "pearson_rho": pearson_res[0],
            "kendall_pval": kendall_res[1],
            "kendall_rho": kendall_res[0],
            "summed_microbe_counts": plate_microbe_df.sum(),
            "percentage_zero_cells": plate_microbe_df.loc[plate_microbe_df == 0].shape[0]/plate_microbe_df.shape[0]
            }
        s = pd.Series(data=d)
        output.append(s)

# print(output)
results_df = pd.concat(output, axis=1).T
results_df = results_df.astype(dtype={
    "taxid": "int64",
    "spearman_rho": "float64",
    "spearman_pval": "float64",
    "percentage_zero_cells": "float64",
    "pearson_pval": "float64",
    "pearson_rho": "float64",
    "kendall_pval": "float64",
    "kendall_rho": "float64",
})

results_df["taxid"] = results_df["taxid"].astype("str")
taxid_df["tax_id"] = taxid_df["tax_id"].astype("str")
taxid_df = taxid_df.set_index("tax_id")
results_df["microbe"] = results_df["taxid"].apply(lambda x: taxid_df.at[x, "name"])

results_df.to_csv(snakemake.output[0], sep="\t", index=False)
# results are per-plate
# now, let's aggregate the results
output = []
grouped = results_df.groupby("microbe")
for name, group in grouped:
    s = aggregate_results(group, name)
    # print(s)
    output.append(s)


agg_results_df = pd.concat(output, axis=1).T.sort_values(by="summed_microbe_count", ascending=False)
print(agg_results_df)
# grouped = results_df.groupby("microbe")
# agg_results_df = grouped.agg(
#     microbe_count=pd.NamedAgg(column="summed_microbe_counts", aggfunc="sum"),
#     spearman_rho=pd.NamedAgg(column="spearman_rho", aggfunc="mean"),
#     spearman_pval=pd.NamedAgg(column="spearman_pval", aggfunc=(lambda x: combine_pvalues(x)[1])),
#     pearson_rho=pd.NamedAgg(column="pearson_rho", aggfunc="mean"),
#     pearson_pval=pd.NamedAgg(column="pearson_pval", aggfunc=(lambda x: combine_pvalues(x)[1])),
#     kendall_rho=pd.NamedAgg(column="kendall_rho", aggfunc="mean"),
#     kendall_pval=pd.NamedAgg(column="kendall_pval", aggfunc=(lambda x: combine_pvalues(x)[1])),
#     percentage_zero_cells=pd.NamedAgg(column="percentage_zero_cells", aggfunc="mean")
# )


agg_results_df.to_csv(snakemake.output[1], sep="\t", index=False)

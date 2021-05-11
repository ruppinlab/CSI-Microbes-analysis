import pandas as pd
from scipy.stats import fisher_exact

def calculate_fisher_exact(x, cells_of_interest, comparison_cells):
    print(x)
    positive_celltype_of_interest = sum(x[cells_of_interest] > 0)
    negative_celltype_of_interest = sum(x[cells_of_interest] == 0)
    positive_comparison_celltype = sum(x[comparison_cells] > 0)
    negative_comparison_celltype = sum(x[comparison_cells] == 0)
    mat = [[positive_celltype_of_interest, negative_celltype_of_interest], [positive_comparison_celltype, negative_comparison_celltype]]
    oddsratio, pvalue = fisher_exact(mat)
    print(oddsratio)
    print(pvalue)
    return oddsratio, pvalue

read_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)

meta_df = pd.read_csv(snakemake.input[1], sep="\t")

celltype_column = snakemake.wildcards["celltype"]
celltype_of_interest = snakemake.wildcards["celltype_of_interest"]
comparison_celltype = snakemake.wildcards["celltype_comparison"]

cells_of_interest = meta_df.loc[meta_df[celltype_column] == celltype_of_interest]["cell"].values

if comparison_celltype == "all":
    comparison_cells = meta_df.loc[meta_df[celltype_column] != celltype_of_interest]["cell"].values
else:
    comparison_cells = meta_df.loc[meta_df[celltype_column] == comparison_celltype]["cell"].values

out = read_df.apply(calculate_fisher_exact, args=(cells_of_interest, comparison_cells), axis=1)
print(out)

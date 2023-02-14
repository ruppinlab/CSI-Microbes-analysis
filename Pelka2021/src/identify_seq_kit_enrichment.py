import pandas as pd
from scipy.stats import hypergeom


# loop through the superkingdom reads
# 2 = "Bacteria"
# 2759 = "Eukaryota"
# 10239 = "Viruses"
files = snakemake.input["microbe_reads"]
output = []
for f in files:
    df = pd.read_csv(f, sep="\t", index_col=0)
    output.append(df)
# concat merges on the index by default
read_df = pd.concat(output, join="outer", axis=1).fillna(0)
read_df = read_df.T

# get the tax_dfs
files = snakemake.input["taxid_maps"]
output = []
for f in files:
    df = pd.read_csv(f, sep="\t", index_col=0)
    output.append(df)

tax_df = pd.concat(output).drop_duplicates()
tax_df = tax_df.loc[tax_df["taxa_level"] == snakemake.wildcards["tax_level"]]
d = dict(zip(tax_df["tax_id"], tax_df.index))
read_df = read_df.rename(columns=d)

meta_df = pd.read_csv(snakemake.input["meta_data"], sep="\t")
meta_df.index = meta_df.apply(lambda x: "{}-{}".format(x["sample"], x["barcode"]), axis=1)
df = read_df.merge(meta_df, left_index=True, right_index=True)

group = "10x_chemistry"#snakemake.wildcards["group"]
group1 = "v2"
group2 = "v3"
min_umis = int(snakemake.wildcards["min_umis"])

df = df.loc[df["Is_Tumor"] == "No"]
group1_df = df.loc[df[group] == group1]
n_group1_cells = group1_df.shape[0]
group2_df = df.loc[df[group] == group2]
n_group2_cells = group2_df.shape[0]
M = n_group1_cells + n_group2_cells
output = []

for taxa in read_df.columns:
    # print(taxa)
    #  # total number of reads
    total_positives = df.loc[df[taxa] >= min_umis].shape[0]
    if total_positives == 0:
        continue
    n_taxa_group1_positives = group1_df.loc[group1_df[taxa] >= min_umis].shape[0]
    n_taxa_group2_positives = group2_df.loc[group2_df[taxa] >= min_umis].shape[0]
    # M - population size - all normal cells
    # n - successes in the population - normal cells with >= min_umis of the particular taxa
    # N - sample size - normal cells sequenced by v2
    # X - number of drawn successes - normal cells sequenced by v2 with >= min_umis of the particular taxa
    # pval = hypergeom.sf(x-1, M, n, N) - from https://alexlenail.medium.com/understanding-and-implementing-the-hypergeometric-test-in-python-a7db688a7458
    pval = hypergeom.sf(n_taxa_group1_positives-1, M, total_positives, n_group1_cells)
    output.append(pd.DataFrame(data={group1: [n_taxa_group1_positives/n_group1_cells], group2: [n_taxa_group2_positives/n_group2_cells], "p-value": [pval], "taxa": [taxa]}))

df = pd.concat(output)

df.sort_values(by="p-value").to_csv(snakemake.output[0], sep="\t", index=False)
    # print(out_df[taxa]/total_reads)
# out_df.name = "Percent.Positive.Cells"
# out_df = out_df.reset_index()
# print(out_df.sort_values(by="Percent.Positive.Cells"))
# sns.boxplot(x="10x_chemistry", y="Percent.Positive.Cells", hue="Is_Tumor", data=out_df)
# plt.savefig(snakemake.output[0])
# repeat but stratified by Is_Tumor == "Yes" vs. Is_Tumor == "No" and 10x_chemistry

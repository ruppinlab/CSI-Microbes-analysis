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

group = snakemake.wildcards["group"]
group1 = snakemake.wildcards["group1"]
group2 = snakemake.wildcards["group2"]

out_df = df.groupby([group])[read_df.columns].sum()
total_reads = out_df.sum(axis=1)
M = total_reads.sum()
#print(total_reads)
total_group1_reads = total_reads[group1]
total_group2_reads = total_reads[group2]
#print(total_group1_reads)
#print(total_group2_reads)
output = []
for taxa in out_df.columns:
    print(taxa)
    #  # total number of reads
    taxa_reads = out_df[taxa].sum()
    taxa_group1_reads = out_df.loc[group1][taxa] # total_reads[snakemake.wildcards["group1"]]
    taxa_group2_reads = out_df.loc[group2][taxa] # total_reads[snakemake.wildcards["group2"]]
    #print(taxa_group1_reads)
    #print(taxa_group2_reads)
    hpd = hypergeom(M, taxa_reads, total_group1_reads)
    p = hpd.pmf(taxa_group1_reads)
    output.append(pd.DataFrame(data={group1: [taxa_group1_reads/total_group1_reads], group2: [taxa_group2_reads/total_group2_reads], "p-value": [p], "taxa": [taxa]}))

df = pd.concat(output)
df.to_csv(snakemake.output[0], sep="\t", index=False)
    # print(out_df[taxa]/total_reads)
# out_df.name = "Percent.Positive.Cells"
# out_df = out_df.reset_index()
# print(out_df.sort_values(by="Percent.Positive.Cells"))
# sns.boxplot(x="10x_chemistry", y="Percent.Positive.Cells", hue="Is_Tumor", data=out_df)
# plt.savefig(snakemake.output[0])
# repeat but stratified by Is_Tumor == "Yes" vs. Is_Tumor == "No" and 10x_chemistry

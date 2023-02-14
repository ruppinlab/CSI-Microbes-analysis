import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


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
# read_df = read_df.drop(columns=["Mycoplasma", "Mycoplasma_wenyonii"], errors="ignore")

meta_df = pd.read_csv(snakemake.input["meta_data"], sep="\t")
meta_df.index = meta_df.apply(lambda x: "{}-{}".format(x["sample"], x["barcode"]), axis=1)
meta_df["Tissue"] = meta_df.apply(lambda x: "Normal" if "N" in (x["sample"].split("-")[0]) else "Tumor", axis=1)
meta_df["cell_population"] = meta_df.apply(lambda x: x["sample"].split("-")[1], axis=1)
df = read_df.merge(meta_df, left_index=True, right_index=True)
min_umis = int(snakemake.wildcards["min_umis"])
# microbe = snakemake.wildcards["microbe_of_interest"]
# pd.set_option('display.max_rows', 500)
# print(df.loc[df.index.str.startswith("C106_N")].shape)
# normal_df = df.loc[df.index.str.startswith("C106_N")][read_df.columns]
# print(normal_df.loc[(normal_df >= min_umis).any(axis=1)])
# print(df.loc[(df.patient == "C106") & (df["Is_Tumor"] == "No")])
# normal_df = df.loc[(df.patient == "C106") & (df["Is_Tumor"] == "No")][read_df.columns]
# print(normal_df.loc[(normal_df >= min_umis).any(axis=1)].shape[0])
#df = df.loc[df["patient"] == "C106"]
# calculate the % of positive cells (depending on the wildcard min_umis) per patient stratified by Is_Tumor == "Yes" vs. Is_Tumor == "No"
#print(read_df)
#print(df[read_df.columns])
#print(df.groupby(["patient", "Is_Tumor", "10x_chemistry"])[read_df.columns])
#df.groupby(["patient", "Is_Tumor", "10x_chemistry"]).apply(lambda x: print(x.loc[(x[read_df.columns] > 0).any(axis=1)]))
#print(df.groupby(["patient", "Is_Tumor", "10x_chemistry"]).apply(lambda x: x.loc[(x[read_df.columns] >= min_umis).any(axis=1)].shape[0]))
#print(df.groupby(["patient", "Is_Tumor", "10x_chemistry"]).apply(lambda x: x.shape[0]))
out_df = df.groupby(["patient", "cell_population", "Tissue"]).apply(lambda x: (x.loc[(x[read_df.columns] >= min_umis).any(axis=1)].shape[0]/x.shape[0])*100)
out_df.to_csv(snakemake.output[1], sep="\t")
name = "% Cells with >= {} UMIs".format(min_umis)
out_df.name = name
out_df = out_df.reset_index()
print(out_df.sort_values(by=name))
sns.boxplot(x="cell_population", y=name, hue="Tissue", data=out_df)
plt.savefig(snakemake.output[0])
# repeat but stratified by Is_Tumor == "Yes" vs. Is_Tumor == "No" and 10x_chemistry

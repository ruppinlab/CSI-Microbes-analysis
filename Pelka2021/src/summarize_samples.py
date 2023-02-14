import pandas as pd
# from scipy.stats import hypergeom


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
patient_df = pd.read_csv("data/patients.tsv", sep="\t")
sample_df = pd.read_csv("data/samples.tsv", sep="\t")
meta_df = meta_df.merge(patient_df[["patient", "MMRStatus"]], on="patient").merge(sample_df[["sample", "ProcessingMethod"]].drop_duplicates(), on="sample")
meta_df.index = meta_df.apply(lambda x: "{}-{}".format(x["sample"], x["barcode"]), axis=1)
df = read_df.merge(meta_df, left_index=True, right_index=True)

min_umis = int(snakemake.wildcards["min_umis"])
# only focus on cells from the tumor
df = df.loc[df["Is_Tumor"] == "Yes"]
# count the number of cells per sample with >= 2 UMIs for each taxa
sample_df = df.groupby(["patient", "MMRStatus", "sample", "ProcessingMethod", "10x_chemistry"]).apply(lambda x: (x[read_df.columns] >= 2).sum())
# focus on taxa with at least 10 cells with >= 2 UMIs
sample_df = sample_df[sample_df.columns[sample_df.sum() > 10]]
# drop Mycoplasma
sample_df = sample_df.drop(columns="Mycoplasma")

sample_df.sum(axis=1).sort_values().iloc[-20:]

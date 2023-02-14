import pandas as pd

df = pd.read_csv("data/units.tsv", sep="\t")
df = df.loc[df["Is_Tumor"] == "Yes"]

tax_level = "genus"
read_file = "output/{}/{}_PathSeq_All_reads.tsv"
output = []
for p in df.patient.unique():
    output.append(pd.read_csv(read_file.format(p, tax_level), sep="\t", index_col=0))

read_df = pd.concat(output, join="outer", axis=1).fillna(0)
read_df = read_df.T

meta_df = pd.read_csv("data/units.tsv", sep="\t")
meta_df.index = meta_df.apply(lambda x: "{}-{}".format(x["sample"], x["barcode"]), axis=1)

tax_map_file = "output/{}/tax_id_map_All_PathSeq.tsv"
output = []
for p in df.patient.unique():
    df = pd.read_csv(tax_map_file.format(p), sep="\t", index_col=0)
    output.append(df)

tax_df = pd.concat(output).drop_duplicates()
tax_df = tax_df.loc[tax_df["taxa_level"] == tax_level]
d = dict(zip(tax_df["tax_id"], tax_df.index))
read_df = read_df.rename(columns=d)

df = meta_df.merge(read_df, left_index=True, right_index=True)
df = df.loc[df["Is_Tumor"] == "Yes"]
read_df = read_df.loc[df.index]

# read_df["max_reads"] = read_df.max(axis=1)
# read_df[["max_reads", "Bacteroides", "Fusobacterium"]].to_csv("output/cohort_reads_for_plot.tsv", sep="\t")
read_df.to_csv(snakemake.output[0], sep="\t")
df = df.rename(columns={"10x_chemistry": "chemistry"})
# df[["patient", "sample", "barcode", "chemistry"]].to_csv("output/cohort_metadata_for_plot.tsv", sep="\t")
df[["patient", "sample", "barcode", "chemistry"]].to_csv(snakemake.output[1], sep="\t")

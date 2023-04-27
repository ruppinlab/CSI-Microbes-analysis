import pandas as pd
from numpy import log2

units = pd.read_csv("data/units.tsv", sep="\t")
units = units.loc[units["Is_Tumor"] == "Yes"]
df_v3 = units.loc[units["10x_chemistry"] == "v3"]
sample_df_v3 = df_v3[["patient", "sample"]].drop_duplicates()
SAMPLE_MICROBE_READ_TABLE = join("output", "{patient}", "{sample}", "{tax_level}_{method}_{kingdom}_reads.tsv")
files = [SAMPLE_MICROBE_READ_TABLE.format(patient=row.patient, sample=row["sample"], tax_level="genus", method="PathSeq", kingdom="All") for _, row in sample_df_v3.iterrows()]


files = snakemake.input["microbe_reads"]
output = []
for f in files:
    df = pd.read_csv(f, sep="\t", index_col=0)
    output.append(df)

# concat merges on the index by default
read_df = pd.concat(output, join="outer", axis=1).fillna(0)

PATIENT_PATHSEQ_TAXID_MAP = join("output", "{patient}", "tax_id_map_{kingdom}_{method}.tsv")
patient_df_v3 = df_v3[["patient"]].drop_duplicates()
tax_files = [PATIENT_PATHSEQ_TAXID_MAP.format(patient=row.patient, method="PathSeq", kingdom="All") for _, row in patient_df_v3.iterrows()]
# add the tax id information
tax_files = snakemake.input["tax_ids"]

output = []
for f in tax_files:
    output.append(pd.read_csv(f, sep="\t"))

tax_map = pd.concat(output).drop_duplicates()

read_df = read_df.merge(tax_map, left_index=True, right_on="tax_id")
read_df = read_df.set_index("name")
read_df = read_df.drop(columns=["tax_id", "taxa_level"])

read_df = read_df.T
genera = read_df.columns
meta_df = pd.read_csv("data/units.tsv", sep="\t")
meta_df = meta_df.loc[(meta_df["10x_chemistry"] == "v3") & (meta_df["Is_Tumor"] == "Yes")]
meta_df.index = meta_df.apply(lambda x: "{}-{}".format(x["sample"], x["barcode"]), axis=1)
df = meta_df[["patient", "sample", "barcode", "celltype1", "celltype2"]].merge(read_df, left_index=True, right_index=True)

df["infection"] = "uninfected"
df.loc[(df[genera] >= 1).any(axis=1), "infection"] = "bystander" # change from bystander
df.loc[(df[genera] >= 2).any(axis=1), "infection"] = "infected"

# drop all the genera columns now
df = df[["patient", "sample", "barcode", "celltype1", "celltype2", "infection"]]

df.groupby(["patient"]).apply(lambda x: x.loc[x.infection == "infected"].shape[0]/x.shape[0]).mean()

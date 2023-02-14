import pandas as pd
from os.path import join


df = pd.read_csv("data/units.tsv", sep="\t")

# let's just focus on tumor samples
df = df.loc[df["Is_Tumor"] == "Yes"]
samples_with_myeloid_cells = df.loc[df["celltype1"] == "Myeloid"][["patient", "sample"]].drop_duplicates()
samples_with_tumor_cells = df.loc[df["celltype1"] == "Epi"][["patient", "sample"]].drop_duplicates()
samples_with_B_cells = df.loc[df["Bcell"] == "Bcell"][["patient", "sample"]].drop_duplicates()
samples_with_TNKILC_cells = df.loc[df["celltype1"] == "TNKILC"][["patient", "sample"]].drop_duplicates()
samples_with_Strom_cells = df.loc[df["celltype1"] == "Strom"][["patient", "sample"]].drop_duplicates()
samples_with_Mast_cells = df.loc[df["celltype1"] == "Mast"][["patient", "sample"]].drop_duplicates()

SAMPLE_hFDR_FISHER_MARKERS = join("output", "results", "{}", "{}", "hFDR-fisher-exact-{}-{}-{}-PathSeq-All-0-1.tsv")

output = []

for _, row in samples_with_myeloid_cells.iterrows():
    res_df = pd.read_csv(SAMPLE_hFDR_FISHER_MARKERS.format(row["patient"], row["sample"], "Myeloid", "Myeloid", "nonMyeloid"), sep="\t")
    res_df["patient"] = row["patient"]
    res_df["sample"] = row["sample"]
    output.append(res_df)

myeloid_res_df = pd.concat(output)
myeloid_res_df["celltype"] = "Myeloid"
myeloid_res_df["celltype_of_interest"] = "Myeloid"
myeloid_res_df["celltype_comparison"] = "nonMyeloid"

output = []

for _, row in samples_with_tumor_cells.iterrows():
    res_df = pd.read_csv(SAMPLE_hFDR_FISHER_MARKERS.format(row["patient"], row["sample"], "Tumor", "Tumor", "nonTumor"), sep="\t")
    res_df["patient"] = row["patient"]
    res_df["sample"] = row["sample"]
    output.append(res_df)

tumor_res_df = pd.concat(output)
tumor_res_df["celltype"] = "Tumor"
tumor_res_df["celltype_of_interest"] = "Tumor"
tumor_res_df["celltype_comparison"] = "nonTumor"

output = []

for _, row in samples_with_B_cells.iterrows():
    res_df = pd.read_csv(SAMPLE_hFDR_FISHER_MARKERS.format(row["patient"], row["sample"], "Bcell", "Bcell", "nonBcell"), sep="\t")
    res_df["patient"] = row["patient"]
    res_df["sample"] = row["sample"]
    output.append(res_df)

Bcell_res_df = pd.concat(output)
Bcell_res_df["celltype"] = "Bcell"
Bcell_res_df["celltype_of_interest"] = "Bcell"
Bcell_res_df["celltype_comparison"] = "nonBcell"

output = []

for _, row in samples_with_TNKILC_cells.iterrows():
    res_df = pd.read_csv(SAMPLE_hFDR_FISHER_MARKERS.format(row["patient"], row["sample"], "TNKILC", "TNKILC", "nonTNKILC"), sep="\t")
    res_df["patient"] = row["patient"]
    res_df["sample"] = row["sample"]
    output.append(res_df)

TNKILC_res_df = pd.concat(output)
TNKILC_res_df["celltype"] = "TNKILC"
TNKILC_res_df["celltype_of_interest"] = "TNKILC"
TNKILC_res_df["celltype_comparison"] = "nonTNKILC"

output = []

for _, row in samples_with_Strom_cells.iterrows():
    res_df = pd.read_csv(SAMPLE_hFDR_FISHER_MARKERS.format(row["patient"], row["sample"], "Strom", "Strom", "nonStrom"), sep="\t")
    res_df["patient"] = row["patient"]
    res_df["sample"] = row["sample"]
    output.append(res_df)

Strom_res_df = pd.concat(output)
Strom_res_df["celltype"] = "Strom"
Strom_res_df["celltype_of_interest"] = "Strom"
Strom_res_df["celltype_comparison"] = "nonStrom"

output = []

for _, row in samples_with_Mast_cells.iterrows():
    res_df = pd.read_csv(SAMPLE_hFDR_FISHER_MARKERS.format(row["patient"], row["sample"], "Mast", "Mast", "nonMast"), sep="\t")
    res_df["patient"] = row["patient"]
    res_df["sample"] = row["sample"]
    output.append(res_df)

Mast_res_df = pd.concat(output)
Mast_res_df["celltype"] = "Mast"
Mast_res_df["celltype_of_interest"] = "Mast"
Mast_res_df["celltype_comparison"] = "nonMast"

df = pd.concat([myeloid_res_df, tumor_res_df, Bcell_res_df, TNKILC_res_df, Strom_res_df, Mast_res_df])

df.to_csv(snakemake.output[0], sep="\t", index=False)

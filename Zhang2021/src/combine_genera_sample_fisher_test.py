import pandas as pd
from os.path import join
from statsmodels.stats.multitest import fdrcorrection


cells = pd.read_csv("data/units.tsv", sep="\t")

# let's just focus on tumor samples
cells = cells.loc[cells["sample"].str.contains("T")]
cells = cells.loc[~cells.patient.isin(["P38", "P39", "P40"])]

# cd45pos_samples = cells.loc[cells["selection"] == "CD45pos"]#[["patient", "sample"]].drop_duplicates()
# cd45neg_samples = cells.loc[cells["selection"] == "CD45neg"]#[["patient", "sample"]].drop_duplicates()
cd45pos_samples = cells.loc[cells["sample"].str.contains("CD45pos")]#[["patient", "sample"]].drop_duplicates()
cd45neg_samples = cells.loc[cells["sample"].str.contains("CD45neg")]#[["patient", "sample"]].drop_duplicates()

tcell_samples = cd45pos_samples.loc[cd45pos_samples["celltype1"] == "Tcell"][["patient", "sample"]].drop_duplicates()
bcell_samples = cd45pos_samples.loc[cd45pos_samples["celltype1"] == "Bcell"][["patient", "sample"]].drop_duplicates()
myeloid_samples = cd45pos_samples.loc[cd45pos_samples["celltype1"] == "Myeloid"][["patient", "sample"]].drop_duplicates()

tumor_samples = cd45neg_samples.loc[cd45neg_samples["Tumor"] == "Tumor"][["patient", "sample"]].drop_duplicates()

stromal_samples = cd45neg_samples.loc[cd45neg_samples["celltype1"] == "Stromal"][["patient", "sample"]].drop_duplicates()


SAMPLE_FISHER_MARKERS = join("output", "{}", "{}", "fisher-exact-{}-{}-{}-genus-PathSeq-All-0-2.tsv")

output = []

for _, row in myeloid_samples.iterrows():
    res_df = pd.read_csv(SAMPLE_FISHER_MARKERS.format(row["patient"], row["sample"], "Myeloid", "Myeloid", "nonMyeloid"), sep="\t")
    res_df["patient"] = row["patient"]
    res_df["sample"] = row["sample"]
    res_df["fdr"] = fdrcorrection(res_df["p.value"], method="i")[1]
    output.append(res_df)

myeloid_res_df = pd.concat(output)
myeloid_res_df["celltype"] = "Myeloid"
myeloid_res_df["celltype_of_interest"] = "Myeloid"
myeloid_res_df["celltype_comparison"] = "nonMyeloid"

output = []

for _, row in bcell_samples.iterrows():
    res_df = pd.read_csv(SAMPLE_FISHER_MARKERS.format(row["patient"], row["sample"], "Bcell", "Bcell", "nonBcell"), sep="\t")
    res_df["patient"] = row["patient"]
    res_df["sample"] = row["sample"]
    res_df["fdr"] = fdrcorrection(res_df["p.value"], method="i")[1]
    output.append(res_df)

Bcell_res_df = pd.concat(output)
Bcell_res_df["celltype"] = "Bcell"
Bcell_res_df["celltype_of_interest"] = "Bcell"
Bcell_res_df["celltype_comparison"] = "nonBcell"

output = []

for _, row in tcell_samples.iterrows():
    res_df = pd.read_csv(SAMPLE_FISHER_MARKERS.format(row["patient"], row["sample"], "Tcell", "Tcell", "nonTcell"), sep="\t")
    res_df["patient"] = row["patient"]
    res_df["sample"] = row["sample"]
    res_df["fdr"] = fdrcorrection(res_df["p.value"], method="i")[1]
    output.append(res_df)

TNKILC_res_df = pd.concat(output)
TNKILC_res_df["celltype"] = "Tcell"
TNKILC_res_df["celltype_of_interest"] = "Tcell"
TNKILC_res_df["celltype_comparison"] = "nonTcell"



output = []

for _, row in tumor_samples.iterrows():
    res_df = pd.read_csv(SAMPLE_FISHER_MARKERS.format(row["patient"], row["sample"], "Tumor", "Tumor", "nonTumor"), sep="\t")
    res_df["patient"] = row["patient"]
    res_df["sample"] = row["sample"]
    res_df["fdr"] = fdrcorrection(res_df["p.value"], method="i")[1]
    output.append(res_df)

tumor_res_df = pd.concat(output)
tumor_res_df["celltype"] = "Tumor"
tumor_res_df["celltype_of_interest"] = "Tumor"
tumor_res_df["celltype_comparison"] = "nonTumor"


output = []

for _, row in stromal_samples.iterrows():
    res_df = pd.read_csv(SAMPLE_FISHER_MARKERS.format(row["patient"], row["sample"], "Stromal", "Stromal", "nonStromal"), sep="\t")
    res_df["patient"] = row["patient"]
    res_df["sample"] = row["sample"]
    res_df["fdr"] = fdrcorrection(res_df["p.value"], method="i")[1]
    output.append(res_df)

stromal_res_df = pd.concat(output)
stromal_res_df["celltype"] = "Stromal"
stromal_res_df["celltype_of_interest"] = "Stromal"
stromal_res_df["celltype_comparison"] = "nonStromal"


df = pd.concat([myeloid_res_df, TNKILC_res_df, Bcell_res_df, tumor_res_df, stromal_res_df])

df.to_csv(snakemake.output[0], sep="\t", index=False)

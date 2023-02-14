import pandas as pd
from os.path import join
from statsmodels.stats.multitest import fdrcorrection


cells = pd.read_csv("data/units.tsv", sep="\t")

# let's just focus on tumor samples
cells = cells.loc[cells["sample"].str.contains("T")]
cells["selection"] = cells["sample"].apply(lambda x: x.split("-")[1])
cells = cells.loc[~cells.patient.isin(["P38", "P39", "P40"])]

cd45pos_samples = cells.loc[cells["selection"] == "CD45pos"]#[["patient", "sample"]].drop_duplicates()
cd45neg_samples = cells.loc[cells["selection"] == "CD45neg"]#[["patient", "sample"]].drop_duplicates()

tcell_samples = cd45pos_samples.loc[cd45pos_samples["celltype1"] == "Tcell"][["patient", "sample"]].drop_duplicates()
bcell_samples = cd45pos_samples.loc[cd45pos_samples["celltype1"] == "Bcell"][["patient", "sample"]].drop_duplicates()
myeloid_samples = cd45pos_samples.loc[cd45pos_samples["celltype1"] == "Myeloid"][["patient", "sample"]].drop_duplicates()

tumor_samples = cd45neg_samples.loc[cd45neg_samples["Tumor"] == "Tumor"][["patient", "sample"]].drop_duplicates()
stromal_samples = cd45neg_samples.loc[cd45neg_samples["celltype1"] == "Stromal"][["patient", "sample"]].drop_duplicates()

CELLTYPE_ENRICHMENT = join("output", "{}", "{}", "fisher-exact-{}-{}-{}-genus-PathSeq-All-0-2.tsv")

# SAMPLE_hFDR_FISHER_MARKERS = join("output", "results", "{}", "{}", "hFDR-fisher-exact-{}-{}-{}-PathSeq-All-0-1.tsv")

output = []

celltype = "Myeloid"
celltype_comparison = "nonMyeloid"
for _, row in myeloid_samples.iterrows():
    res_df = pd.read_csv(CELLTYPE_ENRICHMENT.format(row["patient"], row["sample"], celltype, celltype, celltype_comparison), sep="\t")
    res_df["patient"] = row["patient"]
    res_df["sample"] = row["sample"]
    res_df["celltype_of_interest"] = celltype
    res_df["celltype_comparison"] = celltype_comparison
    res_df["fdr"] = fdrcorrection(res_df["p.value"], method="i")[1]
    output.append(res_df)

celltype = "Tcell"
celltype_comparison = "nonTcell"
for _, row in tcell_samples.iterrows():
    res_df = pd.read_csv(CELLTYPE_ENRICHMENT.format(row["patient"], row["sample"], celltype, celltype, celltype_comparison), sep="\t")
    res_df["patient"] = row["patient"]
    res_df["sample"] = row["sample"]
    res_df["celltype_of_interest"] = celltype
    res_df["celltype_comparison"] = celltype_comparison
    res_df["fdr"] = fdrcorrection(res_df["p.value"], method="i")[1]
    output.append(res_df)


celltype = "Bcell"
celltype_comparison = "nonBcell"
for _, row in bcell_samples.iterrows():
    res_df = pd.read_csv(CELLTYPE_ENRICHMENT.format(row["patient"], row["sample"], celltype, celltype, celltype_comparison), sep="\t")
    res_df["patient"] = row["patient"]
    res_df["sample"] = row["sample"]
    res_df["celltype_of_interest"] = celltype
    res_df["celltype_comparison"] = celltype_comparison
    res_df["fdr"] = fdrcorrection(res_df["p.value"], method="i")[1]
    output.append(res_df)

celltype = "Tumor"
celltype_comparison = "nonTumor"
for _, row in bcell_samples.iterrows():
    res_df = pd.read_csv(CELLTYPE_ENRICHMENT.format(row["patient"], row["sample"], celltype, celltype, celltype_comparison), sep="\t")
    res_df["patient"] = row["patient"]
    res_df["sample"] = row["sample"]
    res_df["celltype_of_interest"] = celltype
    res_df["celltype_comparison"] = celltype_comparison
    res_df["fdr"] = fdrcorrection(res_df["p.value"], method="i")[1]
    output.append(res_df)


df.to_csv(snakemake.output[0], sep="\t", index=False)

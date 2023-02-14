import pandas as pd
from os.path import join
from statsmodels.stats.multitest import fdrcorrection


cells = pd.read_csv("data/units.tsv", sep="\t")
cells = cells.loc[cells["Is_Tumor"] == "Yes"]
cells = cells.loc[cells["10x_chemistry"] == "v3"]
df = cells
samples_with_myeloid_cells = df.loc[df["celltype1"] == "Myeloid"][["patient", "sample"]].drop_duplicates()
samples_with_tumor_cells = df.loc[df["celltype1"] == "Epithelial"][["patient", "sample"]].drop_duplicates()
samples_with_B_cells = df.loc[df["celltype1"] == "B"][["patient", "sample"]].drop_duplicates()
samples_with_TNKILC_cells = df.loc[df["celltype1"] == "TNKILC"][["patient", "sample"]].drop_duplicates()
samples_with_Strom_cells = df.loc[df["celltype1"] == "Stromal"][["patient", "sample"]].drop_duplicates()


CELLTYPE_ENRICHMENT = join("output", "{}", "{}", "fisher-exact-{}-{}-{}-genus-PathSeq-All-0-2.tsv")
# CELLTYPE_ENRICHMENT = join("output", "{}", "{}", "fisher-exact-Tumor-Tumor-nonTumor-genus-PathSeq-All-0-2.tsv")

output = []

celltype = "Myeloid"
celltype_comparison = "nonMyeloid"
for _, row in samples_with_myeloid_cells.iterrows():
    res_df = pd.read_csv(CELLTYPE_ENRICHMENT.format(row["patient"], row["sample"], celltype, celltype, celltype_comparison), sep="\t")
    res_df["patient"] = row["patient"]
    res_df["sample"] = row["sample"]
    res_df["celltype_of_interest"] = celltype
    res_df["celltype_comparison"] = celltype_comparison
    res_df["fdr"] = fdrcorrection(res_df["p.value"], method="i")[1]
    output.append(res_df)

celltype = "Tumor"
celltype_comparison = "nonTumor"
for _, row in samples_with_tumor_cells.iterrows():
    res_df = pd.read_csv(CELLTYPE_ENRICHMENT.format(row["patient"], row["sample"], celltype, celltype, celltype_comparison), sep="\t")
    res_df["patient"] = row["patient"]
    res_df["sample"] = row["sample"]
    res_df["celltype_of_interest"] = celltype
    res_df["celltype_comparison"] = celltype_comparison
    res_df["fdr"] = fdrcorrection(res_df["p.value"], method="i")[1]
    output.append(res_df)

celltype = "TNKILC"
celltype_comparison = "nonTNKILC"
for _, row in samples_with_TNKILC_cells.iterrows():
    res_df = pd.read_csv(CELLTYPE_ENRICHMENT.format(row["patient"], row["sample"], celltype, celltype, celltype_comparison), sep="\t")
    res_df["patient"] = row["patient"]
    res_df["sample"] = row["sample"]
    res_df["celltype_of_interest"] = celltype
    res_df["celltype_comparison"] = celltype_comparison
    res_df["fdr"] = fdrcorrection(res_df["p.value"], method="i")[1]
    output.append(res_df)

celltype = "Stromal"
celltype_comparison = "nonStromal"
for _, row in samples_with_Strom_cells.iterrows():
    res_df = pd.read_csv(CELLTYPE_ENRICHMENT.format(row["patient"], row["sample"], celltype, celltype, celltype_comparison), sep="\t")
    res_df["patient"] = row["patient"]
    res_df["sample"] = row["sample"]
    res_df["celltype_of_interest"] = celltype
    res_df["celltype_comparison"] = celltype_comparison
    res_df["fdr"] = fdrcorrection(res_df["p.value"], method="i")[1]
    output.append(res_df)

celltype = "B"
celltype_comparison = "nonB"
for _, row in samples_with_B_cells.iterrows():
    res_df = pd.read_csv(CELLTYPE_ENRICHMENT.format(row["patient"], row["sample"], celltype, celltype, celltype_comparison), sep="\t")
    res_df["patient"] = row["patient"]
    res_df["sample"] = row["sample"]
    res_df["celltype_of_interest"] = celltype
    res_df["celltype_comparison"] = celltype_comparison
    res_df["fdr"] = fdrcorrection(res_df["p.value"], method="i")[1]
    output.append(res_df)

df = pd.concat(output)

df = df.loc[(df.fdr < .05) & (df["summary.odds.ratio"] > 1) & (df["n.positive.cells.of.interest"] > 1)]

df.to_csv("output/table_s3_sample_genera_enrichments.tsv", sep="\t", index=False)

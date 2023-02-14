import pandas as pd
from os.path import join


cells = pd.read_csv("data/units.tsv", sep="\t")
cells = cells.loc[cells["Is_Tumor"] == "Yes"]
cells = cells.loc[cells["10x_chemistry"] == "v3"]
samples = cells[["patient", "sample"]].drop_duplicates()


CELLTYPE_ENRICHMENT = join("output", "{}", "{}", "binomial-celltype1-genus-PathSeq-All-any-up-2.tsv")
# CELLTYPE_ENRICHMENT = join("output", "{}", "{}", "fisher-exact-Tumor-Tumor-nonTumor-genus-PathSeq-All-0-2.tsv")

output = []

for _, row in samples.iterrows():
    res_df = pd.read_csv(CELLTYPE_ENRICHMENT.format(row["patient"], row["sample"]), sep="\t")
    res_df["patient"] = row["patient"]
    res_df["sample"] = row["sample"]
    output.append(res_df)

df = df.loc[(df.FDR < .05)]

df.to_csv("output/table_s3_sample_genera_enrichments.tsv", sep="\t", index=False)

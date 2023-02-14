import pandas as pd
from os.path import join


cells = pd.read_csv("data/units.tsv", sep="\t")

cells = cells.loc[cells["sample"].str.contains("T")]
cells = cells.loc[~cells.patient.isin(["P38", "P39", "P40"])]

samples = cells[["patient", "sample"]].drop_duplicates()


CELLTYPE_ENRICHMENT = join("output", "{}", "{}", "binomial-celltype1-genus-PathSeq-All-any-up-2.tsv")

output = []

for _, row in samples.iterrows():
    res_df = pd.read_csv(CELLTYPE_ENRICHMENT.format(row["patient"], row["sample"]), sep="\t")
    res_df["patient"] = row["patient"]
    res_df["sample"] = row["sample"]
    output.append(res_df)

df = pd.concat(output)

df = df.loc[(df.FDR < .05)]

df.to_csv(snakemake.output[0], sep="\t", index=False)

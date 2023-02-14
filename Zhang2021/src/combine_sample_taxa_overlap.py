import pandas as pd
from os.path import join


cells = pd.read_csv("data/units.tsv", sep="\t")


samples = cells[["patient", "sample"]].drop_duplicates()


TAXA_OVERLAP = join("output", "taxa_overlap", "{}_{}_genus_PathSeq_All_2_enrichment.tsv")

output = []

for _, row in samples.iterrows():
    res_df = pd.read_csv(TAXA_OVERLAP.format(row["patient"], row["sample"]), sep="\t")
    res_df["patient"] = row["patient"]
    res_df["sample"] = row["sample"]
    output.append(res_df)

df = pd.concat(output)

df.to_csv(snakemake.output[0], sep="\t", index=False)

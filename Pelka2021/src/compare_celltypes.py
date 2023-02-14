import pandas as pd


df = pd.read_csv("output/all/class_PathSeq_All_reads.tsv", sep="\t")

units = pd.read_csv("data/units.tsv", sep="\t")
patients = pd.read_csv("data/patients.tsv", sep="\t")
samples = pd.read_csv("data/samples.tsv", sep="\t")

units = units.merge(patients[["patient", "MMRStatus"]], on="patient")
units = samples[["sample", "Is_Tumor"]].drop_duplicates().merge(units, on="sample")
units.index = units.apply(lambda x: "{}-{}".format(x["sample"], x["barcode"]), axis=1)
units["10x_chemistry"] = units.apply(lambda x: x["sample"].split("_")[-1], axis=1)

df = df.set_index("name")
df = df.T
bacteria = df.columns
species_df = df
df = df.merge(units, left_index=True, right_index=True)

from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import hypergeom

min_umis = 1

for patient in class_df["patient"].unique():
    df = class_df.loc[class_df.patient == patient]
output = []
for b1 in top_bacteria:
    for b2 in top_bacteria:
        if b1 != b2:
            b1_n = df.loc[df[b1] > min_umis].shape[0]
            b2_n = df.loc[df[b2] > min_umis].shape[0]
            overlap_n = df.loc[(df[b1] > min_umis) & (df[b2] > min_umis)].shape[0]
            l = [b1, b2]
            l.sort()
            #print("co-occurrence p-value for {} and {}".format(b1, b2))
            pval = hypergeom.sf(overlap_n-1, df.shape[0], b1_n, b2_n)
            output.append(pd.Series(data={"b1": l[0], "b2": l[1], "pval": pval}))

overlap_df = pd.concat(output, axis=1).T
overlap_df = overlap_df.drop_duplicates(subset=["b1", "b2"])
overlap_df["fdr"] = fdrcorrection(overlap_df["pval"], method="i")[1]
overlap_df = overlap_df.sort_values(by="pval")
print(overlap_df)
overlap_df.to_csv("output/species_overlap_table.tsv", sep="\t", index=False)

min_umis = 0

output = []
for b in bacteria:
    n_cells = df.shape[0]
    n_pos_cells = df.loc[df[b] > min_umis].shape[0]
    TME_cells = df.loc[df["Is_Tumor"] == "No"]
    n_TME_cells = TME_cells.shape[0]
    n_pos_TME_cells = TME_cells.loc[TME_cells[b] > min_umis].shape[0]
    pval = hypergeom.sf(n_pos_TME_cells-1, n_cells, n_pos_cells, n_TME_cells)
    output.append(pd.Series(data={"pval": pval, "bacteria": b}))

all_output_df = pd.concat(output, axis=1).T
all_output_df["fdr"] = fdrcorrection(all_output_df["pval"], method="i")[1]
all_output_df = all_output_df.sort_values(by="pval")
all_output_df.to_csv("output/species_TME_enriched_table.tsv", sep="\t", index=False)

import pandas as pd
from numpy import log2

def get_expected_infected_cells(df, celltype, microbe, cutoff):
    # calculate the expected number of infected cells per cell-type
    n_expected_df = df.groupby(["patient", "sample"]).apply(lambda x: x.loc[(x[microbe] >= cutoff)].shape[0]*(x[celltype].value_counts()/x.shape[0]))
    n_expected_df = n_expected_df.reset_index()
    n_expected_df = n_expected_df.rename(columns={"count": "n_expected_infected"})
    return n_expected_df.groupby([celltype])["n_expected_infected"].sum()


def get_n_infected_cells(df, celltype, microbe, cutoff):
    # calculate the actual number of infected cells per cell-type
    n_infection_df = df.groupby(["patient", "sample", celltype]).apply(lambda x: x.loc[(x[microbe] >= cutoff)].shape[0])
    return n_infection_df.groupby(celltype).sum()


# from os.path import join
# units = pd.read_csv("data/units.tsv", sep="\t")
# units = units.loc[units["Is_Tumor"] == "Yes"]
# df_v3 = units.loc[units["10x_chemistry"] == "v3"]
# sample_df_v3 = df_v3[["patient", "sample"]].drop_duplicates()
# SAMPLE_MICROBE_READ_TABLE = join("output", "{patient}", "{sample}", "{tax_level}_{method}_{kingdom}_reads.tsv")
# files = [SAMPLE_MICROBE_READ_TABLE.format(patient=row.patient, sample=row["sample"], tax_level="genus", method="PathSeq", kingdom="All") for _, row in sample_df_v3.iterrows()]
files = snakemake.input["microbe_reads"]
output = []
for f in files:
    df = pd.read_csv(f, sep="\t", index_col=0)
    output.append(df)

# concat merges on the index by default
read_df = pd.concat(output, join="outer", axis=1).fillna(0)

# PATIENT_PATHSEQ_TAXID_MAP = join("output", "{patient}", "tax_id_map_{kingdom}_{method}.tsv")
# patient_df_v3 = df_v3[["patient"]].drop_duplicates()
# tax_files = [PATIENT_PATHSEQ_TAXID_MAP.format(patient=row.patient, method="PathSeq", kingdom="All") for _, row in patient_df_v3.iterrows()]
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
# drop any genera without at least 2 UMIs in at least one cell
read_df = read_df[read_df.columns[(read_df >= 2).any(axis=0)]]
genera = read_df.columns


meta_df = pd.read_csv("data/units.tsv", sep="\t")
meta_df = meta_df.loc[(meta_df["10x_chemistry"] == "v3") & (meta_df["Is_Tumor"] == "Yes")]
meta_df.index = meta_df.apply(lambda x: "{}-{}".format(x["sample"], x["barcode"]), axis=1)
df = meta_df[["patient", "sample", "barcode", "celltype1", "celltype2"]].merge(read_df, left_index=True, right_index=True)

# drop all the genera columns now
# df = df[["patient", "sample", "barcode", "celltype1", "celltype2", "infection"]]
cutoff = int(snakemake.wildcards["cutoff"])


output = []

for genus in genera:
    expected_infected_cells_per_celltype = get_expected_infected_cells(df, "celltype1", genus, cutoff)
    infected_cells_per_celltype = get_n_infected_cells(df, "celltype1", genus, cutoff)
    # print(genera)
    # print(log2(infected_cells_per_celltype/expected_infected_cells_per_celltype))
    d = {"microbe": genus,
         "expected_infected_cells": expected_infected_cells_per_celltype,
         "actual_infected_cells": infected_cells_per_celltype,
         "log2FC": log2(infected_cells_per_celltype/expected_infected_cells_per_celltype)
        }
    output.append(pd.DataFrame(data=d))

output_df = pd.concat(output).reset_index()
output_df.sort_values(by="log2FC")

# output_df.to_csv("output/n_cells_infected_by_genera_by_celltype_log2FC.tsv", sep="\t", index=False)
output_df.to_csv(snakemake.output[0], sep="\t", index=False)

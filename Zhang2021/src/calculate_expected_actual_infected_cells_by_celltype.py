import pandas as pd
from numpy import log2
from scipy.stats import fisher_exact, combine_pvalues


def get_expected_infected_cells(df, celltype_col):
    # calculate the expected number of infected cells per cell-type
    n_expected_df = df.groupby(["patient", "sample"]).apply(lambda x: x.loc[(x["infection"] == "infected")].shape[0]*(x[celltype_col].value_counts()/x.shape[0]))
    n_expected_df = n_expected_df.reset_index()
    n_expected_df = n_expected_df.rename(columns={"count": "n_expected_infected"})
    return n_expected_df
    # return n_expected_df.groupby(celltype_col)["n_expected_infected"].sum()


def get_n_infected_cells(df, celltype_col):
    # calculate the actual number of infected cells per cell-type
    n_infection_df = df.groupby(["patient", "sample", celltype_col]).apply(lambda x: x.loc[(x["infection"] == "infected")].shape[0])
    return n_infection_df
    # return n_infection_df.groupby(celltype_col).sum()

# first step is to calculate a p-value for each sample for each cell-type
# we were previously calculating this for each genera in each cell-type
# let's calculate this using the fisher.test
# step 1 - calculate a p-value for a given cell-type in a given sample
# generate 2x2 contingency table
def single_sample_celltype_fisher_exact(df, sample, celltype_col, celltype_of_interest):
    sample_df = df.loc[df["sample"] == sample]
    celltype_of_interest_df = sample_df.loc[sample_df[celltype_col] == celltype_of_interest]
    non_celltype_of_interest_df = sample_df.loc[sample_df[celltype_col] != celltype_of_interest]
    # get the number of infected cells from celltype of interest
    a = celltype_of_interest_df.loc[celltype_of_interest_df.infection == "infected"].shape[0]
    b = celltype_of_interest_df.loc[celltype_of_interest_df.infection != "infected"].shape[0]
    c = non_celltype_of_interest_df.loc[non_celltype_of_interest_df.infection == "infected"].shape[0]
    d = non_celltype_of_interest_df.loc[non_celltype_of_interest_df.infection != "infected"].shape[0]
    res = fisher_exact([[a, b], [c, d]], alternative=snakemake.wildcards["direction"])
    return res[1]

def calculate_p_value(df, celltype_col, celltype_of_interest):
    samples = df["sample"].unique()
    # celltype_of_interest = "Myeloid"
    output = []
    for sample in samples:
        pval = single_sample_celltype_fisher_exact(df, sample, celltype_col, celltype_of_interest)
        output.append(pd.Series({"sample": sample, "celltype_of_interest": celltype_of_interest, "pvalue": pval}))
    return pd.concat(output, axis=1).T

def calculate_celltype_pvalues(df, celltype_col, celltype_of_interest):
    expected_infected_cells_per_celltype = get_expected_infected_cells(df, celltype_col)
    celltype_expected_infected_ncells = expected_infected_cells_per_celltype.loc[expected_infected_cells_per_celltype[celltype_col] == celltype_of_interest]
    # infected_cells_per_celltype = get_n_infected_cells(df, celltype_col)
    celltype_pvalues = calculate_p_value(df, celltype_col, celltype_of_interest)
    celltype_pvalue_df = celltype_expected_infected_ncells[["sample", "n_expected_infected"]].merge(celltype_pvalues[["sample", "pvalue"]], on="sample")
    # Stouffer's returns -inf if any p-value is equal to 1 so we need to set .9999 as the maximum value
    celltype_pvalue_df["pvalue"] = celltype_pvalue_df["pvalue"].apply(lambda x: min(x, .9999))
    return combine_pvalues(celltype_pvalue_df["pvalue"].astype("float").values, method="stouffer", weights=celltype_pvalue_df["n_expected_infected"].astype("float").values)[1]



# units = pd.read_csv("data/units.tsv", sep="\t")
# sample_df = units[["patient", "sample"]].drop_duplicates()
# sample_df = sample_df.loc[sample_df["sample"].str.contains("T")]
# sample_df = sample_df.loc[~sample_df.patient.isin(["P38", "P39", "P40"])]
# SAMPLE_MICROBE_READ_TABLE = join("output", "{patient}", "{sample}", "{tax_level}_{method}_{kingdom}_reads.tsv")
# files = [SAMPLE_MICROBE_READ_TABLE.format(patient=row.patient, sample=row["sample"], tax_level="genus", method="PathSeq", kingdom="All") for _, row in sample_df.iterrows()]
files = snakemake.input["microbe_reads"]
output = []
for f in files:
    df = pd.read_csv(f, sep="\t", index_col=0)
    output.append(df)

# concat merges on the index by default
read_df = pd.concat(output, join="outer", axis=1).fillna(0)

# PATIENT_PATHSEQ_TAXID_MAP = join("output", "{patient}", "tax_id_map_{kingdom}_{method}.tsv")
# patient_df = units[["patient"]].drop_duplicates()
# patient_df = patient_df.loc[~patient_df.patient.isin(["P38", "P39", "P40"])]
# tax_files = [PATIENT_PATHSEQ_TAXID_MAP.format(patient=row.patient, tax_level="genus", method="PathSeq", kingdom="All") for _, row in sample_df.iterrows()]
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
meta_df = meta_df.loc[meta_df["sample"].str.contains("T")]
meta_df.index = meta_df.apply(lambda x: "{}-{}".format(x["sample"], x["barcode"]), axis=1)
df = meta_df[["patient", "sample", "barcode", "celltype1", "celltype2"]].merge(read_df, left_index=True, right_index=True)
# microbe_of_interest = "Serratia"
# celltype_of_interest = "celltype1"

df["infection"] = "uninfected"
df.loc[(df[genera] >= int(snakemake.wildcards["cutoff"])).any(axis=1), "infection"] = "infected"

# drop all the genera columns now
df = df[["patient", "sample", "barcode", "celltype1", "celltype2", "infection"]]
df.to_csv("output/infection_phenotype.tsv", sep="\t", index=False)


expected_infected_cells_per_celltype = get_expected_infected_cells(df, "celltype1").groupby("celltype1")["n_expected_infected"].sum()
infected_cells_per_celltype = get_n_infected_cells(df, "celltype1").groupby("celltype1").sum()
output_df = infected_cells_per_celltype.to_frame(name="infected").merge(expected_infected_cells_per_celltype.to_frame(name="expected"), left_index=True, right_index=True)
output_df["log2FC"] = log2(output_df["infected"]/output_df["expected"])
output_df = output_df.reset_index()

myeloid_pvalues = calculate_celltype_pvalues(df, "celltype1", "Myeloid")
stromal_pvalues = calculate_celltype_pvalues(df, "celltype1", "Stromal")
epithelial_pvalues = calculate_celltype_pvalues(df, "celltype1", "Epithelial")
B_pvalues = calculate_celltype_pvalues(df, "celltype1", "B")
TNKILC_pvalues = calculate_celltype_pvalues(df, "celltype1", "TNKILC")
pvalues_df = pd.DataFrame({"pvalue": [myeloid_pvalues, stromal_pvalues, epithelial_pvalues, B_pvalues, TNKILC_pvalues], "celltype1": ["Myeloid", "Stromal", "Epithelial", "B", "TNKILC"]})

output_df = output_df.merge(pvalues_df, on="celltype1")


output_df.to_csv(snakemake.output[0], sep="\t", index=False)

expected_infected_cells_per_celltype = get_expected_infected_cells(df, "celltype2").groupby("celltype2")["n_expected_infected"].sum()
infected_cells_per_celltype = get_n_infected_cells(df, "celltype2").groupby("celltype2").sum()
output_df = infected_cells_per_celltype.to_frame(name="infected").merge(expected_infected_cells_per_celltype.to_frame(name="expected"), left_index=True, right_index=True)
output_df["log2FC"] = log2(output_df["infected"]/output_df["expected"])

DC_pvalues = calculate_celltype_pvalues(df, "celltype2", "DC")
macrophage_pvalues = calculate_celltype_pvalues(df, "celltype2", "Macrophage")
monocyte_pvalues = calculate_celltype_pvalues(df, "celltype2", "Monocyte")
mast_pvalues = calculate_celltype_pvalues(df, "celltype2", "Mast")
pvalues_df = pd.DataFrame({"pvalue": [DC_pvalues, macrophage_pvalues, monocyte_pvalues, mast_pvalues], "celltype2": ["DC", "Macrophage", "Monocyte", "Mast"]})

output_df = output_df.merge(pvalues_df, on="celltype2", how="left")

celltype_map = df[["celltype1", "celltype2"]].drop_duplicates()

celltype_map.merge(output_df, on="celltype2").to_csv(snakemake.output[1], sep="\t", index=False)

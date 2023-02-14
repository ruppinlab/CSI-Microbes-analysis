import pandas as pd
from numpy import log2


def get_expected_infected_cells(df, celltype_col):
    # calculate the expected number of infected cells per cell-type
    n_expected_df = df.groupby(["patient", "sample"]).apply(lambda x: x.loc[(x["infection"] == "infected")].shape[0]*(x[celltype_col].value_counts()/x.shape[0]))
    n_expected_df = n_expected_df.reset_index()
    n_expected_df = n_expected_df.rename(columns={celltype_col: "n_expected_infected", "level_2": celltype_col})
    return n_expected_df.groupby(celltype_col)["n_expected_infected"].sum()


def get_n_infected_cells(df, celltype_col):
    # calculate the actual number of infected cells per cell-type
    n_infection_df = df.groupby(["patient", "sample", celltype_col]).apply(lambda x: x.loc[(x["infection"] == "infected")].shape[0])
    return n_infection_df.groupby(celltype_col).sum()


files = snakemake.input["microbe_reads"]
output = []
for f in files:
    df = pd.read_csv(f, sep="\t", index_col=0)
    output.append(df)

# concat merges on the index by default
read_df = pd.concat(output, join="outer", axis=1).fillna(0)

# add the tax id information
output = []
tax_files = snakemake.input["tax_ids"]
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
df.loc[(df[genera] >= 1).any(axis=1), "infection"] = "bystander" # change from bystander
df.loc[(df[genera] >= 2).any(axis=1), "infection"] = "infected"

# drop all the genera columns now
df = df[["patient", "sample", "barcode", "celltype1", "celltype2", "infection"]]
df.to_csv("output/infection_phenotype.tsv", sep="\t", index=False)


expected_infected_cells_per_celltype = get_expected_infected_cells(df, "celltype1")
infected_cells_per_celltype = get_n_infected_cells(df, "celltype1")
output_df = infected_cells_per_celltype.to_frame(name="infected").merge(expected_infected_cells_per_celltype.to_frame(name="expected"), left_index=True, right_index=True)
output_df["log2FC"] = log2(output_df["infected"]/output_df["expected"])
output_df = output_df.reset_index()

output_df.to_csv(snakemake.output[0], sep="\t", index=False)

expected_infected_cells_per_celltype = get_expected_infected_cells(df, "celltype2")
infected_cells_per_celltype = get_n_infected_cells(df, "celltype2")
output_df = infected_cells_per_celltype.to_frame(name="infected").merge(expected_infected_cells_per_celltype.to_frame(name="expected"), left_index=True, right_index=True)
output_df["log2FC"] = log2(output_df["infected"]/output_df["expected"])

celltype_map = df[["celltype1", "celltype2"]].drop_duplicates()

celltype_map.merge(output_df, on="celltype2").to_csv(snakemake.output[1], sep="\t", index=False)

import pandas as pd


def get_infection_status(x):
    if (x >= 2).any():
        return "infected"
    if (x >= 1).any():
        return "bystander"
    return "uninfected"

files = snakemake.input["microbe_reads"]
output = []
for f in files:
    df = pd.read_csv(f, sep="\t", index_col=0)
    output.append(df)

# concat merges on the index by default
df = pd.concat(output, join="outer", axis=1).fillna(0)

# add the tax id information
output = []
tax_files = snakemake.input["tax_ids"]
for f in tax_files:
    output.append(pd.read_csv(f, sep="\t"))

tax_map = pd.concat(output).drop_duplicates()

df = df.merge(tax_map, left_index=True, right_on="tax_id")
df = df.set_index("name")
df = df.drop(columns=["tax_id", "taxa_level"])

# calculate the number of microbial UMIs per cell
# df = df.T
print(df)
print(df.sum())
df.sum().to_frame("n_microbial_UMIs").to_csv(snakemake.output[0], sep="\t")
df = df.T
df["infection"] = df.apply(get_infection_status, axis=1)
df["infection"].to_csv(snakemake.output[1], sep="\t")

# df.to_csv(snakemake.output[0], sep="\t")

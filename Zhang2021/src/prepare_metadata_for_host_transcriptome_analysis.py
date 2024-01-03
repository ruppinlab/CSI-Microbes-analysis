import pandas as pd


df = pd.read_csv("data/units.tsv", sep="\t")

df = df.loc[df["sample"].str.contains("T")]
# df = df.loc[~df.patient.isin(["P38", "P39", "P40"])]

tax_level = "genus"
read_file = "output/{}/{}_PathSeq_All_reads.tsv"
output = []
for p in df.patient.unique():
    output.append(pd.read_csv(read_file.format(p, tax_level), sep="\t", index_col=0))

read_df = pd.concat(output, join="outer", axis=1).fillna(0)
read_df = read_df.T

meta_df = pd.read_csv("data/units.tsv", sep="\t")
meta_df.index = meta_df.apply(lambda x: "{}-{}".format(x["sample"], x["barcode"]), axis=1)

tax_map_file = "output/{}/tax_id_map_All_PathSeq.tsv"
output = []
for p in df.patient.unique():
    df = pd.read_csv(tax_map_file.format(p), sep="\t", index_col=0)
    output.append(df)

tax_df = pd.concat(output).drop_duplicates()
tax_df = tax_df.loc[tax_df["taxa_level"] == tax_level]
d = dict(zip(tax_df["tax_id"], tax_df.index))
read_df = read_df.rename(columns=d)
df = meta_df.merge(read_df, left_index=True, right_index=True)
df = df.loc[df["sample"].str.contains("T")]
# calculate the total number of microbial UMIs

df["n_microbial_UMIs"] = df[read_df.columns].sum(axis=1)
# binarize infection
df["infection"] = "uninfected"
df.loc[(df[read_df.columns] >= 1).any(axis=1), "infection"] = "bystander" # change from bystander
df.loc[(df[read_df.columns] >= 2).any(axis=1), "infection"] = "infected"
# add infection for other cut-offs
umi_cutoffs = [1,2,3,4,5]

for umi_cutoff in umi_cutoffs:
    infection_column = "infection_{}".format(umi_cutoff)
    print(umi_cutoff)
    print(infection_column)
    df[infection_column] = "uninfected"
    df.loc[(df[read_df.columns] > 0).any(axis=1), infection_column] = "bystander" # change from bystander
    df.loc[(df[read_df.columns] >= umi_cutoff).any(axis=1), infection_column] = "infected"

# special handle for P38/P39/P40 infected cells vs other infected
df["flavo_infection"] = "uninfected"
df.loc[(df[read_df.columns] >= 2).any(axis=1), "flavo_infection"] = "other_infected"
df.loc[df.patient.isin(["P38", "P39", "40"]), "flavo_infection"] = "flavo_infected"

# DE for P38
df["P38_flavo_infection"] = "uninfected"
df.loc[df["patient"] == "P38", "P38_flavo_infection"] = "infected"
df.loc[df["patient"].isin(["P39", "P40"]), "P38_flavo_infection"] = "other"

# DE for P39
df["P39_flavo_infection"] = "uninfected"
df.loc[df["patient"] == "P39", "P39_flavo_infection"] = "infected"
df.loc[df["patient"].isin(["P38", "P40"]), "P39_flavo_infection"] = "other"

# DE for P40
df["P40_flavo_infection"] = "uninfected"
df.loc[df["patient"] == "P40", "P40_flavo_infection"] = "infected"
df.loc[df["patient"].isin(["P38", "P39"]), "P40_flavo_infection"] = "other"

# focus on the top bacteria
genera = ["Parabacteroides", "Pseudomonas", "Prevotella", "Streptococcus",
          "Campylobacter", "Capnocytophaga", "Mycoplasma", "Bacteroides",
          "Leptotrichia", "Alloprevotella", "Sphingomonas", "Phocaeicola",
          "Corynebacterium", "Treponema", "Fusobacterium", "Shigella"]

df = df[["patient", "sample", "celltype1", 'celltype2', "infection", "flavo_infection", "P38_flavo_infection", "P39_flavo_infection", "P40_flavo_infection", 'infection_1', 'infection_2', 'infection_3', 'infection_4', 'infection_5', 'n_microbial_UMIs'] + list(genera)]

# binarize infection for the top genera
for genus in genera:
    infection = "{}_infection".format(genus)
    df[infection] = "uninfected"
    df.loc[(df[genera] >= 1).any(axis=1), infection] = "other"
    df.loc[(df[genus] >= 2), infection] = "infected"


df.to_csv("output/metadata_for_host_transcriptome_analysis.tsv", sep="\t")

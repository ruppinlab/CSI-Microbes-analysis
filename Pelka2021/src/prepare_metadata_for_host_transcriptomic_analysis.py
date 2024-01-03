import pandas as pd

meta_df = pd.read_csv("data/units.tsv", sep="\t")
# only focus on tumors samples and v3 chemistries
meta_df = meta_df.loc[meta_df["Is_Tumor"] == "Yes"]
meta_df = meta_df.loc[meta_df["10x_chemistry"] == "v3"]

tax_level = "genus"
read_file = "output/{}/{}_PathSeq_All_reads.tsv"
output = []
for p in meta_df.patient.unique():
    output.append(pd.read_csv(read_file.format(p, tax_level), sep="\t", index_col=0))

read_df = pd.concat(output, join="outer", axis=1).fillna(0)
read_df = read_df.T

meta_df.index = meta_df.apply(lambda x: "{}-{}".format(x["sample"], x["barcode"]), axis=1)

tax_map_file = "output/{}/tax_id_map_All_PathSeq.tsv"
output = []
for p in meta_df.patient.unique():
    output.append(pd.read_csv(tax_map_file.format(p), sep="\t", index_col=0))

tax_df = pd.concat(output).drop_duplicates()
tax_df = tax_df.loc[tax_df["taxa_level"] == tax_level]
d = dict(zip(tax_df["tax_id"], tax_df.index))
read_df = read_df.rename(columns=d)
df = meta_df.merge(read_df, left_index=True, right_index=True)

# focus only on cells sequenced with v3 from tumor samples
df = df.loc[df["10x_chemistry"] == "v3"]
df = df.loc[df["Is_Tumor"] == "Yes"]

# calculate the total number of microbial UMIs per cell
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

# focus on the top bacteria
genera = ["Parabacteroides", "Pseudomonas", "Prevotella", "Streptococcus",
          "Campylobacter", "Capnocytophaga", "Mycoplasma", "Bacteroides",
          "Leptotrichia", "Alloprevotella", "Sphingomonas", "Phocaeicola",
          "Corynebacterium", "Treponema", "Fusobacterium", "Shigella"]

df = df[["patient", "sample", "celltype1", 'celltype2', 'infection', 'infection_1', 'infection_2', 'infection_3', 'infection_4', 'infection_5', 'n_microbial_UMIs'] + list(genera)]

# change index to be identical to the h5 file
# C154_T_1_1_0_c1_v3-TGCAGATTCACATCAG-1 vs. C103_T_1_1_0_c1_v2_id-AAACCTGCATGCTAGT
df.index = df.index.map(lambda x: "{}_id-{}".format(x.split("-")[0], x.split("-")[1]))


for genus in genera:
    infection = "{}_infection".format(genus)
    df[infection] = "uninfected"
    df.loc[(df[genera] >= 1).any(axis=1), infection] = "other"
    df.loc[(df[genus] >= 2), infection] = "infected"

df.to_csv("output/metadata_for_host_transcriptome_analysis.tsv", sep="\t")
# df.loc[df.celltype1 == "Epi"].to_csv("output/tumor_cells_metadata_for_host_transcriptome_analysis.tsv", sep="\t")

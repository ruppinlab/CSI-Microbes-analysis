import pandas as pd


df = pd.read_csv("data/units.tsv", sep="\t")
# focus only on tumor samples
df = df.loc[df["sample"].str.contains("T")]

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

genera_of_focus = ['Alloprevotella', 'Bacteroides', 'Campylobacter', 'Capnocytophaga',
                   'Corynebacterium', 'Fusobacterium', 'Leptotrichia', 'Mycoplasma',
                   'Parabacteroides', 'Phocaeicola', 'Prevotella', 'Pseudomonas',
                   'Serratia', 'Sphingomonas', 'Streptococcus', 'Treponema']
genera_of_focus = list(set(genera_of_focus).intersection(read_df.columns))

read_df = read_df[genera_of_focus]

df = meta_df.merge(read_df, left_index=True, right_index=True)

df = df.loc[df["sample"].str.contains("T")]

df.index = df.index.str.replace("-", ".")

df[read_df.columns].to_csv("output/cohort_genera_matrix.tsv", sep="\t")
df[meta_df.columns].to_csv("output/cohort_metadata.tsv", sep="\t")

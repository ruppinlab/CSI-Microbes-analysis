import pandas as pd


df = pd.read_csv("data/units.tsv", sep="\t")
# focus only on tumor samples
df = df.loc[df["sample"].str.contains("T")]
df = df.loc[~df.patient.isin(["P38", "P39", "P40"])]

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
read_df = read_df.loc[df.index]


n_infected_cells = (read_df >= 2).sum()
n_infected_cells = n_infected_cells.loc[n_infected_cells > 0]
n_umis = read_df[n_infected_cells.index].sum()
n_infected_cells = n_infected_cells.to_frame("n_infected_cells")
n_patients = [df.loc[df[x] >= 2].patient.nunique() for x in n_infected_cells.index]
n_infected_cells["n_patient"] = n_patients
table_s2_df = n_infected_cells.merge(n_umis.to_frame("n_umis"), left_index=True, right_index=True)
table_s2_df.index.name = "Genera"
table_s2_df.to_csv("output/table_s2.tsv", sep="\t")

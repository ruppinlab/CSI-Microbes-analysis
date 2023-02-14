import pandas as pd


# generate matrix of all cells by all genera
df = pd.read_csv("data/units.tsv", sep="\t")
# df = df.loc[df["Is_Tumor"] == "Yes"]

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
# df = meta_df.merge(read_df, left_index=True, right_index=True)

# filter so we just focus on genera with > 100 UMIs
read_df = read_df[read_df.columns[read_df.sum() > 100]]

# read in human reads from the cd45 negative cells
meta_df = pd.read_csv("data/units.tsv", sep="\t")
cd45neg_human_read_df = pd.read_csv("data/GSE160269_CD45neg_UMIs.txt", sep=" ", index_col=0)
cd45neg_human_read_df.columns = cd45neg_human_read_df.columns.map(lambda x: x.replace("E", "CD45neg")+"-1")
cd45pos_human_read_df = pd.read_csv("data/GSE160269_CD45pos_UMIs.txt", sep=" ", index_col=0)
cd45pos_human_read_df.columns = cd45pos_human_read_df.columns.map(lambda x: x.replace("I", "CD45pos")+"-1")

overlap_cells = list(set(cd45neg_human_read_df.columns).intersection(read_df.columns))

cd45neg_human_microbe_read_df = pd.concat([read_df[overlap_cells], cd45neg_human_read_df[overlap_cells]])
cd45_neg_metadata = meta_df.loc[meta_df.cell.isin(overlap_cells)]

cd45neg_human_microbe_read_df.to_csv("output/cd45neg_human_microbe_read_df.tsv", sep="\t")
cd45_neg_metadata.to_csv("output/cd45neg_human_microbe_read_df.tsv", sep="\t")
# this took forever - let's just focus on tumor cells for now

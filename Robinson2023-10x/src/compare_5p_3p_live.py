import pandas as pd

meta_df = pd.read_csv("data/units.tsv", sep="\t")

tax_level = "genus"
read_file = "output/{}/{}_PathSeq_All_reads.tsv"
read_df = pd.read_csv(read_file.format("P1", tax_level), sep="\t", index_col=0)
read_df = read_df.T

meta_df.index = meta_df.apply(lambda x: "{}-{}".format(x["sample"], x["barcode"]), axis=1)

tax_map_file = "output/{}/tax_id_map_All_PathSeq.tsv"
tax_df = pd.read_csv(tax_map_file.format("P1"), sep="\t", index_col=0)

tax_df = tax_df.loc[tax_df["taxa_level"] == tax_level]
d = dict(zip(tax_df["tax_id"], tax_df.index))
read_df = read_df.rename(columns=d)
df = meta_df.merge(read_df, left_index=True, right_index=True)
df = df.loc[df["sample"] == "SCAF2963_3_Live"]




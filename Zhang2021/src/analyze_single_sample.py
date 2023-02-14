import pandas as pd



read_df = pd.read_csv("output/P40/genus_PathSeq_All_reads.tsv", sep="\t", index_col=0)
tax_df = pd.read_csv("output/P40/tax_id_map_All_PathSeq.tsv", sep="\t")
df = read_df.merge(tax_df[["name", "tax_id"]], left_index=True, right_on="tax_id")
df = df.set_index("name")
df = df.drop(columns="tax_id")
df = df.T
meta_df = pd.read_csv("data/units.tsv", sep="\t")
meta_df.index = meta_df.apply(lambda x: "{}-{}".format(x["sample"], x["barcode"]), axis=1)
df = df.merge(meta_df, left_index=True, right_index=True)
df.groupby(["sample"]).apply(lambda x: x.loc[x["Flavobacterium"] >= 2].shape[0]/x.shape[0])

df.sum(axis=1).sort_values(ascending=False)

(df > 2).sum(axis=1).sort_values(ascending=False)

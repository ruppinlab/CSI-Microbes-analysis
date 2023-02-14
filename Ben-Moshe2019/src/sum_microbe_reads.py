import pandas as pd

read_df = pd.read_csv("output/Pt0/genus_PathSeq_All_reads.tsv", sep="\t", index_col=0)
tax_map = pd.read_csv("output/Pt0/tax_id_map_All_PathSeq.tsv", sep="\t")

read_df = tax_map[["name", "tax_id"]].merge(read_df, left_on="tax_id", right_index=True)
read_df = read_df.set_index("name")
read_df = read_df.drop(columns="tax_id")
meta_df = pd.read_csv("output/Pt0/PathSeq_metadata.tsv", sep="\t")

meta_df = meta_df.loc[meta_df["sample"] == "GSM3454529"]

read_df = read_df[meta_df.cell]

summed_df = read_df.sum(axis=1).sort_values(ascending=False)
# take the top 9 rows
top_df = summed_df.iloc[0:9]
top_df.loc["Other"] = summed_df.iloc[9:].sum()


top_df = top_df.to_frame(name="reads")
top_df = top_df.reset_index()
top_df.to_csv("output/Pt0/summed_microbe_reads.tsv", sep="\t")

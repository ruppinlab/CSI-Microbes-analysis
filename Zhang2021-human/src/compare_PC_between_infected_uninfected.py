import pandas as pd
from scipy.stats import ranksums

celltype = "B"

meta_df = pd.read_csv("data/metadata_for_host_transcriptome_analysis.tsv", sep="\t", index_col=0)
pc_df = pd.read_csv("output/{}_PC_values.tsv".format(celltype), sep="\t")
pc_df.index = pc_df.index.str.replace(".", "-")

df = meta_df.merge(pc_df, left_index=True, right_index=True)

for i in range(1, 15):
    PC = "PC_{}".format(i)
    print(PC)
    print(ranksums(df.loc[df.infection == "infected", PC], df.loc[df.infection == "uninfected", PC]))


import pandas as pd
from scipy.stats import spearmanr


kraken_df = pd.read_csv("output/genus_Kraken_microbe_reads.tsv", sep="\t", index_col=0)
pathseq_df = pd.read_csv("output/genus_PathSeq_microbe_reads.tsv", sep="\t", index_col=0)
kraken_df = kraken_df[pathseq_df.columns]
# sum acrosss samples and see if top microbes are similar
pathseq_df = pathseq_df.sort_index(axis=1)
kraken_df = kraken_df.sort_index(axis=1)
kraken_df.index = kraken_df.index.map(lambda x: x.split(";")[-1].replace("g__", ""))


# pathseq_df = pathseq_df.loc[(pathseq_df > 2).astype(int).sum(axis=1) > 2]
# kraken_df = kraken_df.loc[(kraken_df > 2).astype(int).sum(axis=1) > 2]

top_pathseq_genera = pathseq_df.sum(axis=1).sort_values(ascending=False)[0:20]
top_kraken_genera = kraken_df.sum(axis=1).sort_values(ascending=False)[0:20]
top_genera = set(top_kraken_genera.index).intersection(top_pathseq_genera.index)
rho, pval = spearmanr(kraken_df.loc[top_genera].transpose(), pathseq_df.loc[top_genera].transpose())

pathseq_indices = ["{} (PathSeq)".format(i) for i in pathseq_df.loc[top_genera].index]
kraken_indices = ["{} (Kraken)".format(i) for i in kraken_df.loc[top_genera].index]
new_df = pd.DataFrame(rho, index=kraken_indices+pathseq_indices, columns=kraken_indices+pathseq_indices)
# subset so microbes from PathSeq are rows and microbes from Kraken are columns
new_df = new_df.loc[new_df.index.str.contains("PathSeq")]
new_df = new_df[new_df.columns[new_df.columns.str.contains("Kraken")]]

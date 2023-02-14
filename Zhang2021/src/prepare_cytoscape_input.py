from math import log10

from statsmodels.stats.multitest import fdrcorrection
import pandas as pd


overlap_df = pd.read_csv("output/taxa_overlap/P83_P83T-CD45pos_genus_PathSeq_All_1_enrichment.tsv", sep="\t")
enrichment_df = pd.read_csv("output/P83/P83T-CD45pos/fisher-exact-Myeloid-Myeloid-nonMyeloid-genus-PathSeq-All-0-1.tsv", sep="\t")
enrichment_df["fdr"] = fdrcorrection(enrichment_df["p.value"], method="i")[1]

# let's first create the nodes and their attributes (size = n_infected cells)
node_df = pd.concat([overlap_df[["b1", "b1_n"]].rename(columns={"b1": "genera", "b1_n": "n_infected_cells"}), overlap_df[["b2", "b2_n"]].rename(columns={"b2": "genera", "b2_n": "n_infected_cells"})]).drop_duplicates()
node_df = node_df.merge(enrichment_df[["taxa", "fdr"]], left_on="genera", right_on="taxa").drop(columns="taxa")
node_df["celltype_enrichment"] = ["Myeloid" if x < .05 else "None" for x in node_df.fdr]
node_df[["genera", "n_infected_cells", "celltype_enrichment"]].to_csv("output/cytoscape_input/P83_P83T-CD45pos_node_attributes_min-umis_1.tsv", sep="\t", index=False)

# now, let's generate the network
network_df = overlap_df # pd.concat([overlap_df.rename(columns={"b1": "b2", "b2": "b1", "b1_n": "b2_n", "b2_n": "b1_n"}), overlap_df])

network_df = network_df.rename(columns={"b1" : "SOURCE", "b2": "TARGET"})
# convert the p-value
network_df["log10_fdr"] = network_df["fdr"].apply(lambda x: -log10(x))
network_df = network_df.loc[network_df.fdr < .05]
# need to manually add Leptotrichia
network_df = network_df[["TARGET", "SOURCE", "log10_fdr"]].append(pd.Series(data={"SOURCE": "Leptotrichia", "TARGET": "Leptotrichia", "log10_fdr": 0}), ignore_index=True)
network_df.to_csv("output/cytoscape_input/P83_P83T-CD45pos_network_min-umis_1.tsv", sep="\t")

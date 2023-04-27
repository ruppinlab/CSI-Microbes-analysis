from math import log10

from statsmodels.stats.multitest import fdrcorrection
import pandas as pd


overlap_df = pd.read_csv("output/taxa_overlap/C163_C163_T_1_1_0_c1_v3_genus_PathSeq_All_2_enrichment.tsv", sep="\t")
myeloid_enrichment_df = pd.read_csv("output/C163/C163_T_1_1_0_c1_v3/fisher-exact-Myeloid-Myeloid-nonMyeloid-genus-PathSeq-All-0-2.tsv", sep="\t")
stromal_enrichment_df = pd.read_csv("output/C163/C163_T_1_1_0_c1_v3/fisher-exact-Stromal-Stromal-nonStromal-genus-PathSeq-All-0-2.tsv", sep="\t")

genera_of_interest = ["Bacteroides", "Fusobacterium", "Parabacteroides", "Phocaeicola", "Campylobacter", "Leptotrichia"]
overlap_df = overlap_df.loc[overlap_df["b1"].isin(genera_of_interest) & overlap_df["b2"].isin(genera_of_interest)].copy()
enrichment_df = enrichment_df.loc[enrichment_df["taxa"].isin(genera_of_interest)].copy()

enrichment_df["fdr"] = fdrcorrection(enrichment_df["p.value"], method="i")[1]
overlap_df["fdr"] = fdrcorrection(overlap_df["pval"], method="i")[1]
# let's first create the nodes and their attributes (size = n_infected cells)
node_df = pd.concat([overlap_df[["b1", "b1_n"]].rename(columns={"b1": "genera", "b1_n": "n_infected_cells"}), overlap_df[["b2", "b2_n"]].rename(columns={"b2": "genera", "b2_n": "n_infected_cells"})]).drop_duplicates()
node_df = node_df.merge(enrichment_df[["taxa", "p.value", "fdr"]], left_on="genera", right_on="taxa").drop(columns="taxa")
node_df["celltype_enrichment"] = ["Myeloid" if x < .05 else "None" for x in node_df["p.value"]]
node_df[["genera", "n_infected_cells", "celltype_enrichment"]].to_csv("output/cytoscape_input/C163_node_attributes_min-umis_2.tsv", sep="\t", index=False)

# now, let's generate the network
network_df = pd.concat([overlap_df.rename(columns={"b1": "b2", "b2": "b1", "b1_n": "b2_n", "b2_n": "b1_n"}), overlap_df])

network_df = network_df.rename(columns={"b1" : "SOURCE", "b2": "TARGET"})
# convert the p-value
network_df["log10_fdr"] = network_df["fdr"].apply(lambda x: -log10(x))
network_df = network_df.loc[network_df.fdr < .05]
network_df = network_df.loc[network_df["SOURCE"] > network_df["TARGET"]]
network_df[["TARGET", "SOURCE", "log10_fdr"]].to_csv("output/cytoscape_input/C163_network_min-umis_2.tsv", sep="\t")
# create a file for the node attributes

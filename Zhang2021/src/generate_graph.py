import pandas as pd
from pyvis.network import Network


enrich_df = pd.read_csv("output/combined_genera_enrichment.tsv", sep="\t")
overlap_df = pd.read_csv("output/combined_taxa_overlap.tsv", sep="\t")

# use enrich the count the number of occurrences for each node
bacteria_occurrence = enrich_df[["taxa", "patient"]].drop_duplicates().groupby(["taxa"])["patient"].nunique()
bacteria_occurrence.name = "n_occurrences"
bacteria_occurrence = bacteria_occurrence.reset_index()
enriched_occurrence = enrich_df.loc[enrich_df.fdr < .05][["taxa", "patient"]].drop_duplicates().groupby(["taxa"])["patient"].nunique()
enriched_occurrence.name = "n_enriched_occurrences"
enriched_occurrence = enriched_occurrence.reset_index()
bacteria_nodes = enriched_occurrence.merge(bacteria_occurrence, on="taxa", how="outer").fillna(0)
# now, let's quantify how often a microbe is connected to another microbe
# in how many patients does an edge exist?
overlap_df = overlap_df.loc[overlap_df.fdr < .05]
bacteria_edges = pd.concat([overlap_df[["b1", "b2", "patient"]].drop_duplicates(), overlap_df[["b1", "b2", "patient"]].drop_duplicates().rename(columns={"b2": "b1", "b1": "b2"})])
edges_n_patients = bacteria_edges.groupby(["b1"])["patient"].nunique()
edges_n_patients = edges_n_patients.reset_index()
edges_n_patients = edges_n_patients.rename(columns={"b1": "taxa", "patient": "n_patients"})
edges_n_partners = bacteria_edges.groupby(["b1"])["b2"].nunique()
edges_n_partners = edges_n_partners.reset_index()
edges_n_partners = edges_n_partners.rename(columns={"b1": "taxa", "b2": "n_partners"})
bacteria_edges = edges_n_patients.merge(edges_n_partners)

bacteria_nodes = bacteria_nodes.merge(bacteria_edges, on="taxa", how="outer").fillna(0)
bacteria_nodes = bacteria_nodes.loc[(bacteria_nodes["n_enriched_occurrences"] > 1) | (bacteria_nodes["n_patients"] > 1) | (bacteria_nodes["n_partners"] > 1)]



# df = df.loc[df.fdr < .05]
# bacteria_nodes = pd.concat([df[["b1", "patient"]].drop_duplicates(), df[["b2", "patient"]].drop_duplicates().rename(columns={"b2": "b1"})]).groupby("b1")["patient"].nunique()
# bacteria_nodes = bacteria_nodes.reset_index()
# bacteria_nodes = bacteria_nodes.rename(columns={"b1": "bacteria", "patient": "n_occurrences"})
# bacteria_nodes["celltype_enriched"] = "No"
# bacteria_nodes.loc[bacteria_nodes.bacteria.isin(["Bacteroides", "Parabacteroides", "Campylobacter", "Fusobacterium", "Capnocytophaga", "Treponema"]), "celltype_enriched"] = "Yes"
#
bacteria_edges = overlap_df.groupby(["b1", "b2"])["patient"].nunique().reset_index()
bacteria_edges = bacteria_edges.rename(columns={"patient": "n_occurrences"})
bacteria_edges.loc[bacteria_edges["b1"].isin(bacteria_nodes.taxa) & bacteria_edges["b2"].isin(bacteria_nodes.taxa)]

net = Network()

for _, row in bacteria_nodes.iterrows():
    if row["n_enriched_occurrences"] > 0:
        color = "#dd4b39"
    else:
        color = "#162347"
    net.add_node(row["taxa"], color=color, value=row["n_occurrences"]) # "#dd4b39"


for _, row in bacteria_edges.iterrows():
    net.add_edge(row["b1"], row["b2"], width=row["n_occurrences"])

net.toggle_physics(True)
net.show("mygraph.html")

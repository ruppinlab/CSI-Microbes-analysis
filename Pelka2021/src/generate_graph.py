import pandas as pd
from pyvis.network import Network

c163_df = pd.read_csv("output/taxa_overlap/C163_C163_T_1_1_0_c1_v3_genus_PathSeq_All_1_enrichment.tsv", sep="\t")
c169_df = pd.read_csv("output/taxa_overlap/C169_C169_T_1_1_0_c1_v3_genus_PathSeq_All_1_enrichment.tsv", sep="\t")
c170_df = pd.read_csv("output/taxa_overlap/C170_C170_T_0_0_1_c1_v3_genus_PathSeq_All_1_enrichment.tsv", sep="\t")

c163_df["patient"] = "C163"
c169_df["patient"] = "C169"
c170_df["patient"] = "C170"

df = pd.concat([c163_df, c169_df, c170_df])
df = df.loc[df.fdr < .05]
bacteria_nodes = pd.concat([df[["b1", "patient"]].drop_duplicates(), df[["b2", "patient"]].drop_duplicates().rename(columns={"b2": "b1"})]).groupby("b1")["patient"].nunique()
bacteria_nodes = bacteria_nodes.reset_index()
bacteria_nodes = bacteria_nodes.rename(columns={"b1": "bacteria", "patient": "n_occurrences"})
bacteria_nodes["celltype_enriched"] = "No"
bacteria_nodes.loc[bacteria_nodes.bacteria.isin(["Bacteroides", "Parabacteroides", "Campylobacter", "Fusobacterium", "Capnocytophaga", "Treponema"]), "celltype_enriched"] = "Yes"

bacteria_edges = df.groupby(["b1", "b2"])["patient"].nunique().reset_index()
bacteria_edges = bacteria_edges.rename(columns={"patient": "n_occurrences"})

net = Network()

for _, row in bacteria_nodes.iterrows():
    if row["celltype_enriched"] == "Yes":
        color = "#dd4b39"
    else:
        color = "#162347"
    net.add_node(row["bacteria"], color=color, value=row["n_occurrences"]) # "#dd4b39"


for _, row in bacteria_edges.iterrows():
    net.add_edge(row["b1"], row["b2"], width=row["n_occurrences"])

net.toggle_physics(True)
net.show("mygraph.html")

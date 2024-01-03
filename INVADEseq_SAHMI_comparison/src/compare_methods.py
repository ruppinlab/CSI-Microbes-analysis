import pandas as pd

# sample = "SCAF2961_1_Uninfected"
sample = snakemake.wildcards["sample"]

# CSI-Microbes results
meta_df = pd.read_csv("data/PathSeq_metadata.tsv", sep="\t", index_col=0)
meta_df = meta_df.loc[meta_df["sample"] == sample]
read_df = pd.read_csv("data/genus_PathSeq_All_reads.tsv", sep="\t", index_col=0)
read_df = read_df.T

tax_map_file = "data/tax_id_map_All_PathSeq.tsv"
tax_df = pd.read_csv(tax_map_file, sep="\t", index_col=0)
tax_df = tax_df.loc[tax_df["taxa_level"] == "genus"]
d = dict(zip(tax_df["tax_id"], tax_df.index))
read_df = read_df.rename(columns=d)
df = meta_df.merge(read_df, left_index=True, right_index=True)
df = df.drop(columns=["patient", "sample", "barcode", "celltype1", "exposure", "chemistry"])
df.index = df.index.str.replace("_", "-")

# SAHMI results
sahmi_hits = pd.read_csv("SAHMI_data/{sample}.cell_line_quantile_hits.tsv".format(sample=sample), sep="\t")

sahmi_df = pd.read_csv("SAHMI_data/{sample}.counts.txt".format(sample=sample), sep="\t")
sahmi_df = sahmi_df.loc[sahmi_df["genus"].notna()]
sahmi_df = sahmi_df.loc[sahmi_df["species"].isna()]
sahmi_df = sahmi_df.pivot(columns="genus", index="barcode", values="counts").fillna(0)
sahmi_df.index = sahmi_df.index.map(lambda x: "{}-{}-1".format(sample, x))
sahmi_df.index = sahmi_df.index.str.replace("_", "-")

# invade seq results
invadeseq_df = pd.read_csv("INVADEseq_data/P1-{sample}.gex.filtered_matrix.genus.csv".format(sample=sample))
invadeseq_df.barcode = invadeseq_df.barcode.str.replace("P1-", "")
invadeseq_df.barcode = invadeseq_df.barcode.str.replace("_", "-")
invadeseq_df = invadeseq_df.set_index("barcode")


# CSI-Microbes has 2,864 barcodes
# SAHMI has 13,790 barcodes - only includes non-zero cell barcodes
# INVADE-seq has 855 barcodes - only includes non-zero cell barcodes
# for comparison, let's focus on the barcodes that we determined to be valid

sahmi_df = sahmi_df.loc[sahmi_df.index.isin(df.index)]
invadeseq_df = invadeseq_df.loc[invadeseq_df.index.isin(df.index)]

# we want to compare CSI-Microbes, INVADE-seq, SAHMI with cell line hits and SAHMI without cell line hits
genera_hits = sahmi_hits.loc[(sahmi_hits["rank"] == "G") & (sahmi_hits["rpmm"] > sahmi_hits["0.99"])].name.to_list()
# the total number of microbial reads identified
# the number of microbial reads that map to Fusobacterium
microbial_UMIs = []
microbial_UMIs.append(df.sum().sum())
microbial_UMIs.append(invadeseq_df.sum().sum())
microbial_UMIs.append(sahmi_df.sum().sum())
microbial_UMIs.append(sahmi_df[genera_hits].sum().sum())
Fuso_UMIs = []

if "Fusobacterium" in df.columns:
    Fuso_UMIs.append(df["Fusobacterium"].sum())
else:
    Fuso_UMIs.append(0)

if "Fusobacterium" in invadeseq_df.columns:
    Fuso_UMIs.append(invadeseq_df["Fusobacterium"].sum())
else:
    Fuso_UMIs.append(0)

if "Fusobacterium" in sahmi_df.columns:
    Fuso_UMIs.append(sahmi_df["Fusobacterium"].sum())
else:
    Fuso_UMIs.append(0)
if "Fusobacterium" in genera_hits:
    Fuso_UMIs.append(sahmi_df[genera_hits]["Fusobacterium"].sum())
else:
    Fuso_UMIs.append(0)
Approaches = ["CSI-Microbes", "INVADEseq", "SAHMI_raw", "SAHMI_hits"]

comp_df = pd.DataFrame({"Approaches": Approaches, "Microbial_UMIs": microbial_UMIs, "Fuso_UMIs": Fuso_UMIs})

comp_df.to_csv(snakemake.output[0], sep="\t", index=False)


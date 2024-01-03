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


# invade seq results
invadeseq_df = pd.read_csv("INVADEseq_same_bam/P1-{sample}.gex.filtered_matrix.genus.csv".format(sample=sample))
invadeseq_df.barcode = invadeseq_df.barcode.str.replace("P1-", "")
invadeseq_df.barcode = invadeseq_df.barcode.str.replace("_", "-")
invadeseq_df = invadeseq_df.set_index("barcode")


# CSI-Microbes has 2,864 barcodes
# SAHMI has 13,790 barcodes - only includes non-zero cell barcodes
# INVADE-seq has 855 barcodes - only includes non-zero cell barcodes
# for comparison, let's focus on the barcodes that we determined to be valid

invadeseq_df = invadeseq_df.loc[invadeseq_df.index.isin(df.index)]

# the total number of microbial reads identified
# the number of microbial reads that map to Fusobacterium
microbial_UMIs = []
microbial_UMIs.append(df.sum().sum())
microbial_UMIs.append(invadeseq_df.sum().sum())
Fuso_UMIs = []

if "Fusobacterium" in df.columns:
    Fuso_UMIs.append(df["Fusobacterium"].sum())
else:
    Fuso_UMIs.append(0)

if "Fusobacterium" in invadeseq_df.columns:
    Fuso_UMIs.append(invadeseq_df["Fusobacterium"].sum())
else:
    Fuso_UMIs.append(0)



Approaches = ["CSI-Microbes", "INVADEseq"]

comp_df = pd.DataFrame({"Approaches": Approaches, "Microbial_UMIs": microbial_UMIs, "Fuso_UMIs": Fuso_UMIs})

comp_df.to_csv(snakemake.output[0], sep="\t", index=False)


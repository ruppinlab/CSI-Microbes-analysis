import pandas as pd

HCT_Jurkat_3p_live = pd.read_csv(snakemake.input[0], sep="\t")
HCT_Jurkat_3p_live["cell.type.of.interest"] = "HCT116"
HCT_Jurkat_3p_live["cell.type.comparison"] = "Jurkat"
HCT_Jurkat_3p_live["sample"] = "Live (3')"
#SAMPLE_FISHER_EXACT.format(patient="P1", sample="SCAF2963_3_Live", celltype="celltype1", celltype_of_interest="HCT116", celltype_comparison="Jurkat", tax_level="genus", method="PathSeq", kingdom="All", minprop=0, min_umis=2),
HCT_Jurkat_5p_live = pd.read_csv(snakemake.input[1], sep="\t")
HCT_Jurkat_5p_live["cell.type.of.interest"] = "HCT116"
HCT_Jurkat_5p_live["cell.type.comparison"] = "Jurkat"
HCT_Jurkat_5p_live["sample"] = "Live (5')"
#SAMPLE_FISHER_EXACT.format(patient="P1", sample="SCAF2965_5_Live", celltype="celltype1", celltype_of_interest="HCT116", celltype_comparison="Jurkat", tax_level="genus", method="PathSeq", kingdom="All", minprop=0, min_umis=2),
HCT_THP1_3p_live = pd.read_csv(snakemake.input[2], sep="\t")
HCT_THP1_3p_live["cell.type.of.interest"] = "HCT116"
HCT_THP1_3p_live["cell.type.comparison"] = "THP1"
HCT_THP1_3p_live["sample"] = "Live (3')"
#SAMPLE_FISHER_EXACT.format(patient="P1", sample="SCAF2963_3_Live", celltype="celltype1", celltype_of_interest="HCT116", celltype_comparison="THP1", tax_level="genus", method="PathSeq", kingdom="All", minprop=0, min_umis=2),
HCT_THP1_5p_live = pd.read_csv(snakemake.input[3], sep="\t")
HCT_THP1_5p_live["cell.type.of.interest"] = "HCT116"
HCT_THP1_5p_live["cell.type.comparison"] = "THP1"
HCT_THP1_5p_live["sample"] = "Live (5')"
#SAMPLE_FISHER_EXACT.format(patient="P1", sample="SCAF2965_5_Live", celltype="celltype1", celltype_of_interest="HCT116", celltype_comparison="THP1", tax_level="genus", method="PathSeq", kingdom="All", minprop=0, min_umis=2),
THP1_Jurkat_3p_live = pd.read_csv(snakemake.input[4], sep="\t")
THP1_Jurkat_3p_live["cell.type.of.interest"] = "THP1"
THP1_Jurkat_3p_live["cell.type.comparison"] = "Jurkat"
THP1_Jurkat_3p_live["sample"] = "Live (3')"
#SAMPLE_FISHER_EXACT.format(patient="P1", sample="SCAF2963_3_Live", celltype="celltype1", celltype_of_interest="THP1", celltype_comparison="Jurkat", tax_level="genus", method="PathSeq", kingdom="All", minprop=0, min_umis=2),
THP1_Jurkat_5p_live = pd.read_csv(snakemake.input[5], sep="\t")
THP1_Jurkat_5p_live["cell.type.of.interest"] = "THP1"
THP1_Jurkat_5p_live["cell.type.comparison"] = "Jurkat"
THP1_Jurkat_5p_live["sample"] = "Live (5')"
#SAMPLE_FISHER_EXACT.format(patient="P1", sample="SCAF2965_5_Live", celltype="celltype1", celltype_of_interest="THP1", celltype_comparison="Jurkat", tax_level="genus", method="PathSeq", kingdom="All", minprop=0, min_umis=2),
HCT_THP1_3p_HK = pd.read_csv(snakemake.input[6], sep="\t")
HCT_THP1_3p_HK["cell.type.of.interest"] = "HCT116"
HCT_THP1_3p_HK["cell.type.comparison"] = "THP1"
HCT_THP1_3p_HK["sample"] = "Heat-killed"
#SAMPLE_FISHER_EXACT.format(patient="P1", sample="SCAF2962_2_HK", celltype="celltype1", celltype_of_interest="HCT116", celltype_comparison="THP1", tax_level="genus", method="PathSeq", kingdom="All", minprop=0, min_umis=2),
THP1_Jurkat_3p_HK = pd.read_csv(snakemake.input[7], sep="\t")
THP1_Jurkat_3p_HK["cell.type.of.interest"] = "THP1"
THP1_Jurkat_3p_HK["cell.type.comparison"] = "Jurkat"
THP1_Jurkat_3p_HK["sample"] = "Heat-killed"
#SAMPLE_FISHER_EXACT.format(patient="P1", sample="SCAF2962_2_HK", celltype="celltype1", celltype_of_interest="THP1", celltype_comparison="Jurkat", tax_level="genus", method="PathSeq", kingdom="All", minprop=0, min_umis=2),
df = pd.concat([HCT_THP1_3p_HK, THP1_Jurkat_3p_HK])
df = pd.concat([HCT_Jurkat_3p_live, HCT_Jurkat_5p_live, HCT_THP1_3p_live, HCT_THP1_5p_live, THP1_Jurkat_3p_live, THP1_Jurkat_5p_live, HCT_THP1_3p_HK, THP1_Jurkat_3p_HK])
df.to_csv(snakemake.output[0], sep="\t")


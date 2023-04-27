import pandas as pd


# infected_emptywells_df = pd.read_csv("output/P1/wilcox-Condition-Infected-ERCC-only-genus-PathSeq-spike-All-any-0.5-plate-up-0.5.tsv", sep="\t")
infected_emptywells_df = pd.read_csv(snakemake.input[0], sep="\t")
infected_emptywells_df["celltype.of.interest"] = "Exposed to live Fn"
infected_emptywells_df["celltype.comparison"] = "Empty Wells"
infected_emptywells_df = infected_emptywells_df.drop(columns="AUC.ERCC.only")

# infected_hct116_jurkat_df = pd.read_csv("output/P1/wilcox-Cell_Type_and_Condition-HCT116-infected-Jurkat-infected-genus-PathSeq-spike-All-any-0.5-plate-up-0.5.tsv", sep="\t")
infected_hct116_jurkat_df = pd.read_csv(snakemake.input[1], sep="\t")
infected_hct116_jurkat_df["celltype.of.interest"] = "HCT116 exposed to live Fn"
infected_hct116_jurkat_df["celltype.comparison"] = "Jurkat exposed to live Fn"
infected_hct116_jurkat_df = infected_hct116_jurkat_df.drop(columns="AUC.Jurkat.infected")

df = pd.concat([infected_emptywells_df, infected_hct116_jurkat_df])

df = df[["taxa", "p.value", "FDR", "summary.AUC", "n.reads", "celltype.of.interest", "celltype.comparison"]].sort_values(by="p.value")

# df.to_csv("output/table_S1.tsv", sep="\t", index=False)
df.to_csv(snakemake.output[0], sep="\t", index=False)

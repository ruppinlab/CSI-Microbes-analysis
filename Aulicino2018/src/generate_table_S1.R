
# expand(PATIENT_WILCOX_MARKERS, kingdom=["All"], method="PathSeq",
#            patient="Pt0", celltype=["infected"], lfc=["0.5"], tax_level="genus",
#            celltype_of_interest="infected", celltype_comparison=["bystander"],
#            norm="spike", pvaltype="any", block="plate", direction=["both"], minprop="0.5"),
#     expand(PATIENT_WILCOX_MARKERS, kingdom=["All"], method="PathSeq",
#            patient="Pt0", celltype=["infected"], lfc=["0.5"], tax_level="genus",
#            celltype_of_interest="infected", celltype_comparison=["uninfected"],
#            norm="spike", pvaltype="any", block="plate", direction=["both"], minprop="0.5"),
#     expand(PATIENT_WILCOX_MARKERS, kingdom=["All"], method="PathSeq",
#            patient="Pt0", celltype=["infected"], lfc=["0.5"], tax_level="genus",
#            celltype_of_interest="bystander", celltype_comparison=["uninfected"],
#            norm="spike", pvalt

# infected.bystander.df <- read.table("output/Pt0/wilcox-infected-infected-bystander-genus-PathSeq-spike-All-any-0.5-plate-both-0.5.tsv", header=TRUE)
infected.bystander.df <- read.table(snakemake@input[[1]], header=TRUE)
infected.bystander.df$celltype.of.interest <- "infected"
infected.bystander.df$celltype.comparison <- "bystander"
# infected.uninfected.df <- read.table("output/Pt0/wilcox-infected-infected-uninfected-genus-PathSeq-spike-All-any-0.5-plate-both-0.5.tsv", header=TRUE)
infected.uninfected.df <- read.table(snakemake@input[[2]], header=TRUE)
infected.uninfected.df$celltype.of.interest <- "infected"
infected.uninfected.df$celltype.comparison <- "uninfected"
# bystander.uninfected.df <- read.table("output/Pt0/wilcox-infected-bystander-uninfected-genus-PathSeq-spike-All-any-0.5-plate-both-0.5.tsv", header=TRUE)
bystander.uninfected.df <- read.table(snakemake@input[[3]], header=TRUE)
bystander.uninfected.df$celltype.of.interest <- "bystander"
bystander.uninfected.df$celltype.comparison <- "uninfected"


df <- do.call(rbind, list(infected.bystander.df, infected.uninfected.df, bystander.uninfected.df))
print(df)
df$tax_id <- row.names(df)
df <- df[, c(4,1,5,2,3,6,7)]
row.names(df) <- NULL
#df$taxa <- NULL
df <- df[order(df$p.value, -df$n.reads),]
# write.table(df, "output/table_S1.tsv", sep="\t", row.names=FALSE)
write.table(df, snakemake@output[[1]], sep="\t", row.names=FALSE)

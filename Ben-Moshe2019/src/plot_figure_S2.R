library(scater)
library(scran)
library(ggplot2)
library(ggpubr)
#library(ggforce)
library(ggsci)

pdata <- read.table(snakemake@input[[1]], sep="\t", header=TRUE)
# pdata <- read.table("output/Pt0/GSM3454529/class_PathSeq_Bacteria_metadata.tsv", sep="\t", header=TRUE)
#pdata.exposed <- pdata[pdata$sample == "GSM3454529",]
umi.count.exposed <- read.table(snakemake@input[[2]], sep="\t", header=TRUE, row.names=1)
# umi.count.exposed <- read.table("data/SRPRISM/Pt0/GSM3454529/CB-UMI-count-SL1344.tsv", sep="\t", header=TRUE, row.names=1)
exposed.df <- merge(pdata, umi.count.exposed, by.x = "barcode", by.y = "row.names", all.x = TRUE)
exposed.df$UMI[is.na(exposed.df$UMI)] <- 0
count.df <- exposed.df["UMI"]
rownames(count.df) <- exposed.df$barcode
rownames(exposed.df) <- exposed.df$barcode

#print(exposed.df[exposed.df$UMI > 0,])

sce <- SingleCellExperiment(assays = list(counts = as.matrix(t(count.df))), colData=exposed.df)
#binom <- findMarkers(sce, test="binom", groups=sce$Monocyte, assay.type="counts")
#print(binom[["Monocyte"]])
# let's create a bar plot of the % of cells with one or more Salmonella UMI
monocyte.positive <- sum(sce[,sce$Monocyte == "Monocyte"]$UMI > 0)
monocyte.total <- dim(sce[,sce$Monocyte == "Monocyte"])[[2]]
non.monocyte.positive <- sum(sce[,sce$Monocyte == "nonMonocyte"]$UMI > 0)
non.monocyte.total <- dim(sce[,sce$Monocyte == "nonMonocyte"])[[2]]
monocyte.percentage <- monocyte.positive/monocyte.total
non.monocyte.percentage <- non.monocyte.positive/non.monocyte.total

mat <- matrix(c(monocyte.positive, non.monocyte.positive, monocyte.total, non.monocyte.total), nrow=2, ncol=2, byrow=TRUE)
f.test <- fisher.test(mat)
print(f.test)
#print(monocyte.percentage)
#print(non.monocyte.percentage)
df <- data.frame(celltype=c("monocytes", "non-monocytes"), percentage=c(monocyte.percentage, non.monocyte.percentage))
#data <- data.frame(group1=c("monocytes"), group2=c("non-monocytes"), p=round(binom$Monocyte$p.value, digits=4), y.position=max(monocyte.percentage, non.monocyte.percentage)*1.1)

df$celltype <- factor(df$celltype, levels=c("non-monocytes", "monocytes"))


fig.1d <- ggplot(df, aes(celltype, percentage, fill=celltype)) +
  geom_bar(stat = "identity") +
  theme_pubr(legend="none", base_size=9) + # , base_size= 10 to set font size
  ylab("% of cells >= 1 UMI") +
  xlab(NULL) +
  ggtitle("Salmonella") +
  # scale_color_npg() +
  scale_fill_manual(values=c("#4DBBD599", "#E64B3599")) + # inspired by scale_color_npg
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))

#
# ggbarplot(df, x = "celltype", y = "percentage",
#    fill = "celltype", color = "celltype") +
#    scale_fill_manual(values=c("#e74c3c", "#3498db")) +
#    ylab("Percentage of cells") +
#    xlab("Celltype") +
#    ggtitle("Cells with UMIs mapped to Salmonella") +
#    theme_pubr(legend="none") +
#    stat_pvalue_manual(data, "binomial test p-value = {p}") +
#    scale_y_continuous(labels = scales::percent_format(accuracy = 1))

ggsave(snakemake@output[[1]], width = 2, height = 3.5, units = "in")

library(fgsea)


# read in the pathways - this is all the pathways
pathways <- gmtPathways(snakemake@input[[1]])

# read in the differential abundance
corr.res <- read.table(snakemake@input[[2]])
ranks <- corr.res$rho
names(ranks) <- row.names(corr.res)
ranks <- sort(ranks)

# calculate fGSEA results
fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize=500, eps=0)
fgseaRes <- fgseaRes[order(pval), ]
#print(class(fgseaRes))
#print(snakemake@output[[1]])
write.table(fgseaRes[, c("pathway", "pval", "padj", "log2err", "ES", "NES", "size")], sep="\t", snakemake@output[[1]])

# sig.fgseaRes <- fgseaRes[padj < .05, ]
# sig.neg.fgseaRes <- sig.fgseaRes[NES < 0, ]
# sig.pos.fgseaRes <- sig.fgseaRes[NES > 0, ]

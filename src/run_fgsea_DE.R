library(fgsea)


# read in the pathways - this is all the pathways
pathways <- gmtPathways(snakemake@input[[1]])

# read in the differential abundance
wilcox.res <- read.table(snakemake@input[[2]], header=TRUE, row.names="gene")
#print(t.test.res)
#print(log10(t.test.res[, "p.value"]))
#print(sign(t.test.res[, "summary.logFC"]))
ranks <- wilcox.res[, "summary.AUC"]
#ranks <- log10(t.test.res[, "p.value"])* sign(t.test.res[, "summary.logFC"])
names(ranks) <- row.names(wilcox.res)
ranks <- sort(ranks)

# calculate fGSEA results
fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize=500, scoreType = "pos", eps=0)
fgseaRes <- fgseaRes[order(pval), ]
print(class(fgseaRes))
print(snakemake@output[[1]])
write.table(as.matrix(fgseaRes), snakemake@output[[1]])

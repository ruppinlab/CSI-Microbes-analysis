library(fgsea)


# read in the pathways - this is all the pathways
pathways <- gmtPathways(snakemake@input[[1]])

# read in the differential abundance
res <- read.table(snakemake@input[[2]], header=TRUE, row.names="gene")
#print(t.test.res)
#print(log10(t.test.res[, "p.value"]))
#print(sign(t.test.res[, "summary.logFC"]))
#ranks <- res[, "summary.logFC"]
ranks <- log10(res[, "p.value"])* sign(res[, "summary.logFC"]) *-1
names(ranks) <- row.names(res)
ranks <- sort(ranks)

# calculate fGSEA results
fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize=500, eps=0)
fgseaRes <- fgseaRes[order(pval), ]
print(class(fgseaRes))
print(snakemake@output[[1]])
write.table(as.matrix(fgseaRes), snakemake@output[[1]])

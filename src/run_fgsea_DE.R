library(fgsea)


# read in the pathways - this is all the pathways
pathways <- gmtPathways(snakemake@input[[1]])

# read in the differential abundance
t.test.res <- read.table(snakemake@input[[2]])
ranks <- t.test.res[, "summary.logFC"]
names(ranks) <- row.names(t.test.res)
ranks <- sort(ranks)

# calculate fGSEA results
fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize=500)
fgseaRes <- fgseaRes[order(pval), ]
print(class(fgseaRes))
print(snakemake@output[[1]])
write.table(as.matrix(fgseaRes), snakemake@output[[1]])

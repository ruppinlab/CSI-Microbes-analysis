library(fgsea)


# read in the pathways - this is all the pathways
pathways <- gmtPathways(snakemake@input[[1]])

# read in the differential abundance
res <- read.table(snakemake@input[[2]])
# print(res)
#print(t.test.res)
#print(log10(t.test.res[, "p.value"]))
#print(sign(t.test.res[, "summary.logFC"]))
#ranks <- res[, "summary.AUC"]
ranks <- log10(res[, "p.value"])* apply(res, 1, function(x) ifelse(x["summary.AUC"] > .5, -1, 1))
names(ranks) <- row.names(res)
ranks <- sort(ranks)

# calculate fGSEA results
fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize=500, eps=0)
fgseaRes <- fgseaRes[order(pval), ]
print(class(fgseaRes))
print(snakemake@output[[1]])
write.table(as.matrix(fgseaRes), snakemake@output[[1]])

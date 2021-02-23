
res <- read.table(snakemake@input[[1]], sep="\t", header=TRUE)
tax.map <- read.table(snakemake@input[[2]], sep="\t", header=TRUE)

# rename rownames from tax_id to names
rownames(res) <- lapply(rownames(res), function(x) tax.map[tax.map$tax_id == x, "name"])

write.table(res, snakemake@output[[1]], sep="\t")

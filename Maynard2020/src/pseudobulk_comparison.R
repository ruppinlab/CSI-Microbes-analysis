library(scuttle)
library(scran)

# resource - https://support.bioconductor.org/p/132898/#132956

sce <- readRDS("output/all/tumor_cells.rds")

se <- aggregateAcrossCells(sce, ids = sce$patient)

se$tumor <- "tumor"

se$infected <- ifelse(se$patient %in% c("TH231", "TH236", "TH238", "TH266"), "infected", "control")

# optional filtering step
# se <- se[,se$patient %in% c("TH067", "TH171", "TH179", "TH220", "TH226", "TH231", "TH236", "TH238", "TH248")]

dea <- pseudoBulkDGE(se, label=se$tumor, condition=se$infected, design = ~ 0 + infected, coef=NULL, contrast=c(1, -1))

print(dea$tumor[order(dea$tumor$FDR),])

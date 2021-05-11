
species.df <- read.table(snakemake@input[[1]], header=TRUE)
species.df$tax_level <- "species"
genus.df <- read.table(snakemake@input[[2]], header=TRUE)
genus.df$tax_level <- "genus"
family.df <- read.table(snakemake@input[[3]], header=TRUE)
family.df$tax_level <- "family"
order.df <- read.table(snakemake@input[[4]], header=TRUE)
order.df$tax_level <- "order"
class.df <- read.table(snakemake@input[[5]], header=TRUE)
class.df$tax_level <- "class"
phylum.df <- read.table(snakemake@input[[6]], header=TRUE)
phylum.df$tax_level <- "phylum"

df <- do.call(rbind, list(species.df, genus.df, family.df, order.df, class.df, phylum.df))
#print(df)
df$tax_id <- row.names(df)
row.names(df) <- df$taxa
df$taxa <- NULL
df <- df[order(df$p.value),]
write.table(df, snakemake@output[[1]], sep="\t")

library(structSSI)


class.df <- read.table(snakemake@input[[1]])
order.df <- read.table(snakemake@input[[2]])
family.df <- read.table(snakemake@input[[3]])
genus.df <- read.table(snakemake@input[[4]])
species.df <- read.table(snakemake@input[[5]])
el.df <- read.table(snakemake@input[[6]], header=TRUE)
#print(head(el.df))
# first, identify how many classes
n.classes <- dim(class.df)[1]
# now, let's focus on the classes with FDR < .05
sig.classes <- class.df[class.df$FDR < .05,]
for (sig.class in rownames(sig.classes)){
  print(paste("analyzing", sig.class))

  new.class.df <- class.df[rownames(class.df) == sig.class, ]
  print("new class df")
  print(new.class.df)
  #print("get all nodes where sig.class is the parent")
  # get all the nodes where sig.class is the parent
  class.el <- el.df[el.df$parent == sig.class, ]
  print("class.el")
  print(class.el)
  #print(class.el)
  # get all the p-values involving the orders that are children of sig.class
  new.order.df <- order.df[rownames(order.df) %in% class.el$child, ]
  print("new order df")
  print(new.order.df)
  class.el <- class.el[class.el$child %in% rownames(new.order.df), ]
  print("class.el")
  print(class.el)
  #print(new.order.df)
  # # get all the nodes where the orders are the parent
  # order.el <- el.df[el.df$parent %in% rownames(new.order.df), ]
  # print(order.el)
  # # get all the p-values involving the families that are children of the orders
  # new.family.df <- family.df[rownames(family.df) %in% order.el$child, ]
  # print(new.family.df)
  # # get all the nodes where the families are the parent
  # family.el <- el.df[el.df$parent %in% rownames(new.family.df), ]
  # # get all the p-values involving the genera that are children of the families
  # new.genus.df <- genus.df[rownames(genus.df) %in% family.el$child, ]
  # print(new.genus.df)
  # # get all the nodes where the genera are the parent
  # genus.el <- el.df[el.df$parent %in% rownames(new.genus.df), ]
  # # get all the p-values involving the genera that are children of the families
  # new.species.df <- species.df[rownames(species.df) %in% genus.el$child, ]
  select.df <- rbind(new.class.df, new.order.df)
  print("select df")
  print(select.df)
  vec <- select.df$p.value
  print("vec")
  print(vec)
  names(vec) <- rownames(select.df)
  print("vector of p-values")
  print(vec)
  select.el <- class.el #rbind(class.el)
  mat <- as.matrix(select.el)
  rownames(mat) <- NULL
  colnames(mat) <- NULL
  print("edgelist")
  print(mat)
  if dim(mat[1] > 1){
    
  }
  out <- hFDR.adjust(vec, mat, alpha = 0.05)
  print(out)
}


library(circlize)
library(reshape2)

## REPRESENTING LINKS ##

links = read.table("links.tsv",
  sep='\t', header=TRUE, row.names=1)
# remove (to zero) upper.triangle and diag
links[upper.tri(links)]=0
diag(links)=0
mat = data.matrix(links)
# transform adj. matrix to adj. list
melted_mat = melt(mat)
melted_mat=melted_mat[melted_mat$value!=0,]
# keep count of previously placed links
keep = data.frame(n = names(links), v=0)



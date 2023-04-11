
library(circlize)
library(reshape2)

### REPRESENTING LINKS ###

## Load files necessary
# number of connections
n_par = read.table("tracks.tsv",
  sep='\t', header=TRUE)
# [-7,] drops the ending *total* row
gtab = read.table("circos_tracks.tsv",
  sep='\t', header=TRUE)[-7,]
conn = read.table("links.tsv",
  sep='\t', header=TRUE, row.names=1)
# count total amount of paralogs
n_par=n_par[, c(1,5,6)]
# remove (to zero) upper.triangle and diag
conn[upper.tri(conn)]=0
diag(conn)=0
# percentatges de cada mena de connexió
###conn = conn/sum(conn)
###for (i in c(1:length(conn)))
###  {conn[,i]=conn[,i]/sum(conn[,i])}
###conn[6]=0
# transform adj. matrix to adj. list
mat = melt(data.matrix(conn))
mat=mat[mat$value!=0,]
# transforma factors de chr a vectors de chr:
mat[,4] = levels(mat[,2])[as.integer(mat[,2])]
mat[,5] = levels(mat[,1])[as.integer(mat[,1])]

# keep count of previously placed links
keep = data.frame(n = names(conn), v=0)
links=data.frame()
for (i in c(1:length(mat[,3]))) {
  # for the first chr
  crm_from = mat[i, 4]
  start = keep[keep$n==crm_from, 2]
  end = (mat[i, 3]/n_par[n_par==crm_from,2]) * gtab[gtab ==crm_from, 3] + start
  links[i, 1] = crm_from
  links[i, 2] = start
  links[i, 3] = end
  keep[keep ==crm_from, 2] = end

  # for the second chr
  crm_to = mat[i, 5]
  start = keep[keep ==crm_to, 2]
  end = (mat[i, 3]/n_par[n_par==crm_to, 2]) * gtab[gtab ==crm_to, 3] + start
  links[i, 4] = crm_to
  links[i, 5] = start
  links[i, 6] = end
  keep[keep ==crm_to, 2] = end
}

### REPRESENTING TRACKS ###

xlim=data.frame(crm=gtab[,1] ,start=0, end=gtab$chr_len)

### PLOTTING CIRCOS ###

circos.clear()
circos.initialize(sectors=xlim$crm, xlim=xlim[,c(2,3)])
circos.track(xlim$crm, ylim=c(0,1),
  track.height=.15, #bg.border=NA,
  panel.fun = function(x,y) {
  crm = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  #circos.rect(xlim[1], 0, xlim[2], 1, col = rand_color(1))
  circos.text(mean(xlim), mean(ylim), crm, cex = 0.7, col = "black",
  facing="inside", niceFacing=TRUE)
  })
circos.genomicLink(links[,c(1:3)], links[,c(4:6)],
  col = rand_color(nrow(links), transparency = 0.5), border=NA)


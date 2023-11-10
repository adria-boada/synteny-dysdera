
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

# set the colors for the possibly tandem segments
tandem_col=colorRampPalette(c('blue', 'yellow'))(length(gtab)-3)
#highlights
crm_list=c();start_list=c();end_list=c(); col_list=c()
for (i in c(1:length(gtab[,1]))) {
  chr=gtab[i,1]; start=0
  for (j in c(3:(length(gtab)-1))) {
  start=start+gtab[i,j]
  end=start+gtab[i,j+1]
  #print(start); print(end)
  crm_list=c(crm_list,chr)
  start_list=c(start_list, start)
  end_list=c(end_list, end)
  col_list=c(col_list, tandem_col[j-2])
  }}

# matrix with pair of starts/ends
crm_highlight=matrix(c(crm_list, start_list, end_list, col_list),
  nrow=length(crm_list))
crm_highlight=data.frame(crm_highlight)

## COLOR PALETTE ##

# 8 colors més l'afegit U2. 
paleta_escuer = data.frame(
  sil_crm=c("DsilChrX", "DsilChr1", "DsilChr2", "DsilChr3",
  "DsilChr4", "DsilChr5", "DsilChr6", "DsilChrU1",
  "DsilChrU2"),
  col=c("#f659a5", "#c5090a", "#ff7f07", "#cabf0a",
  "#41a62a", "#4fa1ca", "#4f3b97", "#a29a9c",
  "#a65628"),
  cat_links=c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Scaffold_Agg",
  "na", "na", "na")
  )

  # Test amb un barplot per veure que els colors rutllen i són els adequats:
  if (FALSE) {
  barplot(c(DsilChrX=5, DsilChr1=5, DsilChr2=5, DsilChr3=5, DsilChr4=5, DsilChr5=5, DsilChr6=5, DsilChrU1=5, DsilChrU2=5), col=paleta_escuer)
  }

# set the colors in the `links` data.frame
col_list=c()
for (i in c(1:length(links[,1]))) {
  col_list=c(col_list, paleta_escuer[links[i,1]==paleta_escuer[,3], 2])
}

### PLOTTING CIRCOS ###

# guarda 'Rplot.pdf'
pdf(width=9)

# gràfica
circos.clear()
circos.initialize(sectors=xlim$crm, xlim=xlim[,c(2,3)])
circos.track(xlim$crm, ylim=c(0,1),
  track.height=.15, #bg.border=NA,
  panel.fun = function(x,y) {
  crm = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  start_highl=as.integer(crm_highlight[crm_highlight==crm,2])
  end_highl=as.integer(crm_highlight[crm_highlight==crm,3])
  col_highl=crm_highlight[crm_highlight==crm,4]
  for (i in c(1:length(start_highl))) {circos.rect(start_highl[i], 0, end_highl[i], 1, col=col_highl[i])}
  #circos.rect(xlim[1], 0, xlim[2], 1, col = rand_color(1))
  circos.text(mean(xlim), mean(ylim), crm, cex = 0.7, col = "black",
  facing="inside", niceFacing=TRUE)
  })
circos.genomicLink(links[,c(1:3)], links[,c(4:6)],
  col = col_list, border='black', h.ratio=0.4)

# tanca fitxer
dev.off()


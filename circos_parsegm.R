#! /usr/bin/env Rscript

### MANUAL OPTIONS ###

# Read supplied arguments
args = commandArgs(trailingOnly=TRUE)
# It is necessary to supply 2 arguments:
if (length(args) != 2) {
  stop(paste("Supply exactly 2 argument files:",
             "circos_parsegm.R `segments.txt` `chr_index.idx`",
             sep='\n'),
       call.=FALSE)
}
# Path to segments dataframe file
segments_file = args[1]
# Path to a tabulated index file (2 cols; chr and its length)
index_file = args[2]

library(circlize)      # circular plots
library(colourvalues)  # colouring links depending on variable

### READ DATAFRAME ###

# pairs of segments to link
df = read.table(segments_file,
  sep='\t', header=TRUE, row.names=1)
# index of sequids in the previous file, length-sorted
crm_idx = read.table(index_file,
  sep='\t', header=FALSE, col.names=c('crm', 'end'))
# add a "starting" position for chr
crm_idx$start=0
# and reorder the columns
crm_idx = crm_idx[, c('crm', 'start', 'end')]

# specify the four legend titles
leg_titles = c(
  paste('Ratio of gene number',
  '\n(source vs. destiny segments)'),
  paste('Ratio of mean gene lengths',
  '\n(source vs. destiny segments)'),
  paste('Ratio of mean gene distances',
  '\n(source vs. destiny segments)'),
  'Intrachromosomal'
)

### SUBSETTING CIRCOS LINKS ###
#
# Make sure that the given links start or end in chromosomes which
# are found in the index file. If a segment is located in a chromosome/
# scaffold which is not present in the index file, remove it from the
# dataframe.

# Make sure that the column 'sequid_to' and 'sequid_from'
# are in the column 'crm' for the 'crm_idx' data.frame
df = df[df$sequid_to %in% crm_idx$crm & df$sequid_from %in% crm_idx$crm,]

### COLOUR OF CIRCOS LINKS ###

# specify colour links in data.frame
df_cols = data.frame(
  gn=colour_values(c(1, df[, 7])),
  gl=colour_values(c(1, df[, 8])),
  gd=colour_values(c(1, df[, 9])),
  # add a bool column which tracks whether links are intra or extra chromosomal
  intra=colour_values(c(FALSE, df$sequid_from==df$sequid_to))
)
# specify min and max colour values
col_range = colour_values(c(0,1))
# specify min and max legend values
val_range = data.frame(
      min=c(1, 1, 1, FALSE),
      max=c(
        df$gn[df$gn==max(df$gn)],
        df$gL[df$gL==max(df$gL)],
        df$gD[df$gD==max(df$gD)],
        TRUE)
      )

### PLOTTING CIRCOS FIGURE ###

# # # # # #
# export plot to 'Rplot.pdf'
# increase default width to correctly display legend
pdf(width=9)
# # # # # #

# Repeat for each column of colours
for (i in c(1:4)) {
circos.clear() # good habit to clear previous options
# circos.par(): parameters to solve the overabundance of
# minor scaffold representation in *Dcat*
circos.par(cell.padding = c(0.02, 0,0.02, 0))
circos.initialize(sectors=crm_idx$crm,
                  xlim=crm_idx[, c(2,3)])
circos.track(crm_idx$crm, ylim=c(0,.25),
  track.height=.15, #bg.border=NA,
  panel.fun = function(x,y) {
    crm = CELL_META$sector.index
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    #circos.rect(xlim[1], 0, xlim[2], 1)
    circos.text(mean(xlim), mean(ylim), crm,
      cex=1.2, col='black', facing='inside', niceFacing=TRUE)
  }
)
# connect 1:3 (origin segment) to 4:6 (recipient segment)
circos.genomicLink(df[, c(1:3)], df[, c(4:6)],
  col = df_cols[-1, i], #h.ratio = .4
)
# add a legend with the min and max of misc value,
# and the colours used to represent these.
  legend('bottomleft', title=leg_titles[i], bty='n',
    legend = c(val_range$min[i], val_range$max[i]),
    fill = c(col_range[1], col_range[2]),
    )

} # The loop for creating multiple plots ends

# # # # # #
# shut down plotting device
dev.off()
# # # # # #



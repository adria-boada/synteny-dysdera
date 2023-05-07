#! /usr/bin/env Rscript

library(circlize)      # circular plots
library(colourvalues)  # colouring links depending on variable

### MANUAL OPTIONS ###

# Read supplied arguments
args = commandArgs(trailingOnly=TRUE)
# Path to segments dataframe file (should be first and only argument)
segments_file = args[1]
# manually add chromosome length
# (would be good to add fasta.idx reading capabilities)
crm_idx = data.frame(
  crm=c('Scaffold_1', 'Scaffold_2', 'Scaffold_3', 'Scaffold_4',
        'Scaffold_5_Vmt', 'Scaffold_6', 'Scaffold_7'),
  start=0,
  end=c(434793700, 210901293, 210826211, 205634564,
        158554970, 152324531, 96325250)
)

### READ DATAFRAME ###

df = read.table(segments_file,
  sep='\t', header=TRUE, row.names=1)
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



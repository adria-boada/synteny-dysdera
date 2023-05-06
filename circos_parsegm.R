
library(circlize)
library(reshape2)
library(colourvalues)

### OPTIONS ###

# select a misc. var to represent between...
# gene number col. -> 7
# gene length col. -> 8
# gene distance col. -> 9
misc_var = 9

### READ DATAFRAME ###

df = read.table("segments_linked_pairs.tsv",
  sep='\t', header=TRUE, row.names=1)
# create a separate column with selected misc var...
df$var = df[, misc_var]
if (misc_var==7) {misc_var_title=paste('Ratio of gene number',
                                       '\n(source vs. destiny segments)')}
if (misc_var==8) {misc_var_title=paste('Ratio of mean gene lengths',
                                       '\n(source vs. destiny segments)')}
if (misc_var==9) {misc_var_title=paste('Ratio of mean gene distances',
                                       '\n(source vs. destiny segments)')}

### CREATE COLOUR COLUMN CORRELATED TO MISC RATIOS ###

# create a column with colours associated to
# gene number column
col_ramp = colour_values(c(1, df$var))
df$col = col_ramp[-1]
# check the colours with a barplot
if (FALSE) {barplot(df$var, col=df$col)}
# manually add chromosome length
# (would be good to add fasta.idx reading capabilities)
crm_idx = data.frame(
  crm=c('Scaffold_1', 'Scaffold_2', 'Scaffold_3', 'Scaffold_4',
        'Scaffold_5_Vmt', 'Scaffold_6', 'Scaffold_7'),
  start=0,
  end=c(434793700, 210901293, 210826211, 205634564,
        158554970, 152324531, 96325250)
)

### PLOTTING CIRCOS FIGURE ###

# # # # # #
# export plot to 'Rplot.pdf'
# increade default width to correctly position legend
pdf(width=9)
# # # # # #

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
  col = df[, length(df)], #h.ratio = .4
)
# add a legend with the min and max of misc value,
# and the colours used to represent these.
legend('bottomleft', title=misc_var_title, bty='n',
  legend = c(1, df$var[df$var==max(df$var)]),
  fill = c(col_ramp[1], df$col[df$var==max(df$var)]))

# # # # # #
# shut down plotting device
dev.off()
# # # # # #


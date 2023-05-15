#! /usr/bin/env Rscript


### MANUAL OPTIONS and READING DATAFRAMES ###

# Read supplied arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop(paste("Supply exactly 3 argument files:",
             "circos_parsegm.R `broquil_orthologs.tsv` `sp1_chr.idx` `sp2_chr.idx`",
             "",
             "`sp1_chr.idx` must have the following format:",
             "Species-Legend  sp-brocc.tsv  colours",
             "Scaffold_1      100           #FFFFFF",
             "Scaffold_2      50            #EEEEEE",
             "etc.",
             "`sp2_chr.idx` follows the same format but does not require colours",
             sep='\n'),
       call.=FALSE)
}
# Path to orthologs.tsv dataframe file
orthologs_file = args[1]
# Path to a tabulated index files
index_file_color = args[2]
index_file_cless = args[3]

library(circlize)      # circular plots
library(colourvalues)  # colouring links depending on variable

# pairs and multigenic families of orthologs
df = read.table(orthologs_file,
  sep='\t', header=TRUE, row.names=1)
# index of sequids in the previous file, length-sorted
idx_col = read.table(index_file_color,
  sep='\t', header=FALSE)[, c(1:3)]
idx_lss = read.table(index_file_cless,
  sep='\t', header=FALSE)[, c(1:2)]
# rename idx columns
names(idx_col) = c('crm', 'end', 'col')
names(idx_lss) = c('crm', 'end')
# read both species names ('first' word in tsv header)
sp_names = c(leg_col=idx_col[1,1], leg_lss=idx_lss[1,1],
             df_col=idx_col[1,2], df_lss=idx_lss[1,2])
# check that species in idx[1,2] are found in df$Species
if (!all(sp_names[c(3,4)] %in% unique(df$Species))) {
  stop('ERROR: The species defined in idx[1,2] do not match unique(df$Species)', call.=FALSE)
}
# if there are more than two species in dataframe, warn
if (length(unique(df$Species)) > 2) {
  stop('ERROR: There are more than two species in broccoli df!', call.=FALSE)
}
# add a "starting" position for all chr
idx_col$start = 0
idx_lss$start = 0
# and reorder the columns
highlighting_crm = data.frame(col=idx_col$col[-1])  # store colours in a new var
idx_col = idx_col[, c('crm', 'start', 'end')] # removes colour column
idx_lss = idx_lss[, c('crm', 'start', 'end')]
# paste an identifier to all scaffold
# (to identify dcat scaffold_1 from dtil scaffol_1)
idx_col$crm = paste0(sp_names[3], '@', idx_col$crm)
idx_lss$crm = paste0(sp_names[4], '@', idx_lss$crm)
# the broccoli table also has to be modified...
df$Scaffold[df$Species==sp_names[3]] = paste0(sp_names[3], '@',
  df$Scaffold[df$Species==sp_names[3]])
df$Scaffold[df$Species==sp_names[4]] = paste0(sp_names[4], '@',
  df$Scaffold[df$Species==sp_names[4]])
# remove inputted header from idx_col and idx_lss
idx_col = idx_col[-1,]
idx_lss = idx_lss[-1,]
highlighting_crm$crm = idx_col$crm  # pair colours with crm
# row bind both indexes into a single data.frame
# colour had to be removed in order to have same columns
idx = rbind(idx_col, idx_lss)
# remove all orthologs that are not included in the index
# (if minor scaffolds are not in index, remove orthologs of minor scaffolds from `df`)
df = df[df$Scaffold %in% idx$crm,]


### CREATE GENOMIC LINKS and POINTS DATAFRAMES ###

# resulting df with genomic links
df_links = data.frame()
# resulting df with genomic points
df_points = data.frame()
# keep track of removed and included OGs in analyses:
og.removed = 0 ; og.included = 0
# fill previous df with input df:
for (og in unique(df$OGid)) {
  # select all members of orthogroup
  members = df[df$OGid==og,]
  ##print(members[,c('OGtype', 'Scaffold', 'GeneStart', 'GeneEnd', 'Species', 'OGid')])#DEBUG
  members_col = members[members$Species==sp_names[3],]
  members_lss = members[members$Species==sp_names[4],]

  if (nrow(members_lss)==0 | nrow(members_col)==0) {
    # it is possible to remove all genes from a family while
    # deleting rows located at minor scaffolds
    # just skip these cases...
    og.removed = og.removed+1
    print(paste0('WARNING: ', og, ' has been removed because',
    ' at least one species has all of their orthologs outside provided index sequids.'))
  } else if (nrow(members)==2) {
    og.included = og.included+1
    # two members: we can create a link between coordinates
    # maybe even 3...
    new_row = data.frame(
      sequid_col = members_col$Scaffold,
      start_col = members_col$GeneStart,
      end_col = members_col$GeneEnd,
      sequid_lss = members_lss$Scaffold,
      start_lss = members_lss$GeneStart,
      end_lss = members_lss$GeneEnd
    )
    # append new row to links df:
    df_links = rbind(df_links, new_row)
  } else if (nrow(members_lss)==nrow(members_col)) {
    og.included = og.included+1
    # if they are a many-members orthogroup, plot as genomic points:
    # try to colour each point depending on which species
    # has more genes in the selected orthogroup
    # (compare nrows of `members_lss` versus `members_col`)
      new_rows = data.frame(
        sequid = members$Scaffold,
        start = members$GeneStart,
        end = members$GeneEnd,
        col = 'yellow',
        ogid = members$OGid
      )
      #print(new_rows)#DEBUG
      df_points = rbind(df_points, new_rows)
  } else if (nrow(members_lss) >nrow(members_col)) {
    og.included = og.included+1
    # more genes in species colourless (dcat)
      new_rows = data.frame(
        sequid = c(members_col$Scaffold, members_lss$Scaffold),
        start = c(members_col$GeneStart, members_lss$GeneStart),
        end = c(members_col$GeneEnd, members_lss$GeneEnd),
        col = c(rep('red', nrow(members_col)), rep('green', nrow(members_lss))),
        ogid = c(members_col$OGid, members_lss$OGid) #DEBUG
      )
      #print(new_rows)#DEBUG
      df_points = rbind(df_points, new_rows)
  } else if (nrow(members_lss)< nrow(members_col)) {
    og.included = og.included+1
    # more genes in species colourful (dtil)
      new_rows = data.frame(
        sequid = c(members_col$Scaffold, members_lss$Scaffold),
        start = c(members_col$GeneStart, members_lss$GeneStart),
        end = c(members_col$GeneEnd, members_lss$GeneEnd),
        col = c(rep('green', nrow(members_col)), rep('red', nrow(members_lss))),
        ogid = c(members_col$OGid, members_lss$OGid) #DEBUG
      )
      #print(new_rows)#DEBUG
      df_points = rbind(df_points, new_rows)
  }
}
print(paste('OGs included:', og.included))
print(paste('OGs removed:', og.removed))
print(paste('Overall OGs:', length(unique(df$OGid))))
# add colours to df_links depending on idx_col of 'origin'
df_links$col = 0
for (i in c(1:nrow(highlighting_crm))) {
  df_links$col[df_links$sequid_col==highlighting_crm[i, 2]]=highlighting_crm[i, 1]
}
# create lines instead of genomicPoints()
# the df_points will have to be formatted differently (create df_lines)
df_lines=data.frame()
for (c in idx$crm) {
  e = idx$end[idx$crm == c]
  s = seq(0, e, 1e6)
  s = s[-length(s)]
  for (i in s) {
    # select rows with start between `i` and `i+1e6` and in scaffold `c`
    df_sel = df_points[df_points$start >= i & df_points$start < i+1e6 & df_points$sequid == c,]
    new_rows = data.frame(
      sequid = c,
      start = i,
      end = i+1e6,
      amount_y = nrow(df_sel[df_sel$col=='yellow', ]),
      amount_r = nrow(df_sel[df_sel$col=='red', ]),
      amount_g = nrow(df_sel[df_sel$col=='green', ]))
    df_lines = rbind(df_lines, new_rows)
  }
  df_sel = df_points[i+1e6 <= df_points$start & df_points$start < e & df_points$sequid == c,]
    new_rows = data.frame(
    sequid = c,
    start = i+1e6,
    end = as.numeric(e),
    amount_y = nrow(df_sel[df_sel$col=='yellow', ]),
    amount_r = nrow(df_sel[df_sel$col=='red', ]),
    amount_g = nrow(df_sel[df_sel$col=='green', ]))
df_lines = rbind(df_lines, new_rows)
}
track_ylim = max(c(df_lines$amount_y,
                   df_lines$amount_r,
                   df_lines$amount_g))
print(paste('Most genes in `df_lines.green`:', max(df_lines$amount_g)))
print(paste('Most genes in `df_lines.red`:', max(df_lines$amount_r)))
print(paste('Most genes in `df_lines.yellow`:', max(df_lines$amount_y)))
print(head(n=20, df_lines)) #DEBUG


### PLOTTING CIRCOS FIGURE ###

# # # # # #
# export plot to 'Rplot.pdf'
# increase default width to correctly display legend
pdf(width=9)
# # # # # #

# good habit to clear previous options
circos.clear()
# separate all crms by 1 and species by 5
circos.par(gap.degree = c(rep(1, nrow(idx_col)-1), 5,
  rep(1, nrow(idx_lss)-1), 5))
# ERROR: summation of xdirection is larger than the width
# of some sectors. Remove tiny sectors or reset cell.padding as follows...
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
# initialize with sectors <- chromosome names
circos.initialize(sectors=idx$crm,
  xlim=idx[, c(2,3)])
circos.track(idx$crm, ylim=c(0, track_ylim),
  track.height=.15, #bg.border=NA,
  panel.fun  = function(x,y) {
    # shorthands for a few sector properties
    cllcrm = CELL_META$sector.index
    xlim   = CELL_META$xlim
    ylim   = CELL_META$ylim
    # draw points where members of multigenic orthogroups are at
    # colour points by more/less/equal members in orthogroup for both species
    ##circos.genomicPoints(df_points[df_points$sequid==cllcrm, c('start', 'end')],
    ##  CELL_META$ycenter,  # place points in the Y-center of the cell
    ##  col=df_points[df_points$sequid==cllcrm, 'col'],
    ##  pch=16, cex=2) # make them big and round dots
    # instead of points, draw multigenic orthogroups as
    # repetitive elements would be drawn
    d = df_lines[df_lines$sequid==cllcrm, ]
    circos.genomicLines(d[, c('start', 'end')],
    d[, 'amount_g'], type='segment', col='green')
    circos.genomicLines(d[, c('start', 'end')],
    d[, 'amount_y'], type='segment', col='yellow')
    circos.genomicLines(d[, c('start', 'end')],
    d[, 'amount_r'], type='segment', col='red')
    circos.text(mean(xlim), mean(ylim)+mm_y(9), # add Y space (chr.name above track)
      gsub(".*d_", "", cllcrm),  # crm name that will be displayed
      # the global substitution (gsub) can trim the sp identifier or 'Scaffold' entirely
      cex=1.2, col='black', facing='inside', niceFacing=TRUE)
  }
)
# connect 1:3 (origin segment) to 4:6 (recipient segment)
circos.genomicLink(df_links[, c(1:3)], df_links[, c(4:6)],
  col = df_links[, 'col'], #h.ratio = .4
)
# add a legend with the min and max of misc value,
# and the colours used to represent these.
#  legend('bottomleft', title=leg_titles[i], bty='n',
#    legend = c(val_range$min[i], val_range$max[i]),
#    fill = c(col_range[1], col_range[2]),
#    )

# # # # # #
# shut down plotting device
dev.off()
# # # # # #



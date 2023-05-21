#! /usr/bin/env Rscript


### MANUAL OPTIONS and READING DATAFRAMES ###

# Read supplied files in argument line
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  # help message (if less than the needed args are supplied)
  stop(paste("Supply exactly 3 argument files:",
             "circos_parsegm.R `broquil_orthologs.tsv` `sp1_chr.idx` `sp2_chr.idx`",
             "",
             "`sp1_chr.idx` must have the following format:    ",
             "Species-Legend  sp-brocc.tsv  colours  2ndary-col",
             "Scaffold_1      100           #FFFFFF  #DDDDDD   ",
             "Scaffold_2      50            #EEEEEE  #CCCCCC   ",
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
# adjust window size for counting multigenic families per window
window_size = 4e6 # base-pairs
# remove as many as `remove_top` outliers values
# for amount_g so it is easier to see differences
remove_top=1

library(circlize)      # circular plots
library(colourvalues)  # colouring links depending on variable

# pairs and multigenic families of orthologs
df_input_table = read.table(orthologs_file,
  sep='\t', header=TRUE, row.names=1)
# index of sequids in the previous file, length-sorted
idx_col = read.table(index_file_color,
  sep='\t', header=FALSE)[, c(1:4)]
idx_lss = read.table(index_file_cless,
  sep='\t', header=FALSE)[, c(1:2)]
# rename chromosome_index_dataframe's columns
names(idx_col) = c('crm', 'end', 'col', 'col_2ary')
names(idx_lss) = c('crm', 'end')
# read both species complete names ('first' field in tsv header)
# moreover, read shortened names found in tsv ('second' field in tsv header
sp_names = c(leg_col=idx_col[1,1], leg_lss=idx_lss[1,1],
             df_col=idx_col[1,2], df_lss=idx_lss[1,2])
# check that species in idx[1,2] for both files are found in df$Species
if (!all(sp_names[c(3,4)] %in% unique(df_input_table$Species))) {
  stop('ERROR: The species defined in idx[1,2] do not match unique(df$Species)', call.=FALSE)
}
# if there are more than two species in dataframe, warn
if (length(unique(df_input_table$Species)) > 2) {
  stop('ERROR: There are more than two species in broccoli df!', call.=FALSE)
}
# add a "starting" position for all chr
idx_col$start = 0
idx_lss$start = 0
# store general colours in a new data.frame
# pair these colours with scaffolds
highlighting_crm = data.frame(
  sequid = 0,  # added later in this paragraph
  col = idx_col$col[-1],
  col_2ary = idx_col$col_2ary[-1])
# and sort the idxs columns
idx_col = idx_col[, c('crm', 'start', 'end')] # removes colour column
idx_lss = idx_lss[, c('crm', 'start', 'end')]
# paste an identifier to all scaffold
# (to identify dcat's scaffold_1 from dtil's scaffold_1)
idx_col$crm = paste0(sp_names[3], '@', idx_col$crm)
idx_lss$crm = paste0(sp_names[4], '@', idx_lss$crm)
# add these new scaffold names to the table,
# substituting the previous sequid names...
df_input_table$Scaffold[df_input_table$Species==sp_names[3]] = paste0(sp_names[3], '@',
  df_input_table$Scaffold[df_input_table$Species==sp_names[3]])
df_input_table$Scaffold[df_input_table$Species==sp_names[4]] = paste0(sp_names[4], '@',
  df_input_table$Scaffold[df_input_table$Species==sp_names[4]])
# remove inputted header from idx_col and idx_lss (complete species name...)
idx_col = idx_col[-1,]
idx_lss = idx_lss[-1,]
# row bind both indexes into a single data.frame
# colour had to be removed in order to have same columns
idx = rbind(idx_col, idx_lss)
# remove all orthologs that are not included in the index
# (if minor scaffolds are not in index, remove orthologs of minor scaffolds from `df`)
df_input_table = df_input_table[df_input_table$Scaffold %in% idx$crm,]
# pair colours from `highlighting_crm` with sequids
highlighting_crm$sequid = idx_col$crm


### CREATE GENOMIC LINKS and POINTS DATAFRAMES ###

# function to add members of OG to `df_links` or `df_points`
linking <- function() {
  # track two to one orthologs, so they can be added and
  # removed from CIRCOS plots easily.
  if (nrow(members) == 2) {
    o = '1:1' # two members (1:1)
  } else if (nrow(members_col)==2) {
    o = '1:2' # three members (2 in col, 1 in lss)
  } else { o = '2:1' } # 3 memb. (1 in col, 2 in lss)
  # check the quality of annotation for this orthogroup
  if ('BROC' %in% members$Method) {
    if ('Support' %in% members$BRAKER_status &
        'Good' %in% members$RNAseq_OGstatus) {
      # best quality, BROC with both BRAKER and RNA support
      q = 1
    } else if ('Support' %in% members$BRAKER_status |
        'Good' %in% members$RNAseq_OGstatus) {
      # good quality, one of BRAKER or RNA remains but not both
      q = 2
    } else {q=3}  # BROC without neither BRAKER nor RNA support
  } else {
    if ('Support' %in% members$BRAKER_status |
        'Good' %in% members$RNAseq_OGstatus) {
      # Orthofinder with at least RNA or BRAKER support
      q = 4
    } else {q=5} # worst quality, OF without any support!
  }
    new_row = data.frame(
      sequid_col = members_col$Scaffold,
      start_col = members_col$GeneStart,
      end_col = members_col$GeneEnd,
      sequid_lss = members_lss$Scaffold,
      start_lss = members_lss$GeneStart,
      end_lss = members_lss$GeneEnd,
      n_ogs_col = members_col$OGtype_perspecies,
      n_ogs_lss = members_lss$OGtype_perspecies,
      ogid = og, #DEBUG (og defined in for-loop)
      one_to_one = o,
      quality = q
    )
  # append new row to links df:
  return(rbind(df_links, new_row))}

pointing <- function() {
    # if they are a many-members orthogroup, plot as genomic points:
    # try to colour each point depending on which species
    # has more genes in the selected orthogroup
    # (compare nrows of `members_lss` versus `members_col`)
  
  # track two to one orthologs, so they can be added and
  # removed from CIRCOS plots easily.
  if (nrow(members)==3) {
    to = TRUE # 2:1
  } else {
    to = FALSE # all other cases
  }
  # check the quality of annotation for this orthogroup
  if ('BROC' %in% members$Method) {
    if ('Support' %in% members$BRAKER_status &
        'Good' %in% members$RNAseq_OGstatus) {
      # best quality, BROC with both BRAKER and RNA support
      q = 1
    } else if ('Support' %in% members$BRAKER_status |
        'Good' %in% members$RNAseq_OGstatus) {
      # good quality, one of BRAKER or RNA remains but not both
      q = 2
    } else {q=3}  # BROC without neither BRAKER nor RNA support
  } else {
    if ('Support' %in% members$BRAKER_status |
        'Good' %in% members$RNAseq_OGstatus) {
      # Orthofinder with at least RNA or BRAKER support
      q = 4
    } else {q=5} # worst quality, OF without any support!
  }
  if (nrow(members_lss)==nrow(members_col)) {
    # the same amount of OG members (genes) in both species.
    col = 'yellow'
  } else if (nrow(members_lss) > nrow(members_col)) {
    # more OG members (genes) in species colourless (dcat)
    col = c(rep('red', nrow(members_col)), rep('green', nrow(members_lss)))
  } else if (nrow(members_lss)< nrow(members_col)) {
    # more OG members (genes) in species colourful (dtil)
    col = c(rep('green', nrow(members_col)), rep('red', nrow(members_lss)))
  }
  new_rows = data.frame(
    sequid = c(members_col$Scaffold, members_lss$Scaffold),
    start = c(members_col$GeneStart, members_lss$GeneStart),
    end = c(members_col$GeneEnd, members_lss$GeneEnd),
    col = col,
    ogid = c(members_col$OGid, members_lss$OGid), #DEBUG
    two_to_one = to,
    n_ogs = c(members_col$OGtype_perspecies, members_lss$OGtype_perspecies),
    quality = q
  )
  return(rbind(df_points, new_rows))}

# resulting df with genomic links
df_links = data.frame()
# resulting df with genomic points
df_points = data.frame()
# keep track of removed and included OGs in analyses:
og.removed = 0 ; og.included = 0

# fill previous df with input df:
for (og in unique(df_input_table$OGid)) {
  # select all members of orthogroup
  members = df_input_table[df_input_table$OGid==og,]
  # subset all members from OG, for species `col` and `lss`
  members_col = members[members$Species==sp_names[3],]
  members_lss = members[members$Species==sp_names[4],]
  og.included = og.included+1
  if (nrow(members_lss)==0 | nrow(members_col)==0) {
    # it is possible to remove all genes from a family while
    # deleting rows located at minor scaffolds
    # just skip these cases...
    og.removed = og.removed+1
    og.included = og.included-1
    if (og.removed < 2) {
    print(paste0('WARNING: ', og, ' has been removed because',
    ' at least one species has all of their orthologs outside provided index sequids\n',
    '(probably all genes from OG were in minor scaffolds)'))
    } else {
    print(paste0('WARNING: ', og, ' removed (minor scaffolds)'))
    }
  } else if (nrow(members) == 2) {
    # two members: only add them to `df_links`
    # do not consider them as multigenic (`df_points`)
    df_links = linking()
  } else if (nrow(members) == 3) {
    # three members: add them to both dataframes...
    # they can be included in either data.frame, and called to different
    # desired plots depending on the tag. add tags to compartimentalize
    # plotting (multiple plots from the same data)
    df_links = linking()
    df_points = pointing()
  } else {
    df_points = pointing()
  }
}

# create a column with 'quality colour' depending on 'quality' column value
df_links$qualcol = colour_values(df_links$quality)
df_points$qualcol = colour_values(df_points$quality)

print(paste('OGs included:', og.included))
print(paste('OGs removed (minor scaffolds):', og.removed))
print(paste('Overall OGs:', length(unique(df_input_table$OGid))))

print('--> DEBUG: `df_links` and `df_points` head')
print(head(df_points, n=30))
print(head(df_links, n=30))
print('------------------------------------------')

# add colours to df_links depending on idx_col of 'origin'
# there is no need to manually create a palette, colour_values will take care of it
df_links$col = 0
df_links$col_qual = 0
for (i in c(1:nrow(highlighting_crm))) {
  # add primary colours column
  df_links$col[df_links$sequid_col==highlighting_crm[i, 1] & df_links$one_to_one=='1:1'] = highlighting_crm[i, 2]
  # add secondary colour column
  df_links$col[df_links$sequid_col==highlighting_crm[i, 1] & df_links$one_to_one %in% c('2:1', '1:2')] = highlighting_crm[i, 3]
}
print("--> DEBUG: `highlighting_crm` head()")
print(head(highlighting_crm))
print('------------------------------------')

# reformat `df_points` into `df_lines`
# create genomicLines() instead of genomicPoints()
df_lines=data.frame()
for (c in idx$crm) {
  e = idx$end[idx$crm == c]
  s = seq(0, e, window_size)
  s = s[-length(s)]
  for (i in s) {
    # select rows with start between `i` and `i+window_size` and in scaffold `c`
    df_sel = df_points[df_points$start >= i & df_points$start < i+window_size & df_points$sequid == c,]
    new_rows = data.frame(
      sequid = c,
      start = i,
      end = i+window_size,
      amount_y = nrow(df_sel[df_sel$col=='yellow', ]),
      amount_r = nrow(df_sel[df_sel$col=='red', ]),
      amount_g = nrow(df_sel[df_sel$col=='green', ]),
      # remove 2:1
      toless_amount_y = nrow(df_sel[df_sel$col=='yellow' & df_sel$two_to_one==FALSE, ]),
      toless_amount_r = nrow(df_sel[df_sel$col=='red' & df_sel$two_to_one==FALSE, ]),
      toless_amount_g = nrow(df_sel[df_sel$col=='green' & df_sel$two_to_one==FALSE, ]),
      amount_broc = nrow(df_sel[df_sel$quality %in% c(1,2,3), ]),
      amount_no_broc = nrow(df_sel[df_sel$quality %in% c(4,5), ]) )
    df_lines = rbind(df_lines, new_rows)
  }
  df_sel = df_points[i+window_size <= df_points$start & df_points$start < e & df_points$sequid == c,]
    new_rows = data.frame(
    sequid = c,
    start = i+window_size,
    end = as.numeric(e),
    amount_y = nrow(df_sel[df_sel$col=='yellow', ]),
    amount_r = nrow(df_sel[df_sel$col=='red', ]),
    amount_g = nrow(df_sel[df_sel$col=='green', ]),
    # remove 2:1
    toless_amount_y = nrow(df_sel[df_sel$col=='yellow' & df_sel$two_to_one==FALSE, ]),
    toless_amount_r = nrow(df_sel[df_sel$col=='red' & df_sel$two_to_one==FALSE, ]),
    toless_amount_g = nrow(df_sel[df_sel$col=='green' & df_sel$two_to_one==FALSE, ]),
    amount_broc = nrow(df_sel[df_sel$quality %in% c(1,2,3), ]),
    amount_no_broc = nrow(df_sel[df_sel$quality %in% c(4,5), ]) )
df_lines = rbind(df_lines, new_rows)
}
# compute the range of data (from 0 to max(df))
track_ylim = max(c(df_lines$amount_y,
                   df_lines$amount_r,
                   df_lines$amount_g))
print(paste('Most genes in `df_lines.green`:', max(df_lines$amount_g)))
print(paste('Most genes in `df_lines.red`:', max(df_lines$amount_r)))
print(paste('Most genes in `df_lines.yellow`:', max(df_lines$amount_y)))
print("--> DEBUG: `df_lines` head")
print(head(n=30, df_lines)) #DEBUG
print('--------------------------') #DEBUG


# # # # # #
# export plot to 'Rplot.pdf'
# increase default width to correctly display legend
pdf(width=9)
# # # # # #

### SUMMARY INFORMATION ###

og.total = og.removed+og.included #== length(unique(df_input_table$OGid))
og.removed.perc = round((og.removed/og.total)*100, digits=1)
og.outer = length(unique(df_points[df_points$two_to_one==FALSE, 'ogid']))
og.outer.perc = round((og.outer/og.total)*100, digits=1)
og.links = length(unique(df_links[df_links$one_to_one=='1:1', 'ogid']))
og.links.perc = round((og.links/og.total)*100, digits=1)
og.either = length(unique(df_links[df_links$one_to_one=='2:1' | df_links$one_to_one=='1:2', 'ogid']))
og.either.perc = round((og.either/og.total)*100, digits=1)
plot.new()
legend('center', title='Amount of orthologs included and excluded', bty='n',
  legend = c(
    paste0('OGs removed because all members were in minor scaffolds: ', og.removed, ', ', og.removed.perc, '%'),
    paste0('OGs only as inner links: ', og.links, ', ', og.links.perc, '%'),
    paste0('OGs either as links or in outer windows (2:1): ', og.either, ', ', og.either.perc, '%'),
    paste0('OGs in outer windows (multigenic OGs): ', og.outer, ', ', og.outer.perc, '%'),
    paste0('Total amount of inner links (not OGs but gene pairs): ', nrow(df_links)),
    paste('Most genes in a green window:', max(df_lines$amount_g)),
    paste('Most genes in a red window:', max(df_lines$amount_r)),
    paste('Most genes in a yellow window:', max(df_lines$amount_y)),
    '----------------------',
    paste('Amount of 1:1 OGs as links:', length(unique(df_links[df_links$n_ogs_lss==1 & df_links$n_ogs_col==1, 'ogid']))),
    paste('Amount of 2:1 OGs as links (more Dcat):', length(unique(df_links[df_links$n_ogs_lss==2 & df_links$n_ogs_col==1, 'ogid']))),
    paste('Amount of 1:2 OGs as links (more Dtil):', length(unique(df_links[df_links$n_ogs_lss==1 & df_links$n_ogs_col==2, 'ogid'])))
))


### PLOTTING 'STANDARD' CIRCOS FIGURE ###

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
  xlim=idx[, c('start', 'end')])
circos.track(idx$crm, ylim=c(0, track_ylim),
  track.height=.15, #bg.border=NA,
  panel.fun  = function(x,y) {
    # shorthands for a few sector properties
    cllcrm = CELL_META$sector.index
    xlim   = CELL_META$xlim
    ylim   = CELL_META$ylim
    # instead of points, draw multigenic orthogroups as
    # repetitive elements would be drawn
    # this paragraph incorporates 2:1 OGs as links (remove from points)...
    d = df_lines[df_lines$sequid==cllcrm, ]
    circos.genomicLines(d[, c('start', 'end')],
    d[, 'toless_amount_g'], type='segment', col='green')
    circos.genomicLines(d[, c('start', 'end')],
    d[, 'toless_amount_y'], type='segment', col='yellow')
    circos.genomicLines(d[, c('start', 'end')],
    d[, 'toless_amount_r'], type='segment', col='red')
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
title('Every possible link included (2 to 1 shaded darker)')


### REMOVING *OUTLIER* MULTIGENIC FAMILIES FROM THE PLOT ###
#
# Some outliers make it more difficult to understand the presented plot.
# Let's try to remove them for visualization purposes.
# Remove all rows where amount_g equals its maximum, three times
# (it could remove much more than 3 rows)
removed_max_multigenic_values = df_lines$amount_g[order(df_lines$amount_g, decreasing=TRUE)][c(1:(remove_top+1))]
for (i in c(1:remove_top)) {
  df_lines_n_out = df_lines[! df_lines$amount_g == max(df_lines$amount_g), ]
}
# recompute the new track ylim (bounded by max values of amount of genes)
track_ylim = max(c(df_lines_n_out$amount_y,
                   df_lines_n_out$amount_r,
                   df_lines_n_out$amount_g))
circos.clear()
circos.par(gap.degree = c(rep(1, nrow(idx_col)-1), 5,
  rep(1, nrow(idx_lss)-1), 5))
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors=idx$crm,
  xlim=idx[, c(2,3)])
circos.track(idx$crm, ylim=c(0, track_ylim),
  track.height=.15, #bg.border=NA,
  panel.fun  = function(x,y) {
    # shorthands for a few sector properties
    cllcrm = CELL_META$sector.index
    xlim   = CELL_META$xlim
    ylim   = CELL_META$ylim
    d = df_lines_n_out[df_lines_n_out$sequid==cllcrm, ]
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
# remove 2:1 links from plot (they are already incorporated in outside ring lines)
circos.genomicLink(df_links[, c(1:3)], df_links[, c(4:6)],
  col = df_links[, 'col'], #h.ratio = .4
)
# add a legend with the min and max of misc value,
# and the colours used to represent these.
  legend('bottomleft', title='Gene count of\nremoved outlier windows', bty='n',
    legend = c(paste0('Removed window value: ', removed_max_multigenic_values[-length(removed_max_multigenic_values)]),
    paste0('Next top value which was included: ', removed_max_multigenic_values[length(removed_max_multigenic_values)])
    ),
    cex=.8
    )
title(paste('Removing the top', remove_top,
            'densest multigenic regions (outliers) from outer windows'),
      cex.main=.9)


### REMOVING MULTIGENIC OG THAT DROP TO LINKS ###
# A few OGs have many orthologs (3, 4, etc.) in minor scaffolds
# This script removes these orthologs from the data.frame
# Then, the data.frame might 'become' 1:1 or 2:1 and
# they would be displayed as links, even though the OG is initially 4:2...
# Let's try to remove these cases
circos.clear()
circos.par(gap.degree = c(rep(1, nrow(idx_col)-1), 5,
  rep(1, nrow(idx_lss)-1), 5))
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors=idx$crm,
  xlim=idx[, c(2,3)])
circos.track(idx$crm, ylim=c(0, track_ylim),
  track.height=.15, #bg.border=NA,
  panel.fun  = function(x,y) {
    # shorthands for a few sector properties
    cllcrm = CELL_META$sector.index
    xlim   = CELL_META$xlim
    ylim   = CELL_META$ylim
    d = df_lines_n_out[df_lines_n_out$sequid==cllcrm, ]
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
# remove links that are not 1:1, 2:1, 1:2
d = df_links[
  (df_links$n_ogs_lss == 2 & df_links$n_ogs_col == 1) |
  (df_links$n_ogs_lss == 1 & df_links$n_ogs_col == 2) |
  (df_links$n_ogs_lss == 1 & df_links$n_ogs_col == 1), ]
# connect 1:3 (origin segment) to 4:6 (recipient segment)
circos.genomicLink(d[, c(1:3)], d[, c(4:6)],
  col = d[, 'col'], #h.ratio = .4
)
title(paste('Removing minor scaffold links'),
      cex.main=.9)



### REMOVE 1:1 TO BETTER VIEW 2:1 ###

circos.clear()
# separate all crms by 1 and species by 5
circos.par(gap.degree = c(rep(1, nrow(idx_col)-1), 5,
  rep(1, nrow(idx_lss)-1), 5))
# ERROR: summation of xdirection is larger than the width
# of some sectors. Remove tiny sectors or reset cell.padding as follows...
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
# initialize with sectors <- chromosome names
circos.initialize(sectors=idx$crm,
  xlim=idx[, c('start', 'end')])
circos.track(idx$crm, ylim=c(0, track_ylim),
  track.height=.15, #bg.border=NA,
  panel.fun  = function(x,y) {
    # shorthands for a few sector properties
    cllcrm = CELL_META$sector.index
    xlim   = CELL_META$xlim
    ylim   = CELL_META$ylim
    # instead of points, draw multigenic orthogroups as
    # repetitive elements would be drawn
    # this paragraph incorporates 2:1 OGs as links (remove from points)...
    d = df_lines[df_lines$sequid==cllcrm, ]
    circos.genomicLines(d[, c('start', 'end')],
    d[, 'toless_amount_g'], type='segment', col='green')
    circos.genomicLines(d[, c('start', 'end')],
    d[, 'toless_amount_y'], type='segment', col='yellow')
    circos.genomicLines(d[, c('start', 'end')],
    d[, 'toless_amount_r'], type='segment', col='red')
    circos.text(mean(xlim), mean(ylim)+mm_y(9), # add Y space (chr.name above track)
      gsub(".*d_", "", cllcrm),  # crm name that will be displayed
      # the global substitution (gsub) can trim the sp identifier or 'Scaffold' entirely
      cex=1.2, col='black', facing='inside', niceFacing=TRUE)
  }
)
# select rows not '1:1'
d = df_links[
  (df_links$n_ogs_col==1 & df_links$n_ogs_lss==2) |
  (df_links$n_ogs_col==2 & df_links$n_ogs_lss==1), ]
# colour more in species `col` as bright, more in `lss` as dark
for (i in c(1:nrow(highlighting_crm))) {
  # add bright colours to 1:2 (more in sp_names[1])
  d$col[d$sequid_col==highlighting_crm[i, 1] & d$one_to_one=='1:2'] = highlighting_crm[i, 2]
  # add dark colours to 2:1 (more in sp_names[2])
  d$col[d$sequid_col==highlighting_crm[i, 1] & d$one_to_one=='2:1'] = highlighting_crm[i, 3]
}
# connect 1:3 (origin segment) to 4:6 (recipient segment)
circos.genomicLink(d[, c(1:3)], d[, c(4:6)],
  col = d[, 'col'], #h.ratio = .4
)
title('Removes 1:1; only shows 2:1')
legend('bottomleft', title='Colour for links', bty='n',
  legend = c(paste0('(', nrow(d[d$n_ogs_lss==2, ]), ')', ' Dark links: more in ', sp_names[2]),
    paste0('(', nrow(d[d$n_ogs_col==2, ]), ')', ' Bright links: more in ', sp_names[1])),
  cex=.8
)
rm(d)


### PLOTTING LINKS FROM A SINGLE CHROMOSOME AT A TIME ###

# for each sequid
for (sequid in idx_col$crm) {
# filter out all links that are not in `sequid`
# each iteration will create a plot with links from a single sequid
links = df_links[df_links$sequid_col == sequid, ]
circos.clear()
circos.par(gap.degree = c(rep(1, nrow(idx_col)-1), 5,
  rep(1, nrow(idx_lss)-1), 5))
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors=idx$crm,
  xlim=idx[, c(2,3)])
circos.track(idx$crm, ylim=c(0, track_ylim),
  track.height=.15, #bg.border=NA,
  panel.fun  = function(x,y) {
    # shorthands for a few sector properties
    cllcrm = CELL_META$sector.index
    xlim   = CELL_META$xlim
    ylim   = CELL_META$ylim
    d = df_lines[df_lines$sequid==cllcrm, ]
    circos.genomicLines(d[, c('start', 'end')],
    d[, 'toless_amount_g'], type='segment', col='green')
    circos.genomicLines(d[, c('start', 'end')],
    d[, 'toless_amount_y'], type='segment', col='yellow')
    circos.genomicLines(d[, c('start', 'end')],
    d[, 'toless_amount_r'], type='segment', col='red')
    circos.text(mean(xlim), mean(ylim)+mm_y(9), # add Y space (chr.name above track)
      gsub(".*d_", "", cllcrm),  # crm name that will be displayed
      # the global substitution (gsub) can trim the sp identifier or 'Scaffold' entirely
      cex=1.2, col='black', facing='inside', niceFacing=TRUE)
  }
)
# connect 1:3 (origin segment) to 4:6 (recipient segment)
circos.genomicLink(links[, c(1:3)], links[, c(4:6)],
  col = links[, 'col'], #h.ratio = .4
)
title('Includes 2:1 orthologs (darker links)')
}  # acaba les iteracions de subseccionar les unions cromosomals


### REMOVING 2:1 FROM THE CENTER LINKS ('STANDARD')

circos.clear()
circos.par(gap.degree = c(rep(1, nrow(idx_col)-1), 5,
  rep(1, nrow(idx_lss)-1), 5))
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors=idx$crm,
  xlim=idx[, c(2,3)])
circos.track(idx$crm, ylim=c(0, track_ylim),
  track.height=.15, #bg.border=NA,
  panel.fun  = function(x,y) {
    # shorthands for a few sector properties
    cllcrm = CELL_META$sector.index
    xlim   = CELL_META$xlim
    ylim   = CELL_META$ylim
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
# remove 2:1 links from plot (they are already incorporated in outside ring lines)
d = df_links[df_links$one_to_one=='1:1', ]
circos.genomicLink(d[, c(1:3)], d[, c(4:6)],
  col = d[, 'col'], #h.ratio = .4
)
title('Excluding 2:1 links; incorporating them in "inner-scaffold multigenic windows"')


### COLOUR LINKS AND MULTIGENIC OGs BY QUALITY ###

track_ylim = max(c(df_lines$amount_no_broc,
                   df_lines_n_out$amount_broc))
circos.clear()
circos.par(gap.degree = c(rep(1, nrow(idx_col)-1), 5,
  rep(1, nrow(idx_lss)-1), 5))
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors=idx$crm,
  xlim=idx[, c(2,3)])
circos.track(idx$crm, ylim=c(0, track_ylim),
  track.height=.15, #bg.border=NA,
  panel.fun  = function(x,y) {
    # shorthands for a few sector properties
    cllcrm = CELL_META$sector.index
    xlim   = CELL_META$xlim
    ylim   = CELL_META$ylim
    d = df_lines[df_lines$sequid==cllcrm, ]
    circos.genomicLines(d[, c('start', 'end')],
    d[, 'amount_no_broc'], type='segment', col='red')
    circos.genomicLines(d[, c('start', 'end')],
    d[, 'amount_broc'], type='segment', col='green')
    circos.text(mean(xlim), mean(ylim)+mm_y(9), # add Y space (chr.name above track)
      gsub(".*d_", "", cllcrm),  # crm name that will be displayed
      # the global substitution (gsub) can trim the sp identifier or 'Scaffold' entirely
      cex=1.2, col='black', facing='inside', niceFacing=TRUE)
  }
)
# connect 1:3 (origin segment) to 4:6 (recipient segment)
# remove 2:1 links from plot (they are already incorporated in outside ring lines)
circos.genomicLink(df_links[, c(1:3)], df_links[, c(4:6)],
  col = df_links[, 'qualcol'], #h.ratio = .4
)
# add a legend with the min and max of misc value,
# and the colours used to represent these.
legend('bottomleft', title='Links quality', bty='n',
  legend = c(paste0("(", nrow(df_links[df_links$quality==1, ]),")", ' BROC & RNA & BRAKER'),
             paste0("(", nrow(df_links[df_links$quality==2, ]),")", ' BROC & (RNA|BRAKER)'),
             paste0("(", nrow(df_links[df_links$quality==3, ]),")", ' only BROC'),
             paste0("(", nrow(df_links[df_links$quality==4, ]),")", ' OF & (RNA|BRAKER)'),
             paste0("(", nrow(df_links[df_links$quality==5, ]),")", ' only OF'),
             'Outer multigenic BROC',
             'Outer multigenic OF'),
  fill = c(colour_values(c(1,2,3,4)), 'white', 'green', 'red'),
  cex=.7
  )
title('OG Quality')


### INTERCHROMOSOMAL QUALITY ###

track_ylim = max(c(df_lines$amount_no_broc,
                   df_lines_n_out$amount_broc))
circos.clear()
circos.par(gap.degree = c(rep(1, nrow(idx_col)-1), 5,
  rep(1, nrow(idx_lss)-1), 5))
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors=idx$crm,
  xlim=idx[, c(2,3)])
circos.track(idx$crm, ylim=c(0, track_ylim),
  track.height=.15, #bg.border=NA,
  panel.fun  = function(x,y) {
    # shorthands for a few sector properties
    cllcrm = CELL_META$sector.index
    xlim   = CELL_META$xlim
    ylim   = CELL_META$ylim
    d = df_lines[df_lines$sequid==cllcrm, ]
    circos.genomicLines(d[, c('start', 'end')],
    d[, 'amount_no_broc'], type='segment', col='red')
    circos.genomicLines(d[, c('start', 'end')],
    d[, 'amount_broc'], type='segment', col='green')
    circos.text(mean(xlim), mean(ylim)+mm_y(9), # add Y space (chr.name above track)
      gsub(".*d_", "", cllcrm),  # crm name that will be displayed
      # the global substitution (gsub) can trim the sp identifier or 'Scaffold' entirely
      cex=1.2, col='black', facing='inside', niceFacing=TRUE)
  }
)
# long list to select interchromosomal links...
# remove all links between homologous chromosomes
d = subset(df_links, !(sequid_col == 'Dtil@Scaffold_1' & sequid_lss == 'Dcat@Scaffold_5') &
  !(sequid_col == 'Dtil@Scaffold_2' & sequid_lss == 'Dcat@Scaffold_4') &
  !(sequid_col == 'Dtil@Scaffold_3' & sequid_lss == 'Dcat@Scaffold_2') &
  !(sequid_col == 'Dtil@Scaffold_4' & sequid_lss == 'Dcat@Scaffold_3') &
  !(sequid_col == 'Dtil@Scaffold_5_Vmt' & sequid_lss == 'Dcat@Scaffold_2') &
  !(sequid_col == 'Dtil@Scaffold_6' & sequid_lss == 'Dcat@Scaffold_1_Vmt') &
  !(sequid_col == 'Dtil@Scaffold_7' & sequid_lss == 'Dcat@Scaffold_1_Vmt'))
# connect 1:3 (origin segment) to 4:6 (recipient segment)
# remove 2:1 links from plot (they are already incorporated in outside ring lines)
circos.genomicLink(d[, c(1:3)], d[, c(4:6)],
  col = d[, 'qualcol'], #h.ratio = .4
)
# add a legend with the min and max of misc value,
# and the colours used to represent these.
legend('bottomleft', title='Links quality', bty='n',
  legend = c(paste0("(", nrow(d[d$quality==1, ]),")", ' BROC & RNA & BRAKER'),
             paste0("(", nrow(d[d$quality==2, ]),")", ' BROC & (RNA|BRAKER)'),
             paste0("(", nrow(d[d$quality==3, ]),")", ' only BROC'),
             paste0("(", nrow(d[d$quality==4, ]),")", ' OF & (RNA|BRAKER)'),
             paste0("(", nrow(d[d$quality==5, ]),")", ' only OF'),
             'Outer multigenic BROC',
             'Outer multigenic OF'),
  fill = c(colour_values(c(1,2,3,4)), 'white', 'green', 'red'),
  cex=.7
  )
rm(d) # elimina el data.frame temporal
title('OG Quality (removing links between homologous chromosomes)')

### ONE PLOT PER QUALITY LEVEL (1:5) ###

# prepare titles for each quality value
titles_per_quality = c('BROC & RNA & BRAKER',
  'BROC & (RNA | BRAKER)',
  'only BROC',
  'OF & (RNA | BRAKER)')
# for each quality value
for (q in unique(df_links$quality)) {
circos.clear()
circos.par(gap.degree = c(rep(1, nrow(idx_col)-1), 5,
  rep(1, nrow(idx_lss)-1), 5))
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors=idx$crm,
  xlim=idx[, c(2,3)])
circos.track(idx$crm, ylim=c(0, track_ylim),
  track.height=.15, #bg.border=NA,
  panel.fun  = function(x,y) {
    # shorthands for a few sector properties
    cllcrm = CELL_META$sector.index
    xlim   = CELL_META$xlim
    ylim   = CELL_META$ylim
    d = df_lines[df_lines$sequid==cllcrm, ]
    circos.genomicLines(d[, c('start', 'end')],
    d[, 'toless_amount_g'], type='segment', col='green')
    circos.genomicLines(d[, c('start', 'end')],
    d[, 'toless_amount_y'], type='segment', col='yellow')
    circos.genomicLines(d[, c('start', 'end')],
    d[, 'toless_amount_r'], type='segment', col='red')
    circos.text(mean(xlim), mean(ylim)+mm_y(9), # add Y space (chr.name above track)
      gsub(".*d_", "", cllcrm),  # crm name that will be displayed
      # the global substitution (gsub) can trim the sp identifier or 'Scaffold' entirely
      cex=1.2, col='black', facing='inside', niceFacing=TRUE)
  }
)
# filter out all links that are not of `q` quality
# each iteration will create a plot with links of a single quality
d = df_links[df_links$quality == q, ]
# connect 1:3 (origin segment) to 4:6 (recipient segment)
circos.genomicLink(d[, c(1:3)], d[, c(4:6)],
  col = d[, 'col'], #h.ratio = .4
)
title(titles_per_quality[q])
}  # acaba les iteracions de subseccionar les unions cromosomals


# # # # # #
# shut down plotting device
dev.off()
# # # # # #


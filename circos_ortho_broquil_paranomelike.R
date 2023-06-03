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
window_size = 5e6 # base-pairs
# si la distància entre paràlegs dels 1:2 és inferior a distancia_12...
# considera'ls "propers"
distancia_12 = 1e6
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
# subset species /Dcat/ and /Dtil/
df_input_table['Species' %in% sp_names[c(3,4)]]
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
# recupera els totals i compara'ls amb els obtinguts pel Vadim i la Marta...
cat(paste0(
  'Number of OGids in the original dataframe.tsv: 19 480\n',
  'Of these OGids, 18 617 are partially found in major scaffolds\n',
  'R data.frame has ', length(unique(df_input_table$OGid)), ' unique OGids\n',
  'R data.frame has ', nrow(df_input_table), ' number of rows\n'
  ))



# ortolegs 1:0 inclosos en finestres
df_paranomelike_10 = df_input_table[grepl("^1.0", df_input_table$OGtype) | grepl("^0.1", df_input_table$OGtype), ]
# ortòlegs N:0 inclosos en finestres
df_paranomelike_n0 = df_input_table[grepl("^..0", df_input_table$OGtype) | grepl("^0..", df_input_table$OGtype), ]
# ortolegs 2:1 inclosos com links intraespecífics (entre els 2 però no 1)
df_paranomelike_12 = df_input_table[grepl("^1.2", df_input_table$OGtype) | grepl("^2.1", df_input_table$OGtype), ]

print(head(df_paranomelike_10['OGtype']))#DEBUG
print(head(df_paranomelike_12['OGtype']))#DEBUG
print(head(df_paranomelike_n0['OGtype'], n=30))#DEBUG


### SUBSET GENOMIC LINKS DATAFRAME ###

# crea un data.frame refinat previ al circos,
# connectant els paràlegs dels ortogrups 2:1
df_links_12 = data.frame()
og.removed.links12 = 0
for (og in unique(df_paranomelike_12$OGid)) {
  # select all members of 'paragroup'
  members = df_paranomelike_12[df_paranomelike_12$OGid==og, ]
  # if one paralog from the pair has been removed from the
  # data.frame due to being in a minor scaffold... next for-loop
  if (nrow(members) != 3) {
    cat(paste('WARNING:', og, 'removed from analyses\n'))
    og.removed.links12 = og.removed.links12 + 1
    next
  }
  new_rows = data.frame(
    sequid_doub = members[members$OGtype_perspecies==2, 'Scaffold'],
    start_doub = members[members$OGtype_perspecies==2, 'GeneStart'],
    end_doub = members[members$OGtype_perspecies==2, 'GeneEnd'],
    sequid_sing = members[members$OGtype_perspecies==1, 'Scaffold'][1],
    start_sing = members[members$OGtype_perspecies==1, 'GeneStart'][1],
    end_sing = members[members$OGtype_perspecies==1, 'GeneEnd'][1],
    ogid = og
  )
  # emplena la columna PGtype, que marca la distància entre paràlegs
  if (new_rows$sequid_doub[1] != new_rows$sequid_doub[2]) {
    # la parella de paràlegs no es troba al mateix chr
    PGtype = 4 #='interchromosomal'
  } else if (abs(new_rows$start_doub[1] - new_rows$start_doub[2]) < distancia_12) {
    # els paràlegs es troben propers
    PGtype = 1 #='propers'
  } else {
    # es troben distants però mateix cromosoma
    PGtype = 2 #='distants'
  }
  new_rows$PGtype = PGtype
  df_links_12 = rbind(df_links_12, new_rows)
}
cat(paste('STATUS: as many as', og.removed.links12,
          'OGs were removed from `df_links_12` because at least one member',
          'was found outside indexed scaffolds\n'))
# add colour col depending on PGtype val...
df_links_12$PGtype_col = colour_values(df_links_12$PGtype)
# create unique pairs of these two cols for legend plotting...
legend_pgtype_colours = unique(df_links_12[, c('PGtype', 'PGtype_col')])

print(head(df_links_12))#DEBUG


### SUBSET GENOMIC WINDOWS DATAFRAME

# crea un data.frame refinat previ al circos,
# enregistrant finestres de X Mb on plotejar ortòlegs 1:0 (N:0?)
df_windows_n0 = data.frame()
for (c in idx$crm) {
  e = idx[idx$crm==c, 'end']
  s = seq(0, e, window_size)
  s = s[-length(s)]
  for (i in s) {
    df_in_s = df_paranomelike_n0[df_paranomelike_n0$GeneStart >= i & df_paranomelike_n0$GeneStart < i+window_size & df_paranomelike_n0$Scaffold == c, ]
    amount_10 = nrow(df_in_s[grepl("^1.0", df_in_s$OGtype) | grepl("^0.1", df_in_s$OGtype), ])
    new_rows = data.frame(
      sequid = c,
      start = i,
      end = i+window_size,
      # amount of genes in the selected region where OGtype is 1:0
      amount_10 = amount_10,
      # amount of genes different from 1:0 (n:0)
      amount_n0 = nrow(df_in_s)-amount_10)
    df_windows_n0 = rbind(df_windows_n0, new_rows)
  }
    df_in_s = df_paranomelike_n0[df_paranomelike_n0$GeneStart >= window_size+i & df_paranomelike_n0$GeneStart <= e & df_paranomelike_n0$Scaffold == c, ]
    amount_10 = nrow(df_in_s[grepl("^1.0", df_in_s$OGtype) | grepl("^0.1", df_in_s$OGtype), ])
    new_rows = data.frame(
      sequid = c,
      start = i,
      end = i+window_size,
      # amount of genes in the selected region where OGtype is 1:0
      amount_10 = amount_10,
      # amount of genes different from 1:0 (n:0)
      amount_n0 = nrow(df_in_s)-amount_10)
    df_windows_n0 = rbind(df_windows_n0, new_rows)
}
# add colour codes to amount_10 and more than 10 (amount_n0)
window_colours = data.frame(
  var = c('amount_10', 'amount_n0'),
  col = c('blue', 'green'))

print(head(df_windows_n0[df_windows_n0$amount_n0!=0|df_windows_n0$amount_10!=0,]))#DEBUG
print(head(df_links_12))#DEBUG

# # # # # #
# export plot to 'Rplot.pdf'
# increase default width to correctly display legend
pdf(width=9)
# # # # # #


### SUMMARY INFORMATION ###

og.total = length(unique(df_input_table$OGid))
og.included.windowsn0 = length(unique(df_paranomelike_n0$OGid))
og.included.links12 = length(unique(df_paranomelike_12$OGid)) - og.removed.links12
og.removed.links12.perc = round((og.removed.links12/og.total)*100, digits=1)
og.included.links12.perc = round((og.included.links12/og.total)*100, digits=1)
og.included.windowsn0.perc = round((og.included.windowsn0/og.total)*100, digits=1)
plot.new()
legend('center', title='General summary information', bty='n',
  legend = c(
    paste0('Total amount of unique OGs for both species: ', og.total, ', 100%'),
    paste0('1:2 OGs removed because one of its members was in minor scaffolds: ', og.removed.links12, ', ', og.removed.links12.perc, '%'),
    paste0('1:2 OGs represented as inner links: ', og.included.links12, ', ', og.included.links12.perc, '%'),
    paste0('n:0 OGs represented in outer windows: ', og.included.windowsn0, ', ', og.included.windowsn0.perc, '%'),
    paste0('Total amount of inner links (not OGs but linked pairs of genes): ', nrow(df_links_12)),
    paste('Most genes in a 1:0 window:', max(df_windows_n0$amount_10)),
    paste('Most genes in a >1:0 window:', max(df_windows_n0$amount_n0))
))


### INTENT INVERTIR CROMOSOMES ###
# Els cromosomes no es troben en l'ordre correcte. Hi ha inicis que troben
# el seu homòleg al final del cromosoma d'una altra espècie. D'aquesta manera,
# s'han de capgirar manualment si es disposa de grans blocs colinears que
# ho indiquin...
reverse_coord_windows <- function(dfw, inverted_sequids) {
  # reverse start and end columns of dfw only for selected inverted_sequids
  for (s in inverted_sequids) {
    xrange = as.double(idx[idx$sequid == s, c('start', 'end')])
    # emmascara el df segons fileres `sequid_doub` que haurien d'invertir-se
    mbool = df_windows_n0$sequid == s
    dfw[mbool, 'start'] = xrange[2] - dfw[mbool, 'end'] + xrange[1]
    dfw[mbool, 'end'] = xrange[2] - dfw[mbool, 'start'] + xrange[1]
}}

inverted_crms = c('Dtil@Scaffold_1_X',
              'Dtil@Scaffold_3',
              'Dcat@Scaffold_1_Vmt')
reverse_coord = function(x, xrange) {
  # xrange: starting and end for crm
  # x: coordinate which will be reversed
  return(xrange[2] - x + xrange[1])
}
# invert all relations for crms specified in `inverted_crms`
df_links_12.rev = df_links_12
for (s in inverted_crms) {
  xrange = c(as.numeric(idx[idx$crm == s, 'start']), as.numeric(idx[idx$crm == s, 'end']))
  # emmascara el df segons fileres `sequid_doub` que haurien d'invertir-se
  mbool_doub = df_links_12$sequid_doub == s
  df_links_12.rev$end_doub[mbool_doub] = reverse_coord(df_links_12$start_doub[mbool_doub], xrange)
  df_links_12.rev$start_doub[mbool_doub] = reverse_coord(df_links_12$end_doub[mbool_doub], xrange)
  # emmascara el df segons fileres `sequid_sing` que haurien d'invertir-se
  mbool_sing = df_links_12$sequid_sing == s
  df_links_12.rev$end_sing[mbool_sing] = reverse_coord(df_links_12$start_sing[mbool_sing], xrange)
  df_links_12.rev$start_sing[mbool_sing] = reverse_coord(df_links_12$end_sing[mbool_sing], xrange)
}
# invert all windows for crms specified in `inverted_crms`
df_windows_n0.rev = df_windows_n0
for (s in inverted_crms) {
  xrange = c(as.numeric(idx[idx$crm == s, 'start']), as.numeric(idx[idx$crm == s, 'end']))
  # emmascara el df segons fileres `sequid_doub` que haurien d'invertir-se
  mbool = df_windows_n0$sequid == s
  df_windows_n0.rev$start[mbool] = reverse_coord(df_windows_n0$end[mbool], xrange)
  df_windows_n0.rev$end[mbool] = reverse_coord(df_windows_n0$start[mbool], xrange)
}



### CREATE CIRCOS ###

# create 2 plots, one for Dtil doubles and another for Dcat doubles...
for (pattern in c("^Dtil@", "^Dcat@")) {
# compute the range of data for outer windows (from 0 to max(df))
track_ylim = max(c(df_windows_n0$amount_n0,
                   df_windows_n0$amount_10))
# good habit to clear previous options
circos.clear()
# separate all crms by 1 and species by 5
circos.par(gap.degree = c(rep(1, nrow(idx_col)-1), 5,
  rep(1, nrow(idx_lss)-1), 5))
# ERROR: summation of xdirection is larger than the width
# of some sectors. Remove tiny sectors or reset cell.padding as follows...
circos.par(cell.padding = c(0, 0, 0, 0))
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
    d = df_windows_n0.rev[df_windows_n0.rev$sequid==cllcrm, ]
    circos.genomicLines(d[, c('start', 'end')],
      d[, 'amount_n0'],
      col=window_colours[window_colours$var=='amount_n0', 'col'])
    circos.genomicLines(d[, c('start', 'end')],
      d[, 'amount_10'],
      col=window_colours[window_colours$var=='amount_10', 'col'])
    circos.text(mean(xlim), mean(ylim)+mm_y(9), # add Y space (chr.name above track)
      gsub(".*d_", "", cllcrm),  # crm name that will be displayed
      # the global substitution (gsub) can trim the sp identifier or 'Scaffold' entirely
      cex=1.2, col='black', facing='inside', niceFacing=TRUE)
  }
)
# connect 1:3 (origin segment) to 4:6 (recipient segment)
x = df_links_12.rev[grepl(pattern, df_links_12$sequid_doub), ]
circos.genomicLink(x[, c(1:3)], x[, c(4:6)],
  col = x[, 'PGtype_col'], #h.ratio = .4
)
# texte recompte esdeveniments:
t = paste0('\n',
  'Genes 1:0 for Dtil: ', sum(df_windows_n0[grepl("^Dtil@", df_windows_n0$sequid), 'amount_10']), '\n',
  'Genes 1:0 for Dcat: ', sum(df_windows_n0[grepl("^Dcat@", df_windows_n0$sequid), 'amount_10']), '\n',
  'Genes n:0 for Dtil: ', sum(df_windows_n0[grepl("^Dtil@", df_windows_n0$sequid), 'amount_n0']), '\n',
  'Genes n:0 for Dcat: ', sum(df_windows_n0[grepl("^Dcat@", df_windows_n0$sequid), 'amount_n0']), '\n')
plot_dim = par('usr')
# Coordenades x i y del punt superior a la dreta...
x_topright <- plot_dim[2]
y_topright <- plot_dim[4]
# Espai vertical necessari per al text
text_height <- strheight(t, cex=0.7)
# Coordenades ajustades per col·locar el text a top right
x_text <- x_topright - text_height
y_text <- y_topright - text_height
# afegir el text
text(x_text, y_text, t, cex=0.7, pos=3)

legend('bottomleft', bty='n', cex=.7,
       legend = c(
          paste0('Outer windows of 1:0 (', sum(df_windows_n0$amount_10),' genes)'),
          paste0('Outer windows of >1:0 (', sum(df_windows_n0$amount_n0),' genes)'),
          paste0('Closeby 2:1 paralogs (dist<', distancia_12,') (amount=', nrow(x[x$PGtype==1,]),')' ),
          paste0('Distant 2:1 paralogs (dist>=', distancia_12,') (amount=', nrow(x[x$PGtype==2,]), ')' ),
          paste0('Interchromosomal 2:1 paralogs (amount=', nrow(x[x$PGtype==4,]), ')')
        ),
       fill = c(
          window_colours[window_colours$var=='amount_10', 'col'],
          window_colours[window_colours$var=='amount_n0', 'col'],
          legend_pgtype_colours[legend_pgtype_colours$PGtype==1, 'PGtype_col'], #propers
          legend_pgtype_colours[legend_pgtype_colours$PGtype==2, 'PGtype_col'], #distants
          legend_pgtype_colours[legend_pgtype_colours$PGtype==4, 'PGtype_col']  #interchromosomal
        )
)
species_loop = paste0(strsplit(pattern, '')[[1]][2:5], collapse='')
title(paste0('Links: Pairs of ', species_loop,
             "'s paralogs -- Outer windows: n:0 gene-count"))
}



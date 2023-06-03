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
             "Species-complete     species-shorthand           ",
             "Scaffold_1           100",
             "Scaffold_2           50",
             "etc.",
             "`sp2_chr.idx` follows the same format but does not require colours",
             sep='\n'),
       call.=FALSE)
}
# Path to orthologs.tsv dataframe file
orthologs_file = args[1]
# Path to a tabulated index files
index_file1 = args[2]
index_file2 = args[3]
# Adjust the window size
window_size = 4e6
# OGtype boundary
llindar_casos_especials = 300

library(circlize)
library(colourvalues)
library(gridExtra)
library(grid)
library(gtable)

# table with data about families of orthologs
# row.names=1 specifies that the first col is the pandas index df. col
df_input_table = read.table(orthologs_file,
	sep='\t', header=TRUE, row.names=1)
# index of sequids (their length) in the previous table
idx1 = read.table(index_file1,
	sep='\t', header=FALSE)[,c(1,2)]
idx2 = read.table(index_file2,
	sep='\t', header=FALSE)
# species names incorporated in index
# first row of file is sp. name
sp_names = c(idx1[1,1], idx2[1,1], idx1[1,2], idx2[1,2])
# Paste an species identifier (shorthand in idx) to sequids,
# so scaffold_1 from one species can be fished between scaffold_1 from another species
idx1[c(2:nrow(idx1)), 1] = paste0(sp_names[3], '@', idx1[c(2:nrow(idx1)), 1])
idx2[c(2:nrow(idx2)), 1] = paste0(sp_names[4], '@', idx2[c(2:nrow(idx2)), 1])
idx = rbind(idx1[-1,], idx2[-1,])   # remove species names from general idx
names(idx) = c('sequid', 'end')     # rename columns
idx$start = 0                       # starting position is zero for all sequids
idx$end = as.numeric(idx$end)       # make sure end is numeric
# Change scaffold names from table as it has been done with IDX
df_input_table$Scaffold[df_input_table$Species==sp_names[3]] = paste0(sp_names[3], '@',
  df_input_table$Scaffold[df_input_table$Species==sp_names[3]])
df_input_table$Scaffold[df_input_table$Species==sp_names[4]] = paste0(sp_names[4], '@',
  df_input_table$Scaffold[df_input_table$Species==sp_names[4]])
# remove all orthologs that are not included in the index
# (if minor scaffolds are not in index, remove orthologs of minor scaffolds from `df`)
removed.minor = nrow(df_input_table[!df_input_table$Scaffold %in% idx$sequid, ])
df_input_table = df_input_table[df_input_table$Scaffold %in% idx$sequid, ]

cat(paste0(
  'Number of OGids in the original dataframe.tsv: 19 480 (?)\n',
  'Of these OGids, 18 617 (?) are partially found in major scaffolds\n',
  'R data.frame has ', length(unique(df_input_table$OGid)), ' unique OGids\n',
  'R data.frame has ', nrow(df_input_table), ' number of rows\n',
  'Removed ', removed.minor, ' genes because they were located in minor scaffolds\n'
  ))


### DETECT OWN AND OTHER OGTYPE ###
# Semi-manual (how to split OGtype col has to be specified)
# detection of your own vs their other Ogtype

# for Dtil
m <- df_input_table$Species == "Dtil"
split_ogtypes = strsplit(df_input_table$OGtype[m], '\\.') # split ogtypes by '\\.'
# sapply the second split character to the column `own`
df_input_table$own[m] <- sapply(split_ogtypes, function(x) x[2])
# sapply the first split character to the column `other`
df_input_table$other[m] <- sapply(split_ogtypes, function(x) x[1])

# repeat for Dcat
m <- df_input_table$Species == "Dcat"
split_ogtypes <- strsplit(df_input_table$OGtype[m], '\\.')
df_input_table$own[m] <- sapply(split_ogtypes, function(x) x[1])
df_input_table$other[m] <- sapply(split_ogtypes, function(x) x[2])
df_input_table$own = as.numeric(df_input_table$own)
df_input_table$other = as.numeric(df_input_table$other)

# Quants casos per sobre del llindar
# de casos especials existeixen?
m <- (df_input_table$own >= llindar_casos_especials) | (df_input_table$own >= llindar_casos_especials)
casos_especials = unique(df_input_table[m, 'OGtype'])
metode_casos_especials = unique(df_input_table[m, 'Method'])
casos_especials = paste(metode_casos_especials, ' - ', casos_especials, collapse='\n')
print(casos_especials)


### SLIDING WINDOWS ###
# compute genes found in sliding windows over the genomes
# classify genes depending of OGtype
df_windows = data.frame()
mask_special_cases = logical()
for (c in idx$sequid) {
  e = idx[idx$sequid == c, 'end']  # end for sequid `c`
  # sequence of windows between 0 and e of size `window_size`
  s = seq(0, e, window_size)
  for (i in s) {
    df_in_i = df_input_table[df_input_table$GeneStart >= i &
      df_input_table$GeneStart < i+window_size &
      df_input_table$Scaffold == c, ]
    #print(typeof(e)) ; print(typeof(i)) ; print(typeof(window_size+1)) ; print(typeof(c))#DEBUG
    #print(head(df_in_i[,'OGtype']))#DEBUG
    new_rows = data.frame(
      sequid = c,
      start = i,
      # torna final del cromosoma si la finestra (i+window_size) sobreeix
      end = as.numeric(ifelse((i+window_size) < e, i+window_size, e)),
      amount_1_0 =
        nrow(df_in_i[df_in_i$own==1 & df_in_i$other==0, ]),
      amount_mt1_0 =
        nrow(df_in_i[df_in_i$own >1 & df_in_i$other==0, ]),
      amount_1_1 =
        nrow(df_in_i[df_in_i$own==1 & df_in_i$other==1, ]),
      amount_2_1 =
        nrow(df_in_i[df_in_i$own==2 & df_in_i$other==1, ]),
      amount_mt2_1 =
        nrow(df_in_i[df_in_i$own >2 & df_in_i$other==1, ]),
      amount_1_2 =
        nrow(df_in_i[df_in_i$own==1 & df_in_i$other==2, ]),
      amount_1_mt2 =
        nrow(df_in_i[df_in_i$own==1 & df_in_i$other >2, ]),
      amount_mt1_mt1 =
        nrow(df_in_i[df_in_i$own >1 & df_in_i$other >1, ]),
      amount_total =
        nrow(df_in_i)
      )
    df_windows = rbind(df_windows, new_rows)
    mask_special_cases = c(mask_special_cases,
      any(df_in_i$own >= llindar_casos_especials |
        df_in_i$other >= llindar_casos_especials))
  }}
#print(mask_special_cases)#DEBUG
#print(df_in_i$own)#DEBUG
#print(df_in_i$other)#DEBUG


### INTENT INVERTIR CROMOSOMES ###
# Els cromosomes no es troben en l'ordre correcte. Hi ha inicis que troben
# el seu homòleg al final del cromosoma d'una altra espècie. D'aquesta manera,
# s'han de capgirar manualment si es disposa de grans blocs colinears que
# ho indiquin...
#### Tot i l'esforç, com que no es creen links no feia falta... no?
#### O links des del final fins el principi, en comptes d'unint centres de sectors?
inverted_sequids = c(
  'Dtil@Scaffold_1_X',
  'Dtil@Scaffold_3',
  'Dcat@Scaffold_1_Vmt')
reverse_coord_windows <- function(dfw, inverted_sequids) {
  # reverse start and end columns of dfw only for selected inverted_sequids
  for (s in inverted_sequids) {
    xrange = as.double(idx[idx$sequid == s, c('start', 'end')])
    # emmascara el df segons fileres `sequid_doub` que haurien d'invertir-se
    mbool = dfw[, 'sequid'] == s
    dfw.rev = dfw
    dfw.rev[mbool, 'start'] = xrange[2] - dfw[mbool, 'end'] + xrange[1]
    dfw.rev[mbool, 'end'] = xrange[2] - dfw[mbool, 'start'] + xrange[1]
  }
  return(dfw.rev)
}
df_windows.rev = reverse_coord_windows(df_windows, inverted_sequids)
#print(df_windows.rev)#DEBUG


### OVERALL WINDOW STATS ###
# Extreu les finestres de valor màxim, median, etc. per a cada tipus d'ortòlegs
# Al mateix temps, relativitza entre un rang 0-1 totes les quantitats
df_windows.stats = data.frame()
# Compute stats independently for each species
df_windows.stats.spec1 = data.frame()
df_windows.stats.spec2 = data.frame()
# Compute outliers for each OGtype in `df_windows`
df_windows.outliers = df_windows
df_windows.outliers[,c(4:length(df_windows))] = FALSE

# Calculate outliers with [IQR](https://en.wikipedia.org/wiki/Interquartile_range)
outliers <- function(val_list) {
  q1 = quantile(val_list, 0.25)
  q3 = quantile(val_list, 0.75)
  iqr = q3 - q1
  # Limits of IQR which bound non-outliers
  lower = q1 - 1.5 * iqr
  upper = q3 + 1.5 * iqr
  # Return outliers not bound by previous boundary
  return(val_list[val_list < lower | val_list > upper])
}

# For each column of values in df_windows...
for (i in c(4:length(df_windows))) {
  # Taula overall stats
  cat(paste0("STATUS: Max value for complete df and colname ", names(df_windows)[i], ': ', max(df_windows[,i]), '\n'))
  df_windows.stats[names(df_windows)[i], 'max'] = max(df_windows[,i])
  df_windows.stats.spec1[names(df_windows)[i], 'max'] = max(df_windows[grepl(paste0("^", sp_names[3], "@"), df_windows$sequid), i])
  df_windows.stats.spec2[names(df_windows)[i], 'max'] = max(df_windows[grepl(paste0("^", sp_names[4], "@"), df_windows$sequid), i])
  cat(paste0("STATUS: Mean value for complete df and colname ", names(df_windows)[i], ': ', mean(df_windows[,i]), '\n'))
  df_windows.stats[names(df_windows)[i], 'mean'] = round(mean(df_windows[,i]), digits=3)
  df_windows.stats.spec1[names(df_windows)[i], 'mean'] = round(mean(df_windows[grepl(paste0("^", sp_names[3], "@"), df_windows$sequid), i]), digits=3)
  df_windows.stats.spec2[names(df_windows)[i], 'mean'] = round(mean(df_windows[grepl(paste0("^", sp_names[4], "@"), df_windows$sequid), i]), digits=3)
  cat(paste0("STATUS: Median value for complete df and colname ", names(df_windows)[i], ': ', median(df_windows[,i]), '\n'))
  df_windows.stats[names(df_windows)[i], 'median'] = median(df_windows[,i])
  df_windows.stats.spec1[names(df_windows)[i], 'median'] = median(df_windows[grepl(paste0("^", sp_names[3], "@"), df_windows$sequid), i])
  df_windows.stats.spec2[names(df_windows)[i], 'median'] = median(df_windows[grepl(paste0("^", sp_names[4], "@"), df_windows$sequid), i])
  cat(paste0("STATUS: Sum of these occurrences for complete df and colname ", names(df_windows)[i], ': ', sum(df_windows[,i]), '\n'))
  df_windows.stats[names(df_windows)[i], 'sum'] = sum(df_windows[,i])
  df_windows.stats.spec1[names(df_windows)[i], 'sum'] = sum(df_windows[grepl(paste0("^", sp_names[3], "@"), df_windows$sequid), i])
  df_windows.stats.spec2[names(df_windows)[i], 'sum'] = sum(df_windows[grepl(paste0("^", sp_names[4], "@"), df_windows$sequid), i])
  # Calculate outliers with IQR
  cat(paste0("STATUS: Found ", length(outliers(df_windows[,i])), " outlier windows in complete df and colname ", names(df_windows)[i], '\n'))
  df_windows.stats[names(df_windows)[i], 'n_outliers'] = length(outliers(df_windows[,i]))
  df_windows.stats.spec1[names(df_windows)[i], 'n_outliers'] = length(outliers(df_windows[grepl(paste0("^", sp_names[3], "@"), df_windows$sequid), i]))
  df_windows.stats.spec2[names(df_windows)[i], 'n_outliers'] = length(outliers(df_windows[grepl(paste0("^", sp_names[4], "@"), df_windows$sequid), i]))

  # Marca finestres outliers amb `TRUE` ; d'aquestes en farem un genomicPoints stacked de circlize
  df_windows.outliers[df_windows[,i] %in% outliers(df_windows[,i]), i] = TRUE
  # Computa el punt al mig de la pista
  df_windows.outliers[, 'middle'] = rowMeans(df_windows.outliers[, c('start', 'end')])
  # Finalment, relativitza entre un rang 0-1 totes les quantitats
  # permet comparar tipus ortòlegs entre ells al CIRCOS track (mateix track ylimit)
  df_windows[,i]= df_windows[,i]/max(df_windows[,i])
}

#print(df_windows.outliers)#DEBUG
#print(df_windows)#DEBUG
#print(df_windows.stats)#DEBUG
#print(df_windows.stats.spec1)#DEBUG
#print(df_windows.stats.spec2)#DEBUG


### PRINT GENERAL STATS TO PDF IN DF FORMAT ###

pdf(width=9)
# One title for each drawn table
table_title = textGrob("Genome distribution of OGtypes in sliding windows\n(max members of OGtype in a single window, etc.)", gp=gpar(fontsize=14, fontface='bold'))
# Add more visually appealing rownames
rownames(df_windows.stats) = c('1:0', '>=2:0', '1:1', '2:1', '>=3:1', '1:2', '1:>=3', '>=2:>=2', 'All OGtypes')
rownames(df_windows.stats.spec1) = c('1:0', '>=2:0', '1:1', '2:1', '>=3:1', '1:2', '1:>=3', '>=2:>=2', 'All OGtypes')
rownames(df_windows.stats.spec2) = c('1:0', '>=2:0', '1:1', '2:1', '>=3:1', '1:2', '1:>=3', '>=2:>=2', 'All OGtypes')
# Convert data.frames into "tableGrob" objects (from gridExtra packet)
# ajusta la mida amb `base_size` i els mm de padding de cada cel·la (més ample que alt)
taula1 <- tableGrob(df_windows.stats, theme = ttheme_default(base_size = 11,
                                                core = list(padding=unit(c(10, 1), "mm")) ))
taula2 <- tableGrob(df_windows.stats.spec1, theme = ttheme_default(base_size = 9,
                                                core = list(padding=unit(c(10, 1), "mm")) ))
taula3 <- tableGrob(df_windows.stats.spec2, theme = ttheme_default(base_size = 9,
                                                core = list(padding=unit(c(10, 1), "mm")) ))
# Afegeix titol al peu de cada taula
title_grub = function(title_string, table_grob) {
  title = textGrob(title_string, x=0, hjust=0, gp=gpar(fontface='italic', fontsize=11))
  padding = unit(0.5, "line")
  table_grob = gtable_add_rows(table_grob,
    heights = grobHeight(title)+padding)
  table = gtable_add_grob(table_grob,
    list(title), t=nrow(table_grob), l=2, r=ncol(table_grob))
  return(table)
}
taula1 = title_grub(paste0('Sliding windows of size ', as.character(window_size), ' throughout both species'), taula1)
taula2 = title_grub(paste0('Subsetting ', sp_names[1], ' from complete data.frame'), taula2)
taula3 = title_grub(paste0('Subsetting ', sp_names[2], ' from complete data.frame'), taula3)
# Join tables in the specified matrix (first table spans two rows, repeated one)
grid.arrange(
  taula1, taula2, taula3,
  textGrob(paste(c('Casos especials', 'representats amb', 'punts vermells:', casos_especials), collapse='\n'),
                  gp=gpar(fontsize=8)),
  layout_matrix = rbind(c(1,1,1,4), c(2,2,3,3)),
  top = table_title)


### PLOT THE CIRCOS IN THE FOLLOWING PAGES OF THE SAME PDF ###

lines_colours = colour_values(c(1:(length(df_windows)-3)))
circos.clear()
# separate all crms by 1 and species by 5
circos.par(gap.degree = c(rep(1, nrow(idx1[-1,])-1), 5,
  rep(1, nrow(idx2[-1,])-1), 5))
# ERROR: summation of xdirection is larger than the width
# of some sectors. Remove tiny sectors or reset cell.padding as follows...
circos.par(cell.padding = c(0, 1, 0, 1),
    track.margin=c(0,0))
# initialize with sectors <- chromosome names
circos.initialize(sectors=idx$sequid,
  xlim=idx[, c('start', 'end')])
# create one track per variable/column
for (i in c(4:length(df_windows))) {
circos.genomicTrack(df_windows[, c(1:3, i)],
  track.height=.06, #bg.border=NA,
  panel.fun  = function(region, value, ...) {
    # shorthands for a few sector properties
    cllcrm = CELL_META$sector.index
    xlim   = CELL_META$xlim
    ylim   = CELL_META$ylim
    # instead of points, draw multigenic orthogroups as
    # repetitive elements would be drawn
    # this paragraph incorporates 2:1 OGs as links (remove from points)...
    circos.genomicLines(region, value, col=lines_colours[(i-3)])
    if (i==4) {
    circos.text(mean(xlim), mean(ylim)+mm_y(5), # add Y space (chr.name above track)
      gsub(".*d_", "", cllcrm),  # crm name that will be displayed
      # the global substitution (gsub) can trim the sp identifier or 'Scaffold' entirely
      cex=1.2, col='black', facing='inside', niceFacing=TRUE)
    }
  }
)
}

# list of points (outliers) dataframes; one dataframe per variable
outlier_list = list()
df_no_outliers = data.frame()
for (column in c(4:12)) {
  d = df_windows.outliers[df_windows.outliers[, column], c('sequid', 'start', 'end')]
  d$height = 14 - column
  new_rows = data.frame(
    sequid = unique(df_windows.outliers[!df_windows.outliers$sequid %in% d$sequid, 'sequid']),
    height = 14 - column)
  df_no_outliers = rbind(df_no_outliers, new_rows)
  outlier_list[[as.character(column)]] <- d
}
#print(df_no_outliers)#DEBUG

# place 'special cases' in outlier_list...
d = df_windows.outliers[mask_special_cases, c('sequid', 'start', 'end')]
d$height = 1
if (any(!df_windows.outliers$sequid %in% d$sequid)) {
new_rows = data.frame(
  sequid = unique(df_windows.outliers[!df_windows.outliers$sequid %in% d$sequid, 'sequid']),
  height = 1)
}
outlier_list[[as.character(1)]] <- d
df_no_outliers = rbind(df_no_outliers, new_rows)
lines_colours = c(lines_colours, 'red')# colour special cases red


# plot points where outliers are, stacking variables in the same track
# draw a dashed line for each var. Y-level
# (draw line first so it is placed behind the points)
circos.genomicTrack(outlier_list, stack=TRUE,
  track.height=.12,
  panel.fun = function(region, value, ...) {
    i = getI(...)
    line_height = 11 - i
    circos.lines(CELL_META$cell.xlim, c(line_height, line_height), lty=1, col="#00000040")
    circos.genomicPoints(region, value,
      col=lines_colours[i], pch=16, cex=0.22)
  }
)

# draw in-sector lines that contain no outliers
for (row in c(1:nrow(df_no_outliers))) {
  s = df_no_outliers[row, 'sequid']
  h = df_no_outliers[row, 'height']
  x = unlist(idx[idx$sequid == s, c('start', 'end')])
  #print(s); print(h); print(x) #DEBUG
  circos.lines(sector.index= s,
    y= c(h, h),
    x= x,
    col= '#00000040', lty=3)}

# create single links between centers of homologous chromosomes
# manually specify pairs of homologous chromosomes
order_homol = c(idx$sequid[9], idx$sequid[12],
               idx$sequid[8], idx$sequid[11],
               idx$sequid[8], idx$sequid[10],
               idx$sequid[10])
pos_homol = c(rowMeans(idx[9, c('start', 'end')]),
              rowMeans(idx[12, c('start', 'end')]),
              rowMeans(idx[8, c('start', 'end')]),
              rowMeans(idx[11, c('start', 'end')]),
              rowMeans(idx[8, c('start', 'end')]),
              rowMeans(idx[10, c('start', 'end')]),
              rowMeans(idx[10, c('start', 'end')]))

links = data.frame(
  sequid_til= idx$sequid[c(1:7)],
  pos_til= rowMeans(idx[c(1:7), c('start', 'end')]),
  sequid_cat= order_homol,
  pos_cat= pos_homol)
for (row in c(1:nrow(links))) {
circos.link(links[row, 'sequid_til'], links[row, 'pos_til'],
            links[row, 'sequid_cat'], links[row, 'pos_cat'])
}

legend('bottomleft', bty='n', cex=.7,
       title=paste0("Legend\n", "(same order\n", "as tracks)"),
       legend = c('1:0', '>=2:0', '1:1',
                  '2:1', '>=3:1', '1:2',
                  '1:>=3', '>=2:>=2', 'All OGtypes',
                  paste0('Special cases\n(>=', llindar_casos_especials, ':>=0 or vice versa)')),
       fill = lines_colours)

title(paste0('Gene content per OGtypes in sliding windows of ', window_size))


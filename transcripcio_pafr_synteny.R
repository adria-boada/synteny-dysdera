
## MAIN PAF FILENAME ##

help = paste0("\n1. Supply a PAF file as the first arg. value. ",
              "\n\n2. Supply a TSV (TAB-separated values) file pairing ",
              "chromosome IDs (sequids) and colours ",
              "as the second arg. value. For instance:\n\n",
              "```\n",
              "Dtilchr1\t#ff6dd3\n",
              "DtilchrX\t#ff7777\n",
              "etc...\n",
              "```")

# Load the path to the PAF file with the necessary alignments:
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
  stop(help)
}
path_to_paf = args[1]
path_to_colours = args[2]

## OPTIONAL VARIABLES ##

# Discard mappings with a mapping quality below `min_map_qual`. Mapping quality
# ranges from 0 to 60 (do not input a value outside this threshold).
min_map_qual = 0
# Discard mappings with an alignment length below `min_align_len` (in base
# pairs).
min_align_len = 600
# The following patterns should be found in Qname and Tname strings:
  # Regex patterns which identify only minor scaffolds (all of them):
  pattern_scaff = "Scaffold|ctg"
  # Regex patterns which identify only major chromosomes (all of them):
  pattern_chr = "chr"
# Specify sequids where syntenic relations are arranged orthogonally. Their
# coordinates will be inverted in order for syntenic relations to be arranged
# in parallel.
inverted_sequids = c(
  "DtilchrX",
  "Dtilchr1",
  "Dcatchr2"
)
# Specify species.
species_query = ""
species_query.col = "#eb8f46"
species_target = ""
species_target.col = "#30acac"

## LIBRARIES and ARGUMENTS ##

library(pafr)
library(circlize)
# Read the PAF file
print(paste0("Reading file at ", path_to_paf))
raw_aln_circos = read_paf(path_to_paf)
# Read a file with a specific format which pairs chromosomes with colours.
colour_pairs = read.table(path_to_colours, header=FALSE, sep="\t")
names(colour_pairs) = c("sequid", "colour")

## LOAD AND FILTER DATAFRAME (ie PAF) ##
# The filtered dataframe (aln_circos) excludes minor scaffolds, low mapping
# quality, short alignments and secondary mappings.

# Label chromosomes with "Q." and "T." at the beginning of the sequence ID.
raw_aln_circos$tname = paste0('T.', raw_aln_circos$tname)
raw_aln_circos$qname = paste0('Q.', raw_aln_circos$qname)
# Filtering (keep tabs on initial rows).
initial_rows = nrow(raw_aln_circos)
# Count (later on filter out) alignments with both minor scaffolds.
df_double_minor_scaffold = raw_aln_circos[
  grepl(pattern_scaff, raw_aln_circos$qname) &
  grepl(pattern_scaff, raw_aln_circos$tname), ]
removed.both_min_scaffold = nrow(df_double_minor_scaffold)
# Count (later on filter out) alignments with a single minor scaffold.
df_minor_scaffold = raw_aln_circos[
  grepl(pattern_scaff, raw_aln_circos$qname) |
  grepl(pattern_scaff, raw_aln_circos$tname), ]
df_single_scaffold = df_minor_scaffold[
  grepl(pattern_chr, df_minor_scaffold$qname) |
  grepl(pattern_chr, df_minor_scaffold$tname), ]
removed.single_min_scaffold = nrow(df_single_scaffold)
# Remove mappings which involve any minor scaffold at all.
aln_circos = raw_aln_circos
aln_circos = raw_aln_circos[
  # Sequid names do not (!) contain "Scaffold"
  # (at least one of these names is a chromosome)
  ! grepl(pattern_scaff, raw_aln_circos$qname) &
  ! grepl(pattern_scaff, raw_aln_circos$tname), ]
# Remove low mapping quality.
aln_circos = subset(aln_circos, mapq >= min_map_qual)
df_single_scaffold = subset(df_single_scaffold, mapq >= (min_map_qual))
removed.mapq = (nrow(raw_aln_circos) -
                nrow(subset(raw_aln_circos, mapq >= min_map_qual)))
# Remove short mappings (below specified alignment length)
aln_circos = subset(aln_circos, alen > min_align_len)
df_single_scaffold = subset(df_single_scaffold, alen > (min_align_len))
removed.alilen = (nrow(raw_aln_circos) -
                  nrow(subset(raw_aln_circos, alen > min_align_len)))
# Remove secondary alignments.
aln_circos = filter_secondary_alignments(aln_circos)
removed.secali = (nrow(raw_aln_circos) -
                  nrow(filter_secondary_alignments(raw_aln_circos)))
# Total of rows removed at the end...
removed.all = nrow(raw_aln_circos) - nrow(aln_circos)
included.single_min_scaffold.qualfiltered = nrow(df_single_scaffold)

# Print these previous values to the terminal screen.
  t = initial_rows
  cat(paste0('Complete/raw amount of rows: ', nrow(raw_aln_circos), '\n'))
  cat(paste0('Excluded rows because a minor scaffold maps to another minor scaffold: ',
             round((removed.both_min_scaffold/t)*100, digits=1), ' %\n'))
  cat(paste0('Rows where one minor scaffold maps to a main chromosome: ',
             round((removed.single_min_scaffold/t)*100, digits=1), ' %\n'))
  cat(paste0('Included rows where one minor scaffold maps to a main chromosome: ',
             round((included.single_min_scaffold.qualfiltered/t)*100, digits=1), ' %\n'))
  cat(paste0('Excluded rows because their mapQ. was lower than ',
             min_map_qual, ': ', round((removed.mapq/t)*100, digits=1), ' %\n'))
  cat(paste0('Excluded rows because their alig. length was lower than ',
             min_align_len, ': ', round((removed.alilen/t)*100, digits=1), ' %\n'))
  cat(paste0('Excluded rows because they were secondary mappings: ',
             round((removed.secali/t)*100,digits=1), ' %\n'))
  cat(paste0('The combined percentage of original rows that have been removed is: ',
             round((removed.all/t)*100,digits=1), ' %\n'))

# Read the second argument, which pairs chromosomes with a hexadecimal colour label
# For instance, an example data.frame:
  ##  colour_pairs = data.frame(
  ##    sequid = c("DtilchrX", "Dtilchr1", "Dtilchr2",
  ##               "Dtilchr3", "Dtilchr4", "Dtilchr5",
  ##               "Dtilchr6"),
  ##    colors = c("#ff6dd3", "#ff0008", "#5bc358",
  ##               "#ffbb54", "#ff7777", "#55e8f0",
  ##               "#b9dae9"))
# Create a `colours` column in the main alignment data.frame.
# Initialise as black (unspecified colours in `colour_pairs`
# will be black).
aln_circos$colours <- "black"
# Fill the newly created column `colours`:
for (i in 1:nrow(colour_pairs)) {
  colour = colour_pairs[i, "colour"]
  sequid = colour_pairs[i, "sequid"]
  # Find rows that match `sequid` in either `qname` or `tname` columns
  aln_circos[
             grepl(sequid, aln_circos$qname) |
             grepl(sequid, aln_circos$tname),
           # Substitute the `colour` column of the selected rows by sequid.
           "colours"] = colour
}

# Sort alignments by mapping quality. Consequently, the worst mappings will
# be plotted first, and thus hidden underneath the best mappings, which will
# be plotted afterwards in succession. Secondarily, sort by alignment length.
aln_circos = aln_circos[order(aln_circos$mapq, aln_circos$alen), ]
# To see the differences applied by this setting, create a CIRCOS plot in which
# the worst alignments are plotted above the rest.
worst_aln_first = aln_circos[order(aln_circos$mapq, aln_circos$alen,
                                   decreasing=TRUE), ]

# Search for regions of main sequids which are mapped by minor scaffolds.
query_highlight = df_single_scaffold[grepl("Scaffold|ctg", df_single_scaffold$tname),
                                     c("qname", "qstart", "qend", "alen", "mapq")]
names(query_highlight) = c("sequid", "start", "end", "alen", "mapq")
target_highlight = df_single_scaffold[grepl("Scaffold|ctg", df_single_scaffold$qname),
                                     c("tname", "tstart", "tend", "alen", "mapq")]
names(target_highlight) = c("sequid", "start", "end", "alen", "mapq")
# Concatenate target and query highlighting bed-like data.frame:
outer_highlight = rbind(query_highlight, target_highlight)

## DOTPLOTS by pafr ##

# Function to return a dataframe where column X has been reduced to its last
# character. It allows qname and tname columns to be reduced to single
# characters instead of the full sequid name.
reduce_column = function(df, col) {
  df[, col] = substr(df[, col],
  start=nchar(df[, col]),
  stop=nchar(df[, col]))
  return(df)
}
# Ggplot theme for the dotplots. Remove background lines.
theme = theme(panel.background = element_rect(fill = "white", colour = "grey50"))

print(paste0("Writing dotplots to the PDF `",
             path_to_paf, ".dotplots.pdf`"))
pdf(paste0(path_to_paf, ".dotplots.pdf"))
# Create a dotplot from the filtered PAF file.
# Change qname and tname columns to be the single characters
# at the end of the string of text in those columns:
aln_circos.dotplots = aln_circos
aln_circos.dotplots = reduce_column(aln_circos.dotplots, "qname")
aln_circos.dotplots = reduce_column(aln_circos.dotplots, "tname")
pafr::dotplot(aln_circos.dotplots, label_seqs=TRUE) + theme +
  ggtitle("Filtered PAF mapping")

# Create dotplot from the raw PAF file.
raw_aln_circos.dotplots = raw_aln_circos
raw_aln_circos.dotplots = reduce_column(raw_aln_circos.dotplots, "qname")
raw_aln_circos.dotplots = reduce_column(raw_aln_circos.dotplots, "tname")
pafr::dotplot(raw_aln_circos.dotplots, label_seqs=TRUE) + theme +
  ggtitle("Raw PAF mapping")
dev.off()

## COVERAGE PLOTS by pafr ##

# Create coverage plots to track mapped and not mapped regions.
pdf(paste0(path_to_paf,".coverage.pdf"))
pafr::plot_coverage(aln_circos, fill="qname") +
  ##ggtitle("Filtered PAF mapping") +
  scale_fill_brewer(palette="Set1")
pafr::plot_coverage(aln_circos, fill="mapq")
  ##ggtitle("Filtered PAF mapping")
pafr::plot_coverage(aln_circos, fill="alen")
  ##ggtitle("Filtered PAF mapping")
pafr::plot_coverage(aln_circos, fill="tname", target=FALSE) +
  ##ggtitle("Filtered PAF mapping") +
  scale_fill_brewer(palette="Set1")
pafr::plot_coverage(aln_circos, fill="mapq", target=FALSE)
  ##ggtitle("Filtered PAF mapping")
pafr::plot_coverage(aln_circos, fill="alen", target=FALSE)
  ##ggtitle("Filtered PAF mapping")
## Plot raw coverage... do we really need it?
## pafr::plot_coverage(raw_aln_circos, fill="qname") +
##   ggtitle("Raw PAF mapping")
dev.off()

## "PARALLEL" SYNTENY PLOT ##

# Plot the syntenic relations of chromosomes X. Start by finding their
# sequids in tname and qname columns:
concat_sequid_columns = c(aln_circos$qname, aln_circos$tname)
unique_sequids = unique(concat_sequid_columns)
# Grab all sequids from unique_sequids that end in X (regex "X$")
sexual_sequids = unique_sequids[grepl("X$", unique_sequids)]
pdf(paste0(path_to_paf,".parallel_sex_chr.pdf"))
plot_synteny(aln_circos,
             q_chrom=sexual_sequids[grepl("^Q", sexual_sequids)],
             t_chrom=sexual_sequids[grepl("^T", sexual_sequids)],
             centre=TRUE) + theme_bw()
dev.off()

## PREPARE CHROMOSOMAL INDICES ##

# Create a list with all unique chromosome IDs. in "qname" column.
query_main_chrs = unique(aln_circos[
    ! grepl("Scaffold", aln_circos$qname), "qname"])
target_main_chrs = unique(aln_circos[
    ! grepl("Scaffold", aln_circos$tname), "tname"])

query_idx = data.frame()
target_idx = data.frame()

for (sequid in query_main_chrs) {
  row = data.frame(
    sequid = sequid,
    start = 0,
    end = aln_circos[aln_circos$qname == sequid, 'qlen'][1]
    )
  query_idx = rbind(query_idx, row)
}

for (sequid in target_main_chrs) {
  row = data.frame(
    sequid = sequid,
    start = 0,
    end = aln_circos[aln_circos$tname == sequid, 'tlen'][1]
    )
  target_idx = rbind(target_idx, row)
}

query_idx = query_idx[order(query_idx$end), ]
target_idx = target_idx[order(target_idx$end), ]
index = rbind(query_idx, target_idx)
# Restart index row numbers so they are sorted from one to N.
row.names(index) = NULL

## INVERTING CHROMOSOMES ##

# Chromosomes are not always in the correct order. The beginning of some
# chromosomes pair with the end of their homologous chromosome. This leads to
# "orthogonal" syntenic relations instead of "parallel" syntenic relations.
# (entangled syntenic relations)
reverse_coord_windows <- function(dfw, inverted_sequids) {
  # Reverse start and end columns of dfw only for selected inverted_sequids
  dfw.rev = dfw
  for (s in inverted_sequids) {
    # Range of the selected (s) chromosome which will be inverted
    xrange = as.double(index[
                       grepl(s, index$sequid),
                       c('start', 'end')])
    # Emmascara el df que hauria d'invertir-se
    # Boolean list which masks by sequid (the one which should be inverted).
    Qmbool = grepl(s, dfw$qname)
    Tmbool = grepl(s, dfw$tname)
    dfw.rev[Qmbool, 'qstart'] = xrange[2] - dfw[Qmbool, 'qend'] + xrange[1]
    dfw.rev[Qmbool, 'qend'] = xrange[2] - dfw[Qmbool, 'qstart'] + xrange[1]
    dfw.rev[Tmbool, 'tstart'] = xrange[2] - dfw[Tmbool, 'tend'] + xrange[1]
    dfw.rev[Tmbool, 'tend'] = xrange[2] - dfw[Tmbool, 'tstart'] + xrange[1]
  }
  return(dfw.rev)}

aln_circos = reverse_coord_windows(aln_circos, inverted_sequids)

## PLOT CIRCOS FROM THE GIVEN PAF ##

# Open a PDF output device for the CIRCOS plot.
pdf(paste0(path_to_paf, ".circos.pdf"))
# Plot two CIRCOS; `aln_circos` and `worst_aln_first`.
for (data_circos in list(aln_circos, worst_aln_first)) {

circos.clear()
amount_sequids_query = length(grep("^Q", index$sequid)) -1
amount_sequids_target = length(grep("^T", index$sequid)) -1
# Specify size of the CIRCOS (amount of chromosomes, spaces between chrs.)
circos.par(gap.degree = c(rep(2, amount_sequids_query), 6,
                          rep(2, amount_sequids_target), 6),
           track.margin = c(mm_h(1), mm_h(1))
           )
# Initialise CIRCOS
circos.genomicInitialize(index, plotType=NULL)

# Create the outer tracks (symbolising chromosomes)
circos.track(index$sequid,
             ylim=c(0,1), track.height = mm_h(1),
             cell.padding=c(0, 0, 0, 0),
             #bg.border = NA, # what did this do?
             panel.fun=function(x, y, ...) {
               # shorthands for a few properties
               cell.chr = CELL_META$sector.index
               xlim = CELL_META$xlim
               ylim = CELL_META$ylim
               # Specify the position of the labels:
               circos.text(CELL_META$xcenter, ylim[2] + mm_y(5),
               # Specify the text in the labels:
               substr(cell.chr, start=nchar(cell.chr), stop=nchar(cell.chr)),
               cex = 1, niceFacing = TRUE)
               # Add a genomic axis above the chromosomes track:
               circos.genomicAxis(h = "top")
               # Add outer rectangles that point towards mappings between chrs
               # and minor scaffolds
               h = outer_highlight[outer_highlight$sequid == cell.chr, ]
               circos.genomicPoints(h[, c("start", "end")], 0.5,
                                    col="red", pch=16, cex=.35, #border=NA
                                    )
             })
# Add the empty rectangles below the chromosomes
circos.track(ylim = c(0, 1), track.height = .1)
# Colour the chromosomes
for (seq in index$sequid) {
  # Remove the labels "Q." or "T." at the beginning of the string
  s = substr(seq, start=3, stop=nchar(seq))
  colour = colour_pairs[
    grepl(s, colour_pairs$sequid), "colour"]
  highlight.chromosome(seq, col=colour, track.index = 2)
}

# Create the syntenic links between homologous chromosomes.
circos.genomicLink(data_circos[, c("qname", "qstart", "qend")],
                   data_circos[, c("tname", "tstart", "tend")],
                   col = data_circos[, "colours"])

# Add text labelling each species' genomes
text(-0.85, 0.9, species_query, col=species_query.col,
     cex = 1.1, font = 4)
text(0.9, -0.85, species_target, col=species_target.col,
     cex = 1.1, font = 4)
circos.clear()
}


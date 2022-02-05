
### .paf reader:
# https://github.com/dwinter/pafr
library("pafr")

## Proves amb el testfile de pafr ####

test_alignment <- system.file("extdata", "fungi.paf", package="pafr")
ali <- read_paf(test_alignment)
dotplot(ali)

# Omplir els testfiles amb les diferents columnes:
plot_coverage(ali, fill='alen') 
plot_coverage(ali, fill='dv')
plot_coverage(ali, fill='strand')
# No és gaire informatiu, perque nmatches és dependent de alen.
plot_coverage(ali, fill='nmatch')
# Nova columna: nmatch/base
ali$nmatch_per_base = ali$nmatch/ali$alen
plot_coverage(ali, fill = 'nmatch_per_base')

# Plotejar la divergència seria... 
# Identitat = nmatches/alen = ali$nmatch_per_base
# nmatches + mismatches == alen
(ali$nmatch + ali$NM) == ali$alen # Es cumpleix la condició.
# No hi ha cap entrada on l'equació anterior no es cumpleixi.

 #

# Es pot intentar veure qualitat de mapatge:
plot_coverage(ali, fill='mapq')
# El nombre de (mismatches+gaps) pateix el mateix que nmatch
plot_coverage(ali, fill='NM')
# S'ha de disassociar de alen, per tant nova fila: NM/alen o NM_per_base:
ali$NM_per_base = ali$NM/ali$alen
ali$NM_per_base
plot_coverage(ali, fill = 'NM_per_base')

## Proves amb el cromosoma X ####
##=============================##

chrx <- read_paf("Bioinformatics/testfiles/aln_dysd_x.paf")
dotplot(chrx)

# distribució (histograma) de la mida dels aliniaments
ggplot(chrx, aes(alen)) + 
  geom_histogram(colour="black", fill="steelblue", bins=20) + 
  theme_bw(base_size=16) + 
  ggtitle("Distribution of alignment lengths (test inclòs al paquet R)") +
  scale_x_log10("Alignment-length")

# Fes un subset dels aliniaments més grans, elimina'n els petits.
long_ali <- subset(chrx, alen > .1e4)
plot_synteny(long_ali, q_chrom="DsilChrX", t_chrom="DcatChrX", centre=TRUE)

# Cobertura: parts que solapen entre genomes
plot_coverage(long_ali, fill='qname') +
  scale_fill_brewer(palette="Set1")

# coverage sense filtrar.
plot_coverage(chrx, fill='qname') +
  scale_fill_brewer(palette="Set1")


## Genome-wide analysis ####
##========================##

# fer dotplot total
genome_wide <- read_paf("Escritorio/filtered_whole_tcat.paf")
dotplot(genome_wide)

# distribució (histo  grama) de la mida dels aliniaments
ggplot(genome_wide, aes(alen)) + 
  geom_histogram(colour="black", fill="steelblue", bins=20) + 
  theme_bw(base_size=16) + 
  ggtitle("Distribution of alignment lengths (test inclòs al paquet R)") +
  scale_x_log10("Alignment-length")

# coverage sense filtrar.
plot_coverage(genome_wide, fill='qname') +
  scale_fill_brewer(palette="Set1")

# coverage del query (silvatica en comptes de catalonica):
plot_coverage(genome_wide, target = FALSE, fill = 'tname') +
  scale_fill_brewer(palette = 'Set1')

# Fes un subset dels aliniaments més grans, elimina'n els petits.
long_ali <- subset(genome_wide, alen > 3e4)
plot_synteny(long_ali, q_chrom="DsilChrX", t_chrom="DcatChrX", centre=TRUE)

# Cobertura: parts que solapen entre genomes
plot_coverage(long_ali, fill='qname') +
  scale_fill_brewer(palette="Set1")

long_ali <- subset(genome_wide, alen > .1e4)
plot_synteny(genome_wide, q_chrom="DsilChr1", t_chrom="DcatChrX", centre=TRUE)

# Fes artificialment els aliniaments més grans, per poder visualitzar millor els colors:
col_genome_wide = genome_wide
x=4000

col_genome_wide$tstart = ifelse(genome_wide$alen<x, genome_wide$tstart - x, genome_wide$tstart)
col_genome_wide$tend = ifelse(genome_wide$alen<x, genome_wide$tend + x, genome_wide$tend)

plot_coverage(col_genome_wide, fill='qname') +
  scale_fill_brewer(palette = "Set1")

# Funció que Modifica llargada d'alen per visualitzar facilment plot_coverage()

artificial_coverage=function(df) {
  # Adds alen to end coord and Substracts alen from start coord.
  # Input: pafreader data.frame, with tstart/tend cols.
  
  # Substracts from tstart
  df$tstart = ifelse( (df$tstart - df$alen)<0, 0, df$tstart - df$alen)
  # Adds to tend
  df$tend = ifelse( (df$tend + df$alen)>df$tlen, df$tlen, df$tend + df$alen)
  # Substracts from qstart
  df$qstart = ifelse( (df$qstart - df$alen)<0, 0, df$qstart - df$alen)
  # Adds to qend
  df$qend = ifelse( (df$qend + df$alen)>df$tlen, df$qlen, df$qend + df$alen)
  return(df)
}
  

## Contra el mateix genoma: ####
#==================================##

# Llegeix silvatica contra silvatica.
sil_v_sil = read_paf("sil_v_sil.paf")
# No es veu res d'interès.
dotplot(sil_v_sil)
# 
plot_coverage(sil_v_sil, fill = 'qname')

## No scaffolding for catalonica ####
#==================================##

noscff = read_paf("Escritorio/noscaff_cat_tcat.paf")
dotplot(noscff)

# distribució (histo  grama) de la mida dels aliniaments
ggplot(noscff, aes(alen)) + 
  geom_histogram(colour="black", fill="steelblue", bins=20) + 
  theme_bw(base_size=16) + 
  ggtitle("Distribution of alignment lengths (noscaff_cat_tcat.paf)") +
  scale_x_log10("Alignment-length")

plot_coverage(noscff, fill='qname')  ## + scale_fill_brewer(palette="Set1")



## circos ####
#===========##

#install.packages("circlize")
library("circlize")

# Prova de links: 
# set seed.
set.seed(123)
# generate random bedfiles (features)
bed1 = generateRandomBed(nr = 100)
bed1 = bed1[sample(nrow(bed1), 20), ]
bed2 = generateRandomBed(nr = 100)
bed2 = bed2[sample(nrow(bed2), 20), ]
# Plot an example circular plot via Ideogram.
circos.initializeWithIdeogram()
# Place links ontop of the plot. 
circos.genomicLink(bed1, bed2, col = rand_color(nrow(bed1), transparency = 0.5), 
                   border = NA)


# CREA EL PLOT DE CIRCOS:
# PER COMENÇAR, CARREGA DATA I LLIBRERIES:
{
  library(pafr)
  library(circlize)
  library(RColorBrewer)
# Anem a usar les dades de Dysdera per formatejar-les com a circos:
aln_circos = read_paf("Bioinformatics/testfiles/filtered_whole_tcat.paf")

# S'ha de filtrar els aliniaments o no es veurà res:
aln_circos = subset(aln_circos, mapq > 20)
aln_circos = subset(aln_circos, alen > 6e3)

# Crea dos bedlike dataframes amb les columnes dels aliniaments:
# Es necessita partir el data.frame en dos de nous.
query = aln_circos[c("qname", "qstart", "qend", "qlen")]
target = aln_circos[c("tname", "tstart", "tend", "tlen")]

# Crea una llista ordenada de bedlike dataframes; 
# cada objecte de la llista són tots els alineaments d'un cromosoma 
# (un sol query o un sol target):
query_list = c() ; target_list = c()  # Dos llistes buides
for (chr in sort(unique(query$qname))) {  # Per a cada qname:
  # Extreu les línies de la query que contenen 'qname==chr'
  query_list = c(query_list, chr = list(query[query$qname == chr,]))
  # Extreu les línies de la target que contenen name==chr'
  target_list = c(target_list, chr = list(target[query$qname == chr,]))
}
# Anomena cada objecte de la llista correctament. 
names(query_list) = sort(unique(query$qname))
names(target_list) = sort(unique(query$qname))

# Reordena la llista amb els aliniaments per estètica:
query_list = c(query_list[9], query_list[1:8])
target_list = c(target_list[9], target_list[1:8])

# index silvatica dels principals chr
index_sil_general = data.frame(
  name = c("DsilChrX", "DsilChr1", "DsilChr2", "DsilChr3", "DsilChr4", "DsilChr5", "DsilChr6", "DsilChrU1", "DsilChrU2"),
  start = c(rep(0,9)),
  end = c(317950935, 177171321, 176727214, 174235172, 129237514, 125974146, 80954097, 22260796, 1135477)
)
# index catalonica dels principals chr
index_cat_general = data.frame(
  name = c("DcatChr2", "DcatChrX", "DcatChr1", "DcatChr3", "DcatChr4"),
  start = c(rep(0,5)),
  end = c(827927493, 662930524, 604492468, 431906129, 421735357)
)
# rbind dels dos data.frames:
dysdera_tracks = rbind(index_cat_general, index_sil_general)

}

# COMENÇA A PLOTEJAR CIRCOS:...
{
circos.clear()  # Per si de cas es necessita netejar la pantalla...

# Separació dels tracks:
circos.par("gap.degree" = c(rep(2, 4), 6, rep(2,8), 6))
circos.genomicInitialize(dysdera_tracks, plotType = NULL)

# Complexa funció per fer uns tracks generals de chr:
circos.track(ylim = c(0,1), panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(1), 
              gsub(".*Chr", "", CELL_META$sector.index), cex = 0.6, niceFacing = TRUE)
  # Afageix axis genòmics a sobre del track [INHABILITAT]:
  # circos.genomicAxis(h = "top")
}, track.height = mm_h(1), cell.padding = c(0, 0, 0, 0), bg.border = NA)

# Crea dues cintes, vermella i blava, per a diferenciar les dues espècies:
highlight.chromosome(unique(query$qname),
                     col = rgb(1,0,0,.7))
highlight.chromosome(unique(target$tname),
                     col = rgb(0,0,1,.7))

# Fés els requadres buits de sota els cromosomes.
circos.track(ylim = c(0, 1), track.height = .15)

# Crea una paleta de colors per a pintar el cercle:
palette_cols = brewer.pal(9, "Set1")

# Posa una paleta de colors sobre els cromosomes de silvatica:
x=1
for (chr in names(query_list)) {
  highlight.chromosome(chr,
  col = palette_cols[x],
  track.index = 2)
  x = x+1
}

# Afageix transparencia a la palette pels links [NO FUNCIONA]:
palette_cols = add_transparency(col = palette_cols, transparency = 0.5)

# Links entre les regions aliniades:
# x=1
# for (n in c("x", 1:6, "U1", "U2")) {
#   print(paste0("target.chr", n, " - ", "query.chr", n))
#   circos.genomicLink(eval(as.name(paste0("target.chr", n))), eval(as.name(paste0("query.chr", n))), 
#                      col = palette_cols[x])
#   x = x+1
# }
for (n in 1:length(query_list)) {
  circos.genomicLink(target_list[n], query_list[n],
                     col = palette_cols[n])
}

circos.clear()
text(-0.85, -0.8, "D. catalonica\ngenome", col = "blue", cex = .6)
text(0.9, 0.8, "D. silvatica\ngenome", col = "red", cex = .6)
}

# COMENÇA A PLOTEJAR COVERAGE:
plot_coverage(aln_circos)
plot_coverage(artificial_coverage(aln_circos))



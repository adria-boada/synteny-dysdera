
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
### Creates circos and coverage plots from minimap2 data (i.e., alignment.paf).

# PER COMENÇAR, CARREGA DADES I LLIBRERIES:
{
  library(pafr)
  library(circlize)
  library(RColorBrewer)
 
# Carrega el camí fins als fitxers necessaris i defineix variables:
  
# .paf amb alineament:
  path_to_paf = "Escritorio/Data/Synteny/modified_pafs/filtered_whole_tcat.paf"
  
# Variables per a filtrar el .paf i reduir-ne la mida:
  min_map_qual = 20  # No es mira per sota de 20 de qualitat de mapatge.
  min_align_len = 6e3  # No es mira alineaments que siguin més petits de 6 kb.
  
# Paleta de colors extreta del paper de Paula Escuer:
  paleta_escuer = c(chrx="#f659a5", chr1="#c5090a", chr2="#ff7f07", 
                    chr3="#cabf0a", chr4="#41a62a", chr5="#4fa1ca", 
                    chr6="#4f3b97", chrU1="#a29a9c", chrU2="#a65628")
  # Test amb un barplot per veure que els colors rutllen:
  if (FALSE) {
  barplot(c(chrx=5, chr1=5, chr2=5, chr3=5, chr4=5, chr5=5, chr6=5, chrU1=5, chrU2=5), col=paleta_escuer)
  }
 # Defineix una funció per augmentar la mida dels aliniaments de la funció plot_coverage()

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
 
# Anem a usar les dades de Dysdera de minimap2 per formatejar-les com a circos:
  aln_circos = read_paf(path_to_paf)

# S'ha de filtrar els aliniaments o no es veurà res:
  aln_circos = subset(aln_circos, mapq > min_map_qual)  # Elimina poca qualitat de mapeig
  aln_circos = subset(aln_circos, alen > min_align_len)  # Elimina alineaments curts
 
# Crea dos bedlike dataframes amb les columnes dels aliniaments:
# Es necessita partir el data.frame en dos de nous per fer els links del circos.
  query = aln_circos[c("qname", "qstart", "qend", "qlen")]
  target = aln_circos[c("tname", "tstart", "tend", "tlen")]

# Crea una llista ordenada de bedlike dataframes; 
# cada objecte de la llista són tots els alineaments d'un cromosoma 
# (un sol query o un sol target):
 
query_list = c() ; target_list = c()  # Dos llistes buides
 
for (chr in sort(unique(query$qname))) {  # Per a cada qname:
  # Extreu les línies de la query que contenen 'qname==chr'
  query_list = c(query_list, chr = list(query[query$qname == chr,]))
  # Extreu les línies de la target que contenen 'tname==chr'
  target_list = c(target_list, chr = list(target[query$qname == chr,]))
}
 
# Anomena cada objecte de la llista correctament. 
names(query_list) = sort(unique(query$qname))
names(target_list) = sort(unique(query$qname))

# Reordena la llista amb els aliniaments per estètica:
  # Primer 'X' i després la resta.
query_list = c(query_list[9], query_list[1:8])
target_list = c(target_list[9], target_list[1:8])
 
# Crea uns indexos que separen el cercle en sectors
# Cada sector és un cromosoma amb la mida que indica l'index.

# index silvatica dels principals chr
index_sil_general = data.frame(
  name = c(unique(query$qname)),
  start = c(rep(0, length(unique(query$qname)))),
  end = c(unique(query$qlen))
)
# index catalonica dels principals chr
index_cat_general = data.frame(
  name = c(unique(target$tname)),
  start = c(rep(0, length(unique(target$tname)))),
  end = c(unique(target$tlen))
)
# Ordena els cromosomes de l'index numericament per llargada (end). 
  index_sil_general = index_sil_general[order(-index_sil_general$end), ]
  index_cat_general = index_cat_general[order(-index_cat_general$end), ]
 
# rbind dels dos data.frames:
dysdera_tracks = rbind(index_cat_general, index_sil_general)
}


# HA ACABAT DE CARREGAR DADES

# COMENÇA A PLOTEJAR CIRCOS:...


{
circos.clear()  # Per si de cas es necessita netejar la pantalla...
 # (en cas de que es corri l'script més d'una vegada)

# Separació dels tracks:
circos.par("gap.degree" = c(rep(2, 4), 6, rep(2,7), 4, 6))
 
# Inicialitza amb les dades de Dysdera.
circos.genomicInitialize(dysdera_tracks, plotType = NULL)

# Complexa funció per fer uns tracks generals de chr:
circos.track(ylim = c(0,1), panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(5), 
              gsub(".*Chr", "", CELL_META$sector.index), cex = 1, niceFacing = TRUE)
  # Afageix axis genòmics a sobre del track [INHABILITAT]:
  circos.genomicAxis(h = "top")
}, track.height = mm_h(1), cell.padding = c(0, 0, 0, 0), bg.border = NA)

# Crea dues cintes, vermella i blava, per a diferenciar les dues espècies:
highlight.chromosome(unique(query$qname),
                     col = rgb(1,0,0,.7))
highlight.chromosome(unique(target$tname),
                     col = rgb(0,0,1,.7))
 
# Fés els requadres buits de sota els cromosomes.
circos.track(ylim = c(0, 1), track.height = .1)

# Crea una paleta de colors per a pintar el cercle [OBSOLET]:
# Es necessitava RColorBrewer.
# palette_cols = brewer.pal(9, "Set1")

# Posa una paleta de colors sobre els cromosomes de silvatica:
x=1
for (chr in names(query_list)) {
  highlight.chromosome(chr,
  col = paleta_escuer[x],
  track.index = 2)
  x = x+1
}

# Afageix transparencia a la palette pels links:
  ### [NO FUNCIONA]
paleta_escuer_transp = add_transparency(col = paleta_escuer, transparency = 0.5)

# Links entre les regions aliniades:
# x=1
# for (n in c("x", 1:6, "U1", "U2")) {
#   print(paste0("target.chr", n, " - ", "query.chr", n))
#   circos.genomicLink(eval(as.name(paste0("target.chr", n))), eval(as.name(paste0("query.chr", n))), 
#                      col = paleta_escuer[x])
#   x = x+1
# }
for (n in 1:length(query_list)) {
  circos.genomicLink(target_list[n], query_list[n],
                     col = paleta_escuer_transp[n])
}

circos.clear()
 
# Text que marca quina espècie és quina: 
text(-0.85, -0.9, "D. catalonica\ngenome", col = "blue", cex = 1.1)
text(0.9, 0.9, "D. silvatica\ngenome", col = "red", cex = 1.1)
}


# HA ACABAT LA FIGURA CIRCOS

# FINALMENT, GENERA ELS "PLOT_COVERAGE()" DE PAFR:
plot_coverage(aln_circos)
plot_coverage(artificial_coverage(aln_circos))



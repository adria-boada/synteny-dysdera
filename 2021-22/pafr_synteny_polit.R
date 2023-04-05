
### PAFR_SYNTENY.R
### Creates `circos` and coverage plots from minimap2 data (i.e., alignment.paf).


## CARREGA DADES, OPCIONS I LLIBRERIES ##
{
# lectura i adequament de fitxers PAF (mapatges del minimap2)
  library(pafr)
# crear gràfiques circulars partint de fitxers/data.frames BED-like
  library(circlize)
  
# Carrega el camí fins als fitxers necessaris i defineix variables:
# PAF amb alineament:
  path_to_paf = "Escritorio/Data/Synteny/modified_pafs/filtered_whole_tcat.paf"
  
# Variables per a filtrar el PAF i reduir-ne la mida:
  min_map_qual = 20  # No es mira per sota d'aquesta qualitat de mapatge.
  # El rang de qualitats de mapatge vira entre zero i seixanta. 
  min_align_len = 6e3  # No es mira alineaments que siguin més petits de 6 kb.
  
## PALETA DE COLORS ##
  # 8 colors més l'afegit U2. 
  paleta_escuer = c(chrx="#f659a5", chr1="#c5090a", chr2="#ff7f07", 
                    chr3="#cabf0a", chr4="#41a62a", chr5="#4fa1ca", 
                    chr6="#4f3b97", chrU1="#a29a9c", chrU2="#a65628")
  # Test amb un barplot per veure que els colors rutllen i són els adequats:
  if (FALSE) {
  barplot(c(chrx=5, chr1=5, chr2=5, chr3=5, chr4=5, chr5=5, chr6=5, chrU1=5, chrU2=5), col=paleta_escuer)
  }

# Carrega el mapatge PAF de minimap2:
  aln_circos = read_paf(path_to_paf)

# Filtrar els aliniaments per adequar el fitxer entrada a les teves necessitats.
# Excés d'aliniaments petits dificulten la computació i interpretació.
  aln_circos = subset(aln_circos, mapq >= min_map_qual)  # Elimina poca qualitat de mapeig
  aln_circos = subset(aln_circos, alen > min_align_len)  # Elimina mapatges curts
  aln_circos = filter_secondary_alignments(aln_circos)   # Elimina mapatges secundaris

# Crea dos bedlike dataframes amb les columnes dels aliniaments:
# Es necessita partir el data.frame en dos de nous per fer els links del circos.
  query = aln_circos[c("qname", "qstart", "qend", "qlen")]
  target = aln_circos[c("tname", "tstart", "tend", "tlen")]

# Crea una llista ordenada de bedlike dataframes; 
# cada objecte de la llista són tots els alineaments d'un cromosoma 
# (un sol query o un sol target):
query_list = c() ; target_list = c()  # Dos llistes buides
# Per a cada cromosoma query (qname):
for (chr in sort(unique(query$qname))) {
  # Extreu les línies de la query que contenen 'qname==chr'
  query_list = c(query_list, chr = list(query[query$qname == chr,]))
  # Extreu les línies de la target que contenen 'tname==chr'
  target_list = c(target_list, chr = list(target[query$qname == chr,]))
}
# Anomena cada objecte de la llista correctament. 
names(query_list) = sort(unique(query$qname))
names(target_list) = sort(unique(target$tname))
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

## GRAFICAR CIRCOS DEL PAF INPUT ##

{
circos.clear()  # Per si de cas, s'eliminen les últimes opcions de circos abans de res...
 # (en cas de que es corri l'script interactivament i més d'una vegada)

# Separació dels tracks (gaps entre chr): dos entre cromosomes de la mateixa
# espècie (repetit 4 i 7 vegades) i 6 entre genomes de diferents espècies.
circos.par("gap.degree" = c(rep(2, 4), 6, rep(2,7), 4, 6))

# Inicialitza amb les dades de Dysdera.
circos.genomicInitialize(dysdera_tracks, plotType = NULL)

# Complexa funció per fer uns tracks generals de chr:
circos.track(ylim = c(0,1), panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(5), 
             gsub(".*Chr", "", CELL_META$sector.index), cex = 1, niceFacing = TRUE)
  # Afegeix axis genòmics a sobre del track:
  circos.genomicAxis(h = "top")
}, 
track.height = mm_h(1), cell.padding = c(0, 0, 0, 0), bg.border = NA)

# Crea dues cintes, vermella i blava, per a diferenciar les dues espècies:
highlight.chromosome(unique(query$qname),
                     col = rgb(1,0,0,.7))
highlight.chromosome(unique(target$tname),
                     col = rgb(0,0,1,.7))

# Fés els requadres buits de sota els cromosomes.
circos.track(ylim = c(0, 1), track.height = .1)

# Posa una paleta de colors sobre els cromosomes de silvatica:
x=1
for (chr in names(query_list)) {
  highlight.chromosome(chr,
  col = paleta_escuer[x],
  track.index = 2)
  x = x+1
}

# Afegeix transparència a la palette pels links:
### [NO FUNCIONA per excés de dades a plotejar.]
### Es solapen tants objectes semi-transparents que el color
### final acaba sent igualment sòlid.
paleta_escuer_transp = add_transparency(col = paleta_escuer, transparency = 0.5)

# Links entre les regions aliniades:
for (n in 1:length(query_list)) {
  circos.genomicLink(target_list[n], query_list[n],
                     col = paleta_escuer_transp[n])
}

circos.clear()

# Text que marca quina espècie és quina: 
text(-0.85, -0.9, "D. catalonica\ngenome", col = "blue", cex = 1.1)
text(0.9, 0.9, "D. silvatica\ngenome", col = "red", cex = 1.1)
}


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

## ANNEX: GRÀFIQUES COMPLEMENTÀRIES ##

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


## FIGURES CIRCOS AMB LES CONNEXIONS D'UN SOL CROMOSOMA ##
# Amb la intenció que no es solapin unions i es pugui valorar
# les unions de cromosomes individuals d'interès.
if (FALSE) {  # Desactivar secció de manera predeterminada

# Directori on enviar múltiples figures individualment.
setwd("Escritorio/Results/Synteny/circos/")

# Crea una figura per a cada conjunt de links:
for (link in 1:length(query_list)){
  
  # Crea les figures en format .pdf:
  pdf(file = paste0("individual_circos_", link, ".pdf"), title = "Individual genomic links circos")
  
  # Per si de cas, s'eliminen les últimes opcions de circos abans de res...
  circos.clear()  # (en cas de que es repeteixi l'script)
  
  # Separació dels tracks (gaps entre chr):
  circos.par("gap.degree" = c(rep(2, 4), 6, rep(2,7), 4, 6))
  # Inicialitza amb les dades de Dysdera.
  circos.genomicInitialize(dysdera_tracks, plotType = NULL)

  # Complexa funció per fer uns tracks generals de chr:
  circos.track(ylim = c(0,1), panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(5),
                gsub(".*Chr", "", CELL_META$sector.index), cex = 1, niceFacing = TRUE)
    # Afageix axis genòmics a sobre del track [INHABILITAT?]:
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
  # Actualment usa paleta_escuer. Necessitava RColorBrewer.
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
  ### [NO FUNCIONA per excés de dades a plotejar.]
  paleta_escuer_transp = add_transparency(col = paleta_escuer, transparency = 0.5)
 
  # Links entre les regions aliniades
  # (a partir de la variable link de dalt que engloba aquesta secció)
   {
    circos.genomicLink(target_list[link], query_list[link],
                       col = paleta_escuer_transp[link])
  }
  
  circos.clear()
  
  # Text que marca quina espècie és quina: 
  text(-0.85, -0.9, "D. catalonica\ngenome", col = "blue", cex = 1.1)
  text(0.9, 0.9, "D. silvatica\ngenome", col = "red", cex = 1.1)

  # Tanca device pdf per crear la pròxima figura:
  dev.off()
}
}


## "PLOT_COVERAGE()" DE PAFR ##
if (FALSE) { # Desactivar aquesta secció de manera predeterminada

# Defineix una funció per augmentar la mida dels aliniaments de la funció plot_coverage()
# Amb cromosomes enormes i mapatges curts (comparativament)
# els colors queden difuminats i el gràfic és difícil d'interpretar.
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
 
# No permet personalitzar la paleta de colors 
plot_coverage(aln_circos)
plot_coverage(artificial_coverage(aln_circos))
}


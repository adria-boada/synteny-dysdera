
### PAFR_SYNTENY.R
### Creates `circos` and coverage plots from minimap2 data (i.e., alignment.paf).

# Carrega el camí fins als fitxers necessaris i defineix variables:
# PAF amb alineament:
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
  stop("Supply a PAF file as the first arg. value")
}

## CARREGA DADES, OPCIONS I LLIBRERIES ##

# lectura i adequament de fitxers PAF (mapatges del minimap2)
  library(pafr)
# crear gràfiques circulars partint de fitxers/data.frames BED-like
  library(circlize)
  
  path_to_paf = args[1]
  
# Variables per a filtrar el PAF i reduir-ne la mida:
  min_map_qual = 30  # No es mira per sota d'aquesta qualitat de mapatge.
  # El rang de qualitats de mapatge vira entre zero i seixanta. 
  min_align_len = 6e3  # No es mira alineaments que siguin més petits de 6 kb.
  
## PALETA DE COLORS ESCUER--SILVATICA ##

  # 8 colors més l'afegit U2. 
  paleta_escuer = c(DsilChrX="#f659a5", DsilChr1="#c5090a", DsilChr2="#ff7f07",
                    DsilChr3="#cabf0a", DsilChr4="#41a62a", DsilChr5="#4fa1ca",
                    DsilChr6="#4f3b97", DsilChrU1="#a29a9c",
                    DsilChrU2="#a65628")
  # Test amb un barplot per veure que els colors rutllen i són els adequats:
  if (FALSE) {
  barplot(c(DsilChrX=5, DsilChr1=5, DsilChr2=5, DsilChr3=5, DsilChr4=5, DsilChr5=5, DsilChr6=5, DsilChrU1=5, DsilChrU2=5), col=paleta_escuer)
  }

## CARREGA EL MAPATGE (UNIONS ENTRE CHR) ##

# Cada mapatge es representa unint la parella de
# coordenades dels dos cromosomes al gràfic CIRCOS.

# Carrega el mapatge PAF de minimap2:
  raw_aln_circos = read_paf(path_to_paf)

# Etiqueta cromosomes amb 'Q.' i 'T.' al principi
  raw_aln_circos$tname = paste0('T.', raw_aln_circos$tname)
  raw_aln_circos$qname = paste0('Q.', raw_aln_circos$qname)

# Filtrar els aliniaments per adequar el fitxer entrada a les teves necessitats.
# Excés d'aliniaments petits dificulten la computació i interpretació.

  # M'agradaria començar per un filtre un pèl complex;
  # percentatge de vincles entre cromosomes no homòlegs
  homologous = nrow(raw_aln_circos[
      # X versus X
      (grepl('cat.*X', raw_aln_circos$tname) & grepl('til.*X', raw_aln_circos$qname)) |
      (grepl('til.*X', raw_aln_circos$tname) & grepl('cat.*X', raw_aln_circos$qname)) |
      # til7 and til6 vs cat1
      (grepl('til.*5', raw_aln_circos$tname) & grepl('cat.*2', raw_aln_circos$qname)) |
      (grepl('til.*6', raw_aln_circos$tname) & grepl('cat.*2', raw_aln_circos$qname)) |
      (grepl('cat.*2', raw_aln_circos$tname) & grepl('til.*6', raw_aln_circos$qname)) |
      (grepl('cat.*2', raw_aln_circos$tname) & grepl('til.*5', raw_aln_circos$qname)) |
      # til3 and til5 vs cat2
      (grepl('til.*1', raw_aln_circos$tname) & grepl('cat.*1', raw_aln_circos$qname)) |
      (grepl('til.*4', raw_aln_circos$tname) & grepl('cat.*1', raw_aln_circos$qname)) |
      (grepl('cat.*1', raw_aln_circos$tname) & grepl('til.*4', raw_aln_circos$qname)) |
      (grepl('cat.*1', raw_aln_circos$tname) & grepl('til.*1', raw_aln_circos$qname)) |
      # til4 vs cat3
      (grepl('cat.*3', raw_aln_circos$tname) & grepl('til.*3', raw_aln_circos$qname)) |
      (grepl('til.*3', raw_aln_circos$tname) & grepl('cat.*3', raw_aln_circos$qname)) |
      # til 2 vs cat 4
      (grepl('cat.*4', raw_aln_circos$tname) & grepl('til.*2', raw_aln_circos$qname)) |
      (grepl('til.*2', raw_aln_circos$tname) & grepl('cat.*4', raw_aln_circos$qname))
    , ])

  t = nrow(raw_aln_circos)
  df_double_minor_scaffold = raw_aln_circos[
                                grepl("Scaffold|ctg", raw_aln_circos$qname) &
                                grepl("Scaffold|ctg", raw_aln_circos$tname), ]
  removed.both_min_scaffold = nrow(df_double_minor_scaffold)
  df_minor_scaffold = raw_aln_circos[
                                grepl("Scaffold|ctg", raw_aln_circos$qname) |
                                grepl("Scaffold|ctg", raw_aln_circos$tname), ]
  df_single_scaffold = df_minor_scaffold[
                                grepl('chr', df_minor_scaffold$qname) |
                                grepl('chr', df_minor_scaffold$tname), ]
  removed.single_min_scaffold = nrow(df_single_scaffold)
  aln_circos = raw_aln_circos[ # sequid names do not (!) contain "Scaffold"
                              # (at least one of these names is a chromosome)
                                ! grepl("Scaffold|ctg", raw_aln_circos$qname) &
                                ! grepl("Scaffold|ctg", raw_aln_circos$tname), ]
  # Elimina poca qualitat de mapeig
  aln_circos = subset(aln_circos, mapq >= min_map_qual)
  df_single_scaffold = subset(df_single_scaffold, mapq >= (min_map_qual-5))
  removed.mapq = nrow(raw_aln_circos)-nrow(subset(raw_aln_circos, mapq >= min_map_qual))
  # Elimina mapatges curts
  aln_circos = subset(aln_circos, alen > min_align_len)
  df_single_scaffold = subset(df_single_scaffold, alen > (min_align_len-3500))
  removed.alilen = nrow(raw_aln_circos)-nrow(subset(raw_aln_circos, alen > min_align_len))
  # Elimina mapatges secundaris
  aln_circos = filter_secondary_alignments(aln_circos)
  removed.secali = nrow(raw_aln_circos)-nrow(filter_secondary_alignments(raw_aln_circos))
  # Total ajuntant totes les eliminacions... (tenint en compte solapaments)
  removed.all = nrow(raw_aln_circos)-nrow(aln_circos)

  removed.single_min_scaffold.qualfiltered = nrow(df_single_scaffold)

  # imprimeix valors en pantalla
  cat(paste0('Complete/raw amount of rows: ', nrow(raw_aln_circos), '\n'))
  cat(paste0('Rows mapping pairs of homologous/contiguous chromosomes: ',
             round((homologous/t)*100, digits=1), ' %\n'))
  cat(paste0('Rows mapping one minor scaffold to another minor scaffold (unused): ',
             round((removed.both_min_scaffold/t)*100, digits=1), ' %\n'))
  cat(paste0('Rows mapping one minor scaffold to a main chromosome (unused): ',
             round((removed.single_min_scaffold/t)*100, digits=1), ' %\n'))
  cat(paste0('Rows mapping one minor scaffold to a main chromosome WITH HIGH QUAL. (outer track): ',
             round((removed.single_min_scaffold.qualfiltered/t)*100, digits=1), ' %\n'))
  cat(paste0('Rows removed because their mapQ. was lower than ',
             min_map_qual, ': ', round((removed.mapq/t)*100, digits=1), ' %\n'))
  cat(paste0('Rows removed because their ali. length was lower than ',
             min_align_len, ': ', round((removed.alilen/t)*100, digits=1), ' %\n'))
  cat(paste0('Rows removed because they were secondary mappings: ',
             round((removed.secali/t)*100,digits=1), ' %\n'))
  cat(paste0('The combined percentage of original rows that have been removed is: ',
             round((removed.all/t)*100,digits=1), ' %\n'))

# Emparella colors i cromosomes al df `aln_circos`
  colors_correspondencia = data.frame(
    sequid = c("DtilchrX", "Dtilchr1", "Dtilchr2",
               "Dtilchr3", "Dtilchr4", "Dtilchr5",
               "Dtilchr6"),
    colors = c("#ff6dd3", "#ff0008", "#5bc358",
               "#ffbb54", "#ff7777", "#55e8f0",
               "#b9dae9"))

  # inicialitza columna amb colors
  aln_circos$colors <- NA
  for (i in 1:nrow(colors_correspondencia)) {
    color = colors_correspondencia[i, 'colors']
    sequid = colors_correspondencia[i, 'sequid']
    # omple la columna de colors
    aln_circos[
          grepl(sequid, aln_circos$qname) | # omple fileres on qname == sequid
          grepl(sequid, aln_circos$tname),  # (or tname)
      'colors'] = color                     # omple columna 'colors' amb la variable color
  }

# Ordena els aliniaments segons mapq;
  # d'aquesta manera, els pitjors mapatges quedaràn sotarrats pels millors,
  # i no al revés. Els millors al final (els últims a ser dibuixats)
  aln_circos = aln_circos[order(aln_circos$mapq, aln_circos$alen), ]

  # snippet of the filtered dataframe # debug
  cat('\n')
  #head(as.data.frame(aln_circos[, -24]))
  cat('\n')

# Busca les regions de sequids principals que troben mapatges a scaffolds
  # minoritaris
  query_ressaltat = df_single_scaffold[grepl('Scaffold|ctg', df_single_scaffold$tname),
                      c('qname', 'qstart', 'qend', 'alen', 'mapq')]
  names(query_ressaltat) = c('sequid', 'start', 'end', 'alen', 'mapq')
  target_ressaltat = df_single_scaffold[grepl('Scaffold|ctg', df_single_scaffold$qname),
                      c('tname', 'tstart', 'tend', 'alen', 'mapq')]
  names(target_ressaltat) = c('sequid', 'start', 'end', 'alen', 'mapq')
  ressaltat_exterior = rbind(query_ressaltat, target_ressaltat)

  print('## HEAD() DEL DATAFRAME PER A RESSALTAR EXTERIOR ##')
  head(as.data.frame(ressaltat_exterior))

## CARREGA ÍNDEX DELS CROMOSOMES ##

# Crea uns indexos que separen el cercle en sectors
# Cada sector és un cromosoma amb la mida que indica l'index.

  # query index (remove qnames labeled as 'Scaffold')
main_query_chr = aln_circos[! grepl('Scaffold', aln_circos$qname), 'qname']
main_query_chr = unique(main_query_chr) # unique list of main sequids
main_target_chr = aln_circos[! grepl('Scaffold', aln_circos$tname), 'tname']
main_target_chr = unique(main_target_chr) # unique list of main sequids

index_query = data.frame()
index_target = data.frame()

for (sequid in main_query_chr) {
  row = data.frame(
    sequid = sequid,
    start = 0,
    end = aln_circos[aln_circos$qname == sequid, 'qlen'][1]
    )
  index_query = rbind(index_query, row)
}

for (sequid in main_target_chr) {
  row = data.frame(
    sequid = sequid,
    start = 0,
    end = aln_circos[aln_circos$tname == sequid, 'tlen'][1]
    )
  index_target = rbind(index_target, row)
}

index_query = index_query[order(index_query$end), ]
index_target = index_target[order(index_target$end), ]
dysdera_tracks = rbind(index_query, index_target)

cat('\n')
print('## DYSDERA TRACKS ##')
print(dysdera_tracks)
cat('\n')


### INTENT INVERTIR CROMOSOMES ###
# Els cromosomes no es troben en l'ordre correcte. Hi ha inicis que troben
# el seu homòleg al final del cromosoma d'una altra espècie. D'aquesta manera,
# s'han de capgirar manualment si es disposa de grans blocs colinears que
# ho indiquin...
inverted_sequids = c(
  'DtilchrX',
  'Dtilchr1',
  'Dcatchr2'
)

reverse_coord_windows <- function(dfw, inverted_sequids) {
  # reverse start and end columns of dfw only for selected inverted_sequids
  dfw.rev = dfw
  for (s in inverted_sequids) {
    # range of the selected (s) chromosome which will be inverted
    xrange = as.double(dysdera_tracks[
                       grepl(s, dysdera_tracks$sequid),
                       c('start', 'end')])
    # emmascara el df que hauria d'invertir-se
    Qmbool = grepl(s, dfw$qname)
    Tmbool = grepl(s, dfw$tname)
    dfw.rev[Qmbool, 'qstart'] = xrange[2] - dfw[Qmbool, 'qend'] + xrange[1]
    dfw.rev[Qmbool, 'qend'] = xrange[2] - dfw[Qmbool, 'qstart'] + xrange[1]
    dfw.rev[Tmbool, 'tstart'] = xrange[2] - dfw[Tmbool, 'tend'] + xrange[1]
    dfw.rev[Tmbool, 'tend'] = xrange[2] - dfw[Tmbool, 'tstart'] + xrange[1]
  }
  return(dfw.rev)}

aln_circos = reverse_coord_windows(aln_circos, inverted_sequids)


## INFORMACIÓ RESUM ##

########################## Obre output pdf
pdf(paste0('genplot_', args[1], '.pdf'),
    width=8)
##########################

species_query = dysdera_tracks[grepl('^Q', dysdera_tracks$sequid), 'sequid'][1]
species_query = substr(species_query, start=3, stop=6)
species_target = dysdera_tracks[grepl('^T', dysdera_tracks$sequid), 'sequid'][1]
species_target = substr(species_target, start=3, stop=6)

plot.new()
legend('center', title='General summary information', cex=0.9,
  legend = c(
  paste0('(1) Query: ', species_query, '; (2) Target: ', species_target, '; (3) Total rows: ', nrow(raw_aln_circos)),
  paste0('+ Rows mapping pairs of homologous chromosomes: ',
             round((homologous/t)*100, digits=1), ' %'),
  paste0('+ Rows mapping one minor scaffold to another\nminor scaffold (discarded from plot): ',
             round((removed.both_min_scaffold/t)*100, digits=1), ' %'),
  paste0('+ Rows mapping one minor scaffold to a main\nchromosome (discarded from plot): ',
             round((removed.single_min_scaffold/t)*100, digits=1), ' %'),
  paste0('+ Rows mapping one minor scaffold to a main chromosome\nWITH HIGH QUAL.',
         ' (outmost red/white track,\nbeside coordinate axis): ',
             round((removed.single_min_scaffold.qualfiltered/t)*100, digits=1), ' %'),
  paste0('+ Rows removed because their mapQ. was lower than ',
             min_map_qual, ': ', round((removed.mapq/t)*100, digits=1), ' %'),
  paste0('+ Rows removed because their ali. length was lower than ',
             min_align_len, ': ', round((removed.alilen/t)*100, digits=1), ' %'),
  paste0('+ Rows removed because they were secondary mappings: ',
             round((removed.secali/t)*100,digits=1), ' %\n'),
  paste0('+ The combined percentage of original rows\nthat have been removed is: ',
             round((removed.all/t)*100,digits=1), ' %')
))


## GRAFICAR CIRCOS DEL PAF INPUT ##

circos.clear()  # Per si de cas, s'eliminen les últimes opcions de circos abans de res...
 # (en cas de que es corri l'script interactivament i més d'una vegada)

# Separació dels tracks (gaps entre chr): 2 "unitats" entre cromosomes de la
# mateixa espècie i 6 "unitats" entre genomes de diferents espècies.
query_sequids = length(grep('Q.', dysdera_tracks$sequid)) -1
target_sequids = length(grep('T.', dysdera_tracks$sequid)) -1
circos.par(gap.degree = c(rep(2, query_sequids), 6,
                            rep(2, target_sequids), 6),
           # ajunta primera i segona via/tracks
           track.margin = c(mm_h(1), mm_h(1))
)

# Inicialitza amb les dades de Dysdera.
circos.genomicInitialize(dysdera_tracks,
                         plotType = NULL)

# Complexa funció per fer uns tracks generals de chr:
circos.track(dysdera_tracks$sequid,
  ylim = c(0,1), track.height = mm_h(1),
  cell.padding = c(0, 0, 0, 0), #bg.border = NA,
  panel.fun = function(x, y, ...) {
    # shorthands for a few sector properties
    cllcrm = CELL_META$sector.index
    xlim   = CELL_META$xlim
    ylim   = CELL_META$ylim
    # text que etiqueta als cromosomes...
    circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(5), 
             gsub(".*chr", "", cllcrm), cex = 1, niceFacing = TRUE)
    # Afegeix axis genòmics a sobre del track:
    circos.genomicAxis(h = "top")
    # rectangles que assenyalen mapatges entre chr-scaff.
    r = ressaltat_exterior[ressaltat_exterior$sequid == cllcrm, ]
    circos.genomicPoints(r[, c('start', 'end')], 0.5,
              col='red', pch=16, cex=0.35, #border = NA
    )
})

# Fés els requadres buits de sota els cromosomes.
circos.track(ylim = c(0, 1), track.height = .1)

# Posa una paleta de colors sobre els cromosomes seleccionats
# manualment dins el data.frame `colors_correspondencia`...
for (seq in dysdera_tracks$sequid) {
  # elimina etiqueta Q o T
  s = substr(seq, start=3, stop=nchar(seq))
  color = colors_correspondencia[
                       grepl(s, colors_correspondencia$sequid),
                       'colors']
  highlight.chromosome(seq, col=color, track.index = 2
  )
}

# Afegeix transparència a la palette pels links:
### [NO FUNCIONA per excés de dades a plotejar.]
### Es solapen tants objectes semi-transparents que el color
### final acaba sent igualment sòlid.
### paleta_escuer_transp = add_transparency(col = paleta_escuer, transparency = 0.5)

circos.genomicLink(aln_circos[, c('qname', 'qstart', 'qend')],
                   aln_circos[, c('tname', 'tstart', 'tend')],
                   col = aln_circos[, 'colors'])

# Text que marca la posició de les dues espècies: 
  # Si tilosensis és query, situar el seu nom a la part inferior dreta
  if (any(grepl('^Q.*til', dysdera_tracks$sequid))) {
    text(-0.85, 0.9, "D. catalonica\ngenome", col = "#eb8f46", cex = 1.1, font=2)
    text(0.9, -0.9, "D. tilosensis\ngenome", col = "#30acac", cex = 1.1, font=2)
  # Si tilosensis és target, situar el seu nom a la part superior dreta
  } else if (any(grepl('^T.*til', dysdera_tracks$sequid))) {
    text(-0.85, -0.9, "D. catalonica\ngenome", col = "#eb8f46", cex = 1.1, font=2)
    text(0.9, 0.9, "D. tilosensis\ngenome", col = "#30acac", cex = 1.1, font=2)
}

circos.clear()


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
for (link in 1:length(bed_query)){
  
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
  }, track.height = mm_h(1), cell.padding = c(1, 1, 1, 1), bg.border = NA)
  
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
  for (chr in names(bed_query)) {
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
    circos.genomicLink(bed_target[link], bed_query[link],
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


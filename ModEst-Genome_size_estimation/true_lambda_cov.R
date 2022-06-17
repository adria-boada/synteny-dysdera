
library(fitdistrplus) 
library(truncdist) 
library(splitstackshape) 

# Create a function to find the median without placing values in a vector
# (too long, no mem available)
median_zeroidx_cov <- function(x) {
  # Si la suma del vector és parella:
  if ((sum(x) %% 2) == 0) {
    # el valor medià és una mitjana de dos valors
    med_val = c((sum(x)/2)+1, (sum(x)/2)-1)
    solution = 0
    # Per als valors +/- 1:
    for (i in med_val) {
      a = 0
      line = 0
      # A quins valors discrets pertany?
      while (a < i) {line = line+1 ; a = a+sum(x[line])}
      solution = solution + line
    }
    # Com que és zero-indexed, resta dos a sol.
    return((solution-2)/2)
  }
  
  # Si la suma és senar:
  else {
    med_val = ( sum(x)+1 )/2
    a = 0
    line = 0
    while (a < med_val) {line = line+1 ; a = a+sum(x[line])}
    # Com que el cov és zero-idx, resta u a sol.
    return(line - 1)
  }
}

# Distribució Poisson truncada: http://freerangestats.info/blog/2018/03/20/truncated-poisson

# transform Qualimap output to R-object.
# diferencia respecte samtools: inclou zero coverage. Realment necessari?


  # Objectes de prova curts:
    obj <- c(1,1, 2,2, 3,3)
    obj <- c(rep(1,10), rep(2,10), rep(3,10))


  # LAB DATASET PATHS 
    # Qualistats inclou zero:
    obj <- read.table("Escritorio/genome_size/raw_data_qualimapReport/coverage_histogram.txt", header = TRUE)
    # exclou zero:
    obj <- obj[-1, ]
    
    # Carrega una femella de la Paula.
    obj <- read.table("Escritorio/Results/Gsize_Dscabricula/paula181_gstats.txt", skip=4, col.names = c("cov", "freq"))
    

  # HOME PC DATA PATHS
    obj <- read.table("Escriptori/u2_gstats.txt", skip = 5, col.names = c("cov", "freq"))
    obj <- read.table("Escriptori/u2_qualistats.txt", col.names = c("cov", "freq"))
    
    obj <- read.table("Escriptori/paula194_gstats.txt", skip = 4, col.names = c("cov", "freq"))
    obj <- read.table("Escriptori/paula194_qualistats.txt", col.names = c("cov", "freq"))
    
    obj <- read.table("Escriptori/deng_gstats.txt", skip = 4, col.names = c("cov", "freq"))
    
    obj <- read.table("Escriptori/Dscabricula_gstats.txt", skip = 5, col.names = c("cov", "freq"))
    
    # VADIM:
    obj <- read.table("Baixades/Deng_19mer_rawreads.histo", col.names = c("cov", "freq"))
    obj <- read.table("Baixades/Dsca_19mer_out.histo", col.names = c("cov", "freq"))
    obj <- read.table("Baixades/DsilPEsg3680196_raw_17_out.histo", col.names = c("cov", "freq"))

    
    plot(obj, xlim = c(1, 20), type = "l", lty=1, col = "green")
    legend(x = "topright", legend = c("Denghofii", "Dsilvatica", "Dscabricula"), fill = c("red", "green", "blue"))
    lines(obj, lty=1, col = "red")
    lines(obj, lty=3, col = "green")
    
  # Crear una taula més petita (mides genòmiques enormes pengen el programa):
  # Disminueix la mida del dataset però manté els percentatges de cada variable:
    quantitat_objectes_nova_taula = 10000  # Transforma de N entrades a q_o_n_t entrades.
    # factor de re-escalada.
    factor = sum(obj$freq)/quantitat_objectes_nova_taula
    # Transforma freqüència x al 10% al 10% de q_o_n_t.
    obj$freq = obj$freq/factor

  # Com que el coverage és una variable discreta, s'ha de manejar de la següent forma:
    obj <- expandRows(obj, "freq")  # Questa funció no treballa bé amb datasets grans.
    obj <- as.vector(obj$cov) 
    summary(obj)

# funció mediana de coverage del vector de freqüències.
median_zeroidx_cov(obj$freq)

# define function for mode:
# el valor de la 1era columna que te el màxim a la 2na columna.
mode <- function(x) {x[,1][x[,2]==max(x[,2])]}

# funció de la moda per a un objecte vectorial, amb valors discrets.
mode <- function(obj) {uniqv <- unique(obj) ; uniqv[which.max(tabulate(match(obj, uniqv)))]}

# Agafa un valor 5 per sota de mòdul. Si és més petit que zero, agafa el zero. 
mini <- function(x) ifelse(test = mode(x)-5>=0, yes = mode(x)-5, no = 1)
# Agafa un valor per sobre de 5 del mòdul.
maxi <- function(x) mode(x) + 5

dtruncated_poisson <- function(x, lambda) {dtrunc(x, "pois", a=min, b=max, lambda=lambda)} 
ptruncated_poisson <- function(q, lambda) {ptrunc(q, "pois", a=min, b=max, lambda=lambda)}

dtruncated_poisson(obj$freq, mode(obj))

unique(dtrunc(obj, "pois", a = mini(obj), b=maxi(obj), lambda=mode(obj)))

fitdist(obj, "truncated_poisson", start = list(lambda = mode(obj)))


 # Funcionen les aproximacions del vector a vectors més petits per calcular la Poisson?

complet = as.vector( c(rep(1,30), rep(2,20), rep(3,25), rep(4,30), rep(5,50), rep(6,60), rep(7,50), rep(8,45), rep(9,40), rep(10,30), rep(11,20), rep(12,10) ) )
partit = as.vector( c(rep(1,7.5), rep(2,5), rep(3,6.25), rep(4,7.5), rep(5,12.5), rep(6,15), rep(7,12.5), rep(8,11.25), rep(9,10), rep(10,7.5), rep(11,5), rep(12,2.5) ) )

 # Intenta fer un fit de distribució al model.
fitdistr(complet, "Poisson")
fitdistr(partit, "Poisson")

dtrunc(complet, "pois", a=mini(complet), b=maxi(complet))







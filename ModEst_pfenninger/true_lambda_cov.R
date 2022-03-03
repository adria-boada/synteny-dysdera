
library(fitdistrplus) 
library(truncdist) 
library(splitstackshape) 

# Distribució Poisson truncada: http://freerangestats.info/blog/2018/03/20/truncated-poisson

# transform Qualimap output to R-object.
# diferencia respecte samtools: inclou zero coverage.
obj <- read.table("Escritorio/Dscabricula_gstats.txt", skip = 5, col.names = c("cov", "freq"))
obj <- read.table("Escriptori/genome_size_dsilchru2/raw_data_qualimapReport/coverage_histogram.txt", header = TRUE)

# path pel pc del lab:
# inclou zero:
obj <- read.table("Escritorio/genome_size/raw_data_qualimapReport/coverage_histogram.txt", header = TRUE)
# exclou zero:
obj <- obj[-1, ]

 # Carrega una femella de la Paula.
obj <- read.table("Escritorio/Results/Gsize_Dscabricula/paula181_gstats.txt", skip=4, col.names = c("cov", "freq"))

 # objectes de prova:
obj <- c(1,1, 2,2, 3,3)
obj <- c(rep(1,10), rep(2,10), rep(3,10))

 # Crear una taula més petita:
factor = sum(obj$freq)/10000
obj$freq = obj$freq/factor


obj <- expandRows(obj, "freq") 
obj <- as.vector(obj$cov) 
summary(obj)

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

# mediana de coverage del vector de freqüències.
median_zeroidx_cov(obj$freq)

# define function for mode:
# el valor de la 1era columna que te el màxim a la 2na columna.
mode <- function(x) {x[,1][x[,2]==max(x[,2])]}

# funció de la moda per a un objecte vectorial amb els valors quantitzats.
mode <- function(obj) {uniqv <- unique(obj) ; uniqv[which.max(tabulate(match(obj, uniqv)))]}

# Agafa un valor 5 per sota de mòdul. Si és més petit que zero, agafa'l a ell.
mini <- function(x) ifelse(test = mode(x)-5>=0, yes = mode(x)-5, no = 0)
maxi <- function(x) mode(x) + 5

dtruncated_poisson <- function(x, lambda) {dtrunc(x, "pois", a=min, b=max, lambda=lambda)} 
ptruncated_poisson <- function(q, lambda) {ptrunc(q, "pois", a=min, b=max, lambda=lambda)}

dtruncated_poisson(obj$freq, mode(obj))

fitdist(obj$freq, "pois", start = list(lambda = mode(obj)))


 # Funcionen les aproximacions del vector a vectors més petits per calcular la Poisson?

complet = as.vector( c(rep(1,30), rep(2,20), rep(3,25), rep(4,30), rep(5,50), rep(6,60), rep(7,50), rep(8,45), rep(9,40), rep(10,30), rep(11,20), rep(12,10) ) )
partit = as.vector( c(rep(1,7.5), rep(2,5), rep(3,6.25), rep(4,7.5), rep(5,12.5), rep(6,15), rep(7,12.5), rep(8,11.25), rep(9,10), rep(10,7.5), rep(11,5), rep(12,2.5) ) )

 # Intenta fer un fit de distribució al model.
fitdistr(complet, "Poisson")
fitdistr(partit, "Poisson")

dtrunc(complet, "pois", a=mini(complet), b=maxi(complet))







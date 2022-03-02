
library(fitdistrplus) 
library(truncdist) 
library(splitstackshape) 

# transform Qualimap output to R-object.
# diferencia respecte samtools: inclou zero coverage.
obj <- read.table("Escriptori/u2_gstats.txt", skip = 5)
obj <- read.table("Escriptori/genome_size_dsilchru2/raw_data_qualimapReport/coverage_histogram.txt", header = TRUE)

# path pel pc del lab:
# inclou zero:
obj <- read.table("Escritorio/genome_size/raw_data_qualimapReport/coverage_histogram.txt", header = TRUE)
# exclou zero:
obj <- obj[-1, ]


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

# Agafa un valor 5 per sota de mòdul. Si és més petit que zero, agafa'l a ell.
min <- ifelse(test = mode(obj)-5>=0, yes = mode(obj)-5, no = 0)
max <- mode(obj) + 5

dtruncated_poisson <- function(x, lambda) {dtrunc(x, "pois", a=min, b=max, lambda=lambda)} 
ptruncated_poisson <- function(q, lambda) {ptrunc(q, "pois", a=min, b=max, lambda=lambda)}

dtruncated_poisson(obj$freq, mode(obj))

fitdist(obj$freq, "pois", start = list(lambda = mode(obj)))



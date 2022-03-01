
library(fitdistrplus) 
library(truncdist) 
library(splitstackshape) 

# transform Qualimap output to R-object.
# diferencia respecte samtools: inclou zero coverage.
obj <- read.table("Escriptori/u2_gstats.txt", skip = 4)
obj <- read.table("Escriptori/genome_size_dsilchru2/raw_data_qualimapReport/coverage_histogram.txt", header = TRUE)

obj <- expandRows(obj, "freq") 
obj <- as.vector(obj$cov) 
summary(obj) 
#define function for mode 
mode <- function(obj) {uniqv <- unique(obj) uniqv[which.max(tabulate(match(obj, uniqv)))]} 
#mode <- function(x)
#mode===>  
  {obj[,1][obj[,2]==max(obj[,2])]}


min <- mode - 5 
max <- mode + 5 
dtruncated_poisson <- function(x, lambda) {dtrunc(x, "pois", a=min, b=max, lambda=lambda)} 
ptruncated_poisson <- function(q, lambda) {ptrunc(q, "pois", a=min, b=max, lambda=lambda)} 
fitdist(obj, "pois", start = list(lambda = mode))
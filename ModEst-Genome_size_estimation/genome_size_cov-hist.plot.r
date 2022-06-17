
# Captura arguments del command-line:
args = commandArgs(trailingOnly=TRUE)

# Nom de la gràfica: extreure filename dels args.
filenme = basename(args[1])

# Llegeix la taula a partir de la segona fila
# (a la primera fila hi ha el "genome size")
x=read.table(args[1], skip = 5)

# Imprimeix màxims locals:
max_local = which(diff(sign(diff(x[,2])))==-2)+1
# Posa-hi a la grafica els primers 16 maxs ([1]+[2:16])
title_max = paste0(max_local[1])
for (i in max_local[2:16]) title_max = paste0(title_max, ", ", i)

print(paste0("These are the local maxima: ", title_max))

# Grafica:
pdf(paste0(filenme, "_coverage_hist.pdf"))
plot(x[,1], x[,2], log="x", type="l",
	xlab="Coverage", ylab="Count", 
	main = c(paste0("Max_local: ", title_max), paste0("file: ", filenme) ) 
)

dev.off()

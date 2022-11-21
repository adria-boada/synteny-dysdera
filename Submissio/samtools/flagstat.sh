
#! /bin/bash

# Exporta els camins fins a la ultima versio del soft.
# Es posa al davant del $PATH perque prengui prioritat.
export PATH="/users-d3/adria.boada/bin/samtools-1.14:$PATH"
# Comprova-ho...
which samtools

# Counts the number of alignments for each bit-wise FLAG type.
samtools flagstat $1


#! /bin/bash

# INTENCIÓ:
# --------
#
# Filtrar les lectures que no siguin primàries de l'aliniament (eliminar flag
# 256, que son els 'not primary alignment').
# Per fer-ho amb samtools: opció -F (filtrar fora) 256 (lectures no primàries).
#
# --------

# Inicialitzem algunes opcions...

# PATH
source /users-d3/adria.boada/.bashrc
# Comprova que uses la versió adequada del soft. que necessites...
echo "Samtools (tested with version 1.14):"
echo "  $(which samtools)"
echo

# INPUT PARSING
while [ "$1" != "" ]; do
  case $1 in
      -t )    shift
              threads=$1
              shift
              ;;
      -o )    shift
              out_primary_ali=$1
              shift
              ;;
      -i )    shift
              inpbam=$1
              shift
  esac
done

# DEFAULT THREADS
# (if the var is not set as an arg.)
if [ -z "$threads" ] ; then
    threads=6
fi
echo "Threads used: $threads"

# -F: exclude reads with any flag of 256 or 2048 (their sum)
samtools view -@ $threads -bF 2304 $inpbam > $out_primary_ali


#! /bin/bash

# INTENCIÓ:
# --------
#
# Filtrar les lectures que no siguin primàries de l'aliniament (eliminar flag
# 256, que son els 'not primary alignment').
# Per fer-ho amb samtools: opció -F (filtrar fora) 256 (lectures no primàries).
#
# Per a còrrer l'script: 
# `filter_alignment.sh file1.bam file2.bam ... fileN.bam [-t $threads]`
#
# Output names: 'primary_file1.bam'
# Run with -cwd (current working directory) being the folder storing the
# bamfiles, otherwise the path might get jammed by the name change.
# (i.e. primary_folder/bamfile.bam -> Error)
#
# --------

# Inicialitzem algunes opcions...

# PATH
source /users-d3/adria.boada/.bashrc
# Comprova que uses la versió adequada del soft. que necessites...
echo "Samtools (tested with version 1.14):"
echo "  $(which samtools)"
echo
echo "Current working directory:"
echo "  $(pwd)"

# INPUT PARSING
while [ "$1" != "" ]; do
  case $1 in
      -t )    shift
              threads=$1
              shift
              ;;
      * )     ARGS+=" $1"
              shift
  esac
done

# DEFAULT THREADS
# (if the var is not set as an arg.)
if [ -z "$threads" ] ; then
    threads=6
fi
echo "Threads used: $threads"

for fn in $(echo "$ARGS") ; do
  echo -e "Input: $fn\nOutput: primary_$fn"
  # -F: exclude reads with any flag of 256 or 2048 (their sum)
  samtools view -@ $threads -bF 2304 $fn > primary_${fn}
done


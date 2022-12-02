#! /bin/bash

# INTENCIÓ:
# --------
#
# A partir de múltiples bamfiles imprimeix un report (en format markdown)
# del coverage mitjà, del nombre de reads al mapatge, de la seva llargada, etc.
# Al final fés una gràfica histograma a partir del 'dataframe' de freqüències
# de posicions(?) amb tant de cobriment,
#
# La suma de freq*covg hauria de ser igual al total de bases mapades?
#
# --------

# Inicialitzem algunes opcions...
# PATH
source /users-d3/adria.boada/.bashrc
# Comprova que uses la versió adequada del soft. que necessites...
echo "samtools path: $(which samtools)" # samtools-1.14
echo "bwa path: $(which bwa)"						# /soft/bwa/bwa
echo "$(python3 --version)"							# python3.9.7 (at least >=3.6?)

# INPUT PARSING
while [ "$1" != "" ]; do
  case $1 in
      -t )    shift
              threads=$1
              ;;
      -r )    shift
              report_fname=$1
              ;;
      * )     ARGS+=" $1"
  esac
  shift
done

# DO NOT RUN WITHOUT A GIVEN OUTPUT FILENAME
if [ -z "$report_fname" ] ; then
    echo "No output filename was given with the '-r' option."
    echo "It should be a markdown file where results will be written."
    echo
    echo "Example:"
    echo "${0##*/} -r 'out-report-filename.md' -t 'threads' BAM-1 BAM-2 ... BAM-N"
    exit 1
fi

# DEFAULT THREADS IF LEFT UNSET IN ARGS
if [ -z "$threads" ] ; then
    threads=11
fi
echo "threads used: $threads"

# PER A CADA BAM-FILE:
# Extreu les banderes del mapatge
# (ha anat bé, el mapatge? percentatges mapats?)
samtools flagstat ${fn}.bam > ${fn}.stats.tmp

# Continua amb l'anàlisi de l'histograma de coverage (samtools stats):
# Extreu els stats de samtools i guarda'ls en un arxiu.
samtools stats --threads "$threads" ${fn}.bam >> ${fn}.stats.tmp


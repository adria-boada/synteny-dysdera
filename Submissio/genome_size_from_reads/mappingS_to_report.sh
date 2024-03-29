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
echo "Samtools (tested with version 1.14):"
echo "  $(which samtools)"
echo "Python3 (tested with version 3.9.7)"
echo "  $(python3 --version)"
echo "Make sure that the Shebang of the python3 file is congruent with its path..."
echo "  $(which python3)"
# Per a còrrer un script.py, activa'l amb chmod i busca'l al path amb
# `which script.py` per assegurar que usa la direcció adequada. Evita
# a tot cost direccions absolutes.
echo "Homemade script nº 1 (calculates means):"
echo "  $(which covfreq_tomeans.py)"
echo "Homemade script nº 2 (parses alignment info.):"
echo "  $(which samstats_ofinterest.sh)"
echo

# INPUT PARSING
while [ "$1" != "" ]; do
  case $1 in
      -t )    shift
              threads=$1
              shift
              ;;
      -r )    shift
              report_fname=$1
              shift
              ;;
      # La palanca '-h' activa la impressio
      # d'un histograma inicial [INHABILITAT]
      -h )    shift
              histo='TRUE'
              ;;
      * )     ARGS+=" $1"
              shift
  esac
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

# DEFAULT THREADS
# (if the var is not set as an arg.)
if [ -z "$threads" ] ; then
    threads=11
fi
echo "Threads used: $threads"
echo "The file '$report_fname' will be overwritten"
echo
echo '' > $report_fname

# PER A CADA BAM-FILE:
# echo $ARGS fa clars els espais pel bon funcionament del for-loop.
for fn in $(echo "$ARGS") ; do
    # Elimina la extensió final (.bam) del fitxer:
    fn=${fn%\.*}
    echo "Etapa: $fn"
    # Extreu les 'banderes' del mapatge
    # (ha anat bé, el mapatge? -> percentatge de lectures mapades)
    samtools flagstat --threads "$threads" ${fn}.bam > ${fn}.stats.tmp
    # Continua amb l'anàlisi de l'histograma de coverage (samtools stats):
    # Extreu els stats de samtools i guarda'ls en un arxiu.
    samtools stats --threads "$threads" ${fn}.bam >> ${fn}.stats.tmp

    # Append els resultats de .stats.tmp al report, 
    # en ordre i en format markdown.
    echo -e "\n## INPUT-FILE: ${fn}\n" >> $report_fname
    echo -e "### General mapping stats\n" >> $report_fname
    samstats_ofinterest.sh ${fn}.stats.tmp >> $report_fname

    # Retalla dades d'interés (read bp cuml. sum).
    readcumsum=$(tail $report_fname | grep "Raw reads' cum" | cut -d' ' -f6)

    # Modifica el format de l'histograma cru per adequar-lo als
    # requeriments del guió python3 que en calcula la freqüència mitjana.
    # El fitxer amb l'histograma adequat és ${fn}_covg_histogram.csv
    echo "Coverage,Frequency" > ${fn}_covg_histogram.csv
    # Recull la secció 'COVERAGE' dels stats i retalla'n les columnes:
    grep "^COV" ${fn}.stats.tmp | cut -f3- |
    	  head -n-1 | # Elimina la última fila (reads >1000 profunditat)
        # No els podem incorporar a l'anàlisi pq no es pot multiplicar
        # profunditat de '>1000' per cap freqüència (no és un número,
        # és un rang que inclou totes les freqs superior a 1000).
        tr '\t' ',' >> ${fn}_covg_histogram.csv # subst. TABS per comes
    # Empra el guió de python3 per fer l'anàlisis estadístic.
    echo -e "\n### Coverage frequencies\n" >> $report_fname
    covfreq_tomeans.py ${fn}_covg_histogram.csv >> $report_fname

    # Retalla dades d'interés (observed coverage)
    obscov=$(tail $report_fname | grep "Mean: " | cut -d' ' -f3)

    # Resumeix les dades més importants, per a una ràpida obtenció:
    # '\numprint{}' és una funció de LaTeX per presentar nombres grans.
    gensize=$(python3 -c "print($readcumsum / $obscov)")
    echo ''; echo "| $(basename $fn) | \\numa{$readcumsum} | x\\numa{$obscov} | \\numa{$gensize} |" >> $report_fname
    # Necessitem pandas per còrrer anàlisis...!!!
    # Aquesta mena de fer taules és el pitjor tuestón.

    # Un cop els anàlisis acaben, elimina els fitxers temporals:
    ##rm --force ${fn}_covg_histogram.csv
    ##rm --force ${fn}.stats.tmp
    # Recorda que més tard vols fer el plotting amb pandas
    # dels histogrames un sobre l'altre.

done

# HISTOGRAMA PER A CADA SET DE DADES BAM
# Un sol plot amb tots els histogrames recollits.
# Possiblement millor si el plot es troba a l'inici (pag.1) del report...
# Associar una palanca per incloure o excloure la figura...


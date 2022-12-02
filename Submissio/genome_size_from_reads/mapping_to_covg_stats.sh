
#! /bin/bash

### Inicialitzar algunes opcions...
# Threads
threads=11
echo "threads used: $threads"

	# Exporta els camins fins a la ultima versio del soft.
	# Es posa al davant del $PATH perque prengui prioritat.
	#export PATH="/soft/bwa:$PATH"
	#export PATH="/users-d3/adria.boada/bin/samtools-1.14:$PATH"
# El paràgraf tabulat de dalt l'estalvies carregant el .bashrc
# (on ja apuntes els camins de tots els programes).
source /users-d3/adria.boada/.bashrc
# Comprova que uses la versió adequada del soft. que necessites...
echo "samtools path: $(which samtools)" #samtools-1.14
echo "bwa path: $(which bwa)"		#/soft/bwa/bwa
echo "$(python3 --version)"		#python3.9.7 (at least >=3.6?)

# Ref assembly.
ref_asm=$2

# Read-file(s) queried to map (one file for single-end or two files for paired-ends)
reads=$3
mates=$4

# Prepare filenames for the output.
# The filename gets the first four characters from $ref_asm and $lec1.
#~name_ref=${ref_asm##*/}
#~name_read=${reads##*/}
#~fn="RF_${name_ref:0:4}_LC_${name_read:0:4}"
fn=$1  # easier

# Launch bwa-mem and transmutate its SAM output to BAM thanks to samtools.
# Make sure you have updated samtools' software and you also have it incorporated in your $PATH variable.
# Test it with `which samtools` and `echo $PATH | grep samtools`

### Opcions bwa-mem:
# -a    : Output all found alignments for single-end or unpaired paired-end reads.
#         These alignments will be flagged as secondary alignments.
# -c INT: Discard a `MEM` if it has more than `INT` occurrences in the genome.
#         It is an insensitive parameter. Default: 10000.
# -t INT: Quantes CPUs tens pensades emprar?
# -o FILE: Nom del fitxer on redireccionar l'output (format SAM).

# Opcions samtools: 
# No es pot fer directament amb `sort`, has de passar per `view`?
# `sort` no accepta [in.sam], tant sols [in.bam].
# -F    : Descarta els aliniaments/mapatges marcats amb la 'bandera' 4 (flag de segment unmapped).
#         Una f minúscula voldria dir 'quedar-se', majúscula 'descarta'.
#         La 'bandera' 4: lectures que no han mapat.
# -b    : Transforma a [out.bam], igual que '--bam'.
# -u    : Transforma a [out.bam] uncompressed (recomanat per a pipelines)
# -@    : Quantes CPUs tens pensades emprar?
# -o    : Sense especificar res, la opció predeterminada és `stdout`.
# -T    : Fitxer FASTA de referència.
# El guionet '-' al final de la ordre de samtools li diu que llegeixi stdin.

### Comença el mapatge amb bwa... 
bwa mem -a -t $threads $ref_asm "$reads" "$mates" |
	# Possibilitat d'evaluar la canonada de bwa gracies a tee:
	#tee tmp.sam |
	samtools view -T $ref_asm -@ "$threads" -u - |
	samtools sort -@ "$threads" -o ${fn}.bam -

# Extreu les banderes del mapatge
# (ha anat bé, el mapatge? percentatges mapats?)
samtools flagstat ${fn}.bam > ${fn}.stats.tmp

### Continua amb l'anàlisi de l'histograma de coverage (samtools stats):
# Extreu els stats de samtools i guarda'ls en un arxiu.
samtools stats --threads "$threads" ${fn}.bam >> ${fn}.stats.tmp

# Obtén i imprimeix les dades d'interés:
# Assegurar-se de que els fitxers $reads i $mates
# son fitxers comprimits (*gz, comanda zcat).
# Per aconseguir nombre de reads,
# es divideix el nombre de línies del
# fitxer fq.gz per quatre
# (cada read és format per quatre línies)

echo "## General mapping stats"
echo
# Treure'ls dels stats generals:
qr=$(grep "^SN" ${fn}.stats.tmp|cut -f2-|grep "raw total sequences")
echo "+ Number of raw reads queried: $(echo $qr|cut -f4 -d' ')"
prim=$(grep "primary mapped" ${fn}.stats.tmp|cut -f1 -d' '|tr -d '(')
echo "+ Primary alignments: $prim"
nprim=$(grep "^SN" ${fn}.stats.tmp|cut -f2-|grep "non-primary alignments")
echo "+ Non-primary alignments: $(echo $nprim|cut -f3 -d' ')"
supp=$(grep "^SN" ${fn}.stats.tmp|cut -f2-|grep "supplementary alignments")
echo "+ Supplementary alignments: $(echo $supp|cut -f3 -d' ')"
nuc=$(grep "^SN" ${fn}.stats.tmp|cut -f2-|grep "total length:")
echo "+ Raw read length: $(echo $nuc|cut -f3 -d' ')"
nuc=$(grep "^SN" ${fn}.stats.tmp|cut -f2-|grep "(cigar):")
echo "+ Bases mapped (more accurate): $(echo $nuc|cut -f4 -d' ')"
err=$(grep "^SN" ${fn}.stats.tmp|cut -f2-|grep "error rate")
echo "+ Error rate: $(echo $err|cut -f3 -d' ')"
insert=$(grep "^SN" ${fn}.stats.tmp|cut -f2-|grep "insert size")
echo "+ Insert size and its standard deviation: $(echo $insert|cut -f4 -d' ') ~~$(echo $insert|cut -f9 -d' ')"
perc=$(grep "primary mapped" ${fn}.stats.tmp|cut -f6 -d' '|tr -d '(')
echo "* Percentage of primary reads mapped: $perc"
echo
echo "-- Table with a priori and a posteriori expected coverages for 3Gb and 1.5Gb--"
echo

echo "## Coverage histogram stats"
echo
echo "Coverage,Frequency" >${fn}_covg_histogram.csv
grep "^COV" ${fn}.stats.tmp | cut -f3- |
	# Modifica l'histograma cru per formatar-lo segons les necessitats
	# del guió de python3...
	head -n-1 | # Elimina la última fila (reads >1000 profunditat)
	# No els podem incorporar a l'anàlisi pq no es pot multiplicar
	# profunditat de '>1000' per cap freqüència (no és un número,
	# és un rang que inclou totes les freqs superior a 1000).
	tr '\t' ',' >>${fn}_covg_histogram.csv # subst. TABS per comes
# Empra el guió de python3 per fer l'anàlisis estadístic.
# (imprimeix als mateixos logs, directament).
python3 /users-d3/adria.boada/home/Submissio/covfreq_tomeans.py ${fn}_covg_histogram.csv
echo

# Once the analyses are finished, remove tmp files...
#~rm ${fn}_covg_histogram.csv


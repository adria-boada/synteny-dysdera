
#! /bin/bash/

### cami fins samtools versio 1.14
export PATH="/users-d3/adria.boada/bin/samtools-1.14:$PATH"
# comprova-ho...
which samtools

threads=6

### Gestio ARGS/ arguments del guio.
if [ -z "$1" ] ; then
	echo "Input a sorted BAM to analyze as the first argument."
	exit 1
else
	input_bam="$1"
fi

### Per a l`output, anomena`l tallant l`extensio i el cami de l`input. 
name_bam=${input_bam%.bam} ; name_bam=${name_bam##*/}
# No ens fa falta si es guarda stdout a logs.
# (opció -o de la comanda qsub).

### Extreu els stats de samtools i guarda'ls en un arxiu.
samtools stats --threads "$threads" $input_bam > ${name_bam}.stats.tmp

### Imprimeix el nom de l'arxiu d'entrada:
# (no es perdi sobre qui fas l'anàlisi)
echo " # INPUT: ${input_bam##*/}"
echo

### Extreu el nombre de `mapped nucleotides` i el nombre de `total sequenced nucleotides`
seq_nuc=$(cat ${input_bam}.stats.tmp | grep "^SN" | cut -f2- | grep "total length:")
#elimina espais entremig, separa el número de dins la línia:
seq_nuc=$(echo $seq_nuc|cut -d' ' -f3)
map_nuc=$(cat ${input_bam}.stats.tmp | grep "^SN" | cut -f2- | grep "(cigar):")
#elimina espais entremig, separa el número de dins la línia:
map_nuc=$(echo $map_nuc|cut -d' ' -f4)

### Retorna la mida del mapatge i sequenciació:
echo "SUM OF SEQUENCED READS NUCLEOTIDES: ${seq_nuc%#*}"
echo "SUM OF MAPPED BASES AGAINST THE REF.GEN.: ${map_nuc%#*}"
echo "Per aconseguir la mida dels fitxers *.fastq:"
echo "echo \$(zcat *.fq.gz | wc -l) / 4 |bc"
echo "Compta el nombre de línies del fitxer i divideix per quatre (cada quatre linies forma un read).

### Fòrmula per aproximar genome size.
# Sembla acabar sovint per sobre de la inferència per kmers.
echo -e "Formula that estimates genome size:"
echo "genome size = (bases mapped) / (modal coverage)"

### Imprimeix dades bàsiques. 
cat ${input_bam}.stats.tmp | grep "^SN" | cut -f2- 
echo # espai per separar info general de l'histo-coverage.

### Extreu els punts d'un histograma del coverage:
echo ; echo "Coverage histogram"
cat ${input_bam}.stats.tmp | grep "^COV" | cut -f2-


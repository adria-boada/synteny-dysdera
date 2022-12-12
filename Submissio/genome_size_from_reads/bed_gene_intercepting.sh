#! /bin/bash

# INTENCIÓ:
# ---------
#
# A partir de parelles de coordenades BUSCO i aliniaments BAM, busca
# els reads dins el BAM que intersecten amb les coordenades BUSCO dins el BED.
#
# ---------

# Inicialitzem opcions:

# PATH
source /users-d3/adria.boada/.bashrc
echo "Bedtools (tested with v2.30.0)"
echo "$(which bedtools) $(bedtools --version|cut -f2 -d' ')"
echo

# INPUT PARSING
# Mentre existeixi argument nº 1:
while [ "$1" != "" ]; do
	# Revisa si és un bam o un bed file:
	# `shift` mou els arguments una posició endavant.
	# ($1 'desapareix', i $2 passa a ser $1)
	case $1 in
		*bam )	bamfile=$1
						shift
						;;
		*bed )	bedfile=$1
						shift
						;;
		# Si no és cap dels dos, surt del guió:
		* )
			echo "A file with neither bam nor bed extension was detected"
			echo "Make sure the arguments are correct"
			exit 1
	esac
	# Si aconsegueixes una parella de bed i bam, (bed && bam)
	# executa la intersecció entre ells dos.
	if [[ "$bedfile" && "$bamfile" ]]; then
		# Aquest nom no és capaç de diferenciar la intersecció
		# del mateix BAM amb diversos BEDs (un pel primer
		# exó, un pel gen sencer, el mateix BAM, output sobrescrit)
		# solució: córrer en dos etapes i reanomenar manualment.
		newfile_name="BUSCO_intersct_${bamfile##*/}"
		echo "Creating $newfile_name"
		echo "from BED $bedfile"
		echo "and BAM $bamfile"
		echo
		bedtools intersect -a $bamfile \
			-b $bedfile > $newfile_name
		bedfile=""
	fi
done


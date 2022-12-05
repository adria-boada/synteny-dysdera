#! /bin/bash

# OBJECTIVE
# ---------
#
# Create bwa mem mappings from many pairs of FQ--FASTA files,
# each with a newly given name (by a third variable with BAM extension).
#
# Example:
# ./bwa_mem.sh reads.fq mates.fq assembly.fasta outname.bam \
#    outname.bam mates.fq assembly.fasta reads.fq
# The above command would create two mappings with filenames `outname.bam`
#
# ---------

# Initialize options...

# PATH
source /users-d3/adria.boada/.bashrc
# Comprova que uses la versió adequada del soft. que necessites...
echo "Samtools (tested with version 1.14):"
echo "  PATH: $(which samtools)"
echo "  ver.: $(samtools --version)"
echo "Burrows-Wheeler Aligner (bwa)"
echo "  PATH: $(which bwa)"
echo

# INPUT PARSING
while [ "$1" != "" ]; do
	case $1 in
		-t )
			shift
			threads=$1
			;;
		# Enumerar possibles extensions de reads and mates
		*fq | *fq.gz | *fastq.gz | *fastq )
			if [ -z "$reads" ]; then
				reads=$1
			else
				mates=$1
			fi
			;;
		# Enumerar possibles extensions d'assemblies
		*fasta | *fas | *fa | *fa.gz | *fasta.gz )
			ref_asm=$1
			;;
		# Anomena al bam output
		*bam )
			output=$1
			;;
		* )
			echo "A file with an unrecognised extension was detected"
			echo "Make sure the arguments or the script are correct"
			exit 1
	esac
	# shift moves to the next argv. ($1 <- $2 <- $3 etc.)
	shift

	if [ -z "$threads" ]; then
		echo "Amount of threads has not been set as the first variable."
		echo "Remember to set threads with the '-t INT' option."
		exit 1
	fi

	# MAIN ALIGNMENT SCRIPT
	# Si totes les vars necessàries per fer alignment han sigut definides:
	if [[ "$reads" && "$mates" && "$ref_asm" && "$output" ]]; then
		# Endavant.
		echo "Starting alignment..."
		echo "Reads: $reads"
		echo "Mates: $mates"
		echo "Reference assembly: $ref_asm"
		echo "Output name: $output"
		echo

	bwa mem -a -t $threads $ref_asm $reads $mates |
		# Possibilitat d'evaluar la canonada de bwa gracies a tee:
		# canalitza i al mateix temps escriu una còpia del que surt de bwa mem
		# al fitxer tmp.sam.
		#tee tmp.sam |
		# Compress from SAM to BAM format (-u)
		samtools view -T $ref_asm -@ $threads -u - |
		# Sort the BAM file.
    samtools sort -@ $threads -o $output -
	# Delete reads variable...
	reads=""
	mates=""
	output=""
	fi
done

# Launch bwa-mem and transmutate its SAM output to BAM thanks to samtools.
# Make sure you have updated samtools' software and you also have it incorporated in your $PATH variable.
# Test it with `which samtools` and `echo $PATH | grep samtools`

# Opcions bwa-mem:
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


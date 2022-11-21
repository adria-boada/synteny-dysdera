
# !/bin/bash

# Threads
threads=11

# Exporta els camins fins a la ultima versio del soft.
# Es posa al davant del $PATH perque prengui prioritat.
export PATH="/soft/bwa:$PATH"
export PATH="/users-d3/adria.boada/bin/samtools-1.14:$PATH"
which samtools
which bwa

# Ref assembly.
ref_asm=$1

# Read-file(s) queried to map (one file for single-end or two files for paired-ends)
reads=$2
mates=$3

# Prepare filenames for the output.
# The filename gets the first four characters from $ref_asm and $lec1.
name_ref=${ref_asm##*/}
name_read=${reads##*/}
fn="RF_${name_ref:0:4}_LC_${name_read:0:4}"

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

# -F'4' method... envia a un fitxer temporal (tmp.sam).
bwa mem -a -t $threads $ref_asm "$reads" "$mates" |
	# Possibilitat d'evaluar la canonada de bwa gracies a tee:
	#tee tmp.sam |
	samtools view -F 4 -T $ref_asm -@ "$threads" -u - |
        samtools sort -@ "$threads" -o ${fn}.bam -


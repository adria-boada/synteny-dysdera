#! /bin/bash -x
#$ -cwd
#$ -V
#$ -q h12.q
#$ -l hostname="hercules12"
#$ -N "mmap_sam_asm"

# add threads in command line
# -pe ompi128h12 24
# minimap takes one IO thread: qsub N, minimap N-1

# command:
# qsub -pe ompi128h12 26 Submissio/minimap2/minimap_sam_asm.sh sil cat sil_v_cat 24

# samtools and minimap path to exec.
samT="Programes/samtools-1.14/samtools"
mmp="/users-d3/adria.boada/bin/minimap2-2.24_x64-linux/minimap2"

#3: Output name (*.sam/bam)

#1: Target sequence
#2: Query sequence

#4: Number of threads

echo "ORDRE: "
echo "~/Programes/minimap2-2.24_x64-linux/minimap2 -t $4 -ax asm20 $1 $2 > $3.sam"

$mmp -t $4 -ax asm20 $1 $2 | $samT sort -o $3.bam

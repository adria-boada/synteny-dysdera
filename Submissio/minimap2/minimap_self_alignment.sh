#! /bin/bash -x
#$ -cwd
#$ -V
#$ -q h12.q
#$ -l hostname="hercules12"
#$ -N "mmap_self"

# add threads in command line
# -pe ompi128h12 24
# minimap takes one IO thread: qsub N, minimap N-1

mmp="/users-d3/adria.boada/bin/minimap2-2.24_x64-linux/minimap2"

#3: Output name (*.paf)

#1: Target sequence
#2: Query sequence

#4: Number of threads

echo "ORDRE: "
echo "~/Programes/minimap2-2.24_x64-linux/minimap2 -t $3 -DP -k19 -w19 -m200 $1 $1 > $2.paf"

$mmp -t $3 -DP -k19 -w19 -m200 $1 $1 > $2.paf

# minimap2 -cx asm20 --cs $1 $2 > $3    # opcio completa

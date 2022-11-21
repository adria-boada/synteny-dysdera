#! /bin/bash -x
#$ -cwd
#$ -V
#$ -q h12.q
#$ -l hostname="hercules12"
#$ -N "mmap_pafasm"

# add threads in command line
# -pe ompi128h12 24
# minimap takes one IO thread: qsub N, minimap N-1

mmp="/users-d3/adria.boada/bin/minimap2-2.24_x64-linux/minimap2"

#3: Number of threads

#1: Target sequence
#2: Query sequence


echo "ORDRE: "
echo "$mmp -t $4 -cx asm20 $1 $2 > target_${outT1}.paf"

 # Retalla el cami i extensio dels fitxers.
ouT1=${1%%.*} ; outT1=${ouT1##*/}
ouT2=${2%%.*} ; outT2=${ouT2##*/}

echo "Paf-file names: target_${outT1}.paf and target_${outT2}.paf"

 # Target $1
$mmp -t $3 -cx asm20 $1 $2 > target_${outT1}.paf
 # Target $2
$mmp -t $3 -cx asm20 $2 $1 > target_${outT2}.paf


# minimap2 -cx asm20 --cs $1 $2 > $3    # opcio completa

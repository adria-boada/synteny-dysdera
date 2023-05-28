#! /bin/bash -x

source /users-d3/adria.boada/.bashrc

echo ...

which minimap2

# asm20
minimap2 -t $threads -cx asm20 --cs $fasta1 $fasta2 > ${title1}_asm20.paf
minimap2 -t $threads -cx asm20 --cs $fasta2 $fasta1 > ${title2}_asm20.paf

# asm10
minimap2 -t $threads -cx asm10 --cs $fasta1 $fasta2 > ${title1}_asm10.paf
minimap2 -t $threads -cx asm10 --cs $fasta2 $fasta1 > ${title2}_asm10.paf


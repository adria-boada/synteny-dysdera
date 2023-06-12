#! /bin/bash -x

# source PATH to minimap2
source /users-d3/adria.boada/.bashrc
echo "minimap2 version number: $(minimap2 --version)"
echo "minimap2 path: $(which minimap2)"

# define vars
fasta1=$1      # dtil fasta
fasta2=$2      # dcat fasta
threads=$3     # threads
title1="tDcat_qDtil"
title2="tDtil_qDcat"

# asm20
echo "Computing $fasta1 vs $fasta2 with settings asm20"
minimap2 -t $threads -cx asm20 --cs $fasta1 $fasta2 > ${title2}_asm20.paf
echo "Computing $fasta2 vs $fasta1 with settings asm20"
minimap2 -t $threads -cx asm20 --cs $fasta2 $fasta1 > ${title1}_asm20.paf

# asm10
echo "Computing $fasta1 vs $fasta2 with settings asm10"
minimap2 -t $threads -cx asm10 --cs $fasta1 $fasta2 > ${title2}_asm10.paf
echo "Computing $fasta2 vs $fasta1 with settings asm10"
minimap2 -t $threads -cx asm10 --cs $fasta2 $fasta1 > ${title1}_asm10.paf


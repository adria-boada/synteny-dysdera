#!/bin/bash
#
# qsub -cwd -V -N 'coverage_"region"' -j y -q h0809.q -l hostname=hercules09 -pe ompi24 8 -b y bash Submissio/samt_coverage.sh "region"

source .bash_profile

~/Programes/samtools-1.14/samtools coverage -r $1 -b Data/bamlist_Deng.txt

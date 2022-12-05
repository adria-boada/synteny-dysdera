#! /bin/bash

# Indexes the FASTA reference assemblies.
# It is required in order for bwa mem to work properly.
# Call the script as follows:
# bwa_index.sh ref1.fa ref2.fa ... refN.fa
# It will index each ref.fa passed as argv.

# PATH
echo "Burrows-Wheeler Aligner (bwa)"
echo "  PATH: $(which bwa)"
echo

# INDEXING
while [ "$1" != "" ]; do
	echo "Indexing $1"
	bwa index $1
	shift
done


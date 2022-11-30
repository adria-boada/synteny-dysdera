
#! /bin/bash

source /users-d3/adria.boada/.bashrc
echo "bamtools path: $(which bamtools)"

echo "arg-one: bam mapping/alignment that should be indexed"
echo "  --->  $1"

# index bamfile.
bamtools index -in $1

# randomly sample fragments from it.
bamtools random -in $1 -out "sampled100_$1" -n 100


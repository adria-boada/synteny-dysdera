
#! /bin/bash

source /users-d3/adria.boada/.bashrc
echo "bamtools path: $(which bamtools)"

echo "arg-one: bam mapping/alignment that should be indexed"
echo "  --->  $1"

# index bamfile.
bamtools index -in $1

# randomly sample fragments from it,
# creating a test-file.bam in the process.
echo "a random sample of 100 lines will be performed and the result"
echo "will be dumped to 'sampled100rando_$1'"
bamtools random -in $1 -out "sampled100rando_$1" -n 100


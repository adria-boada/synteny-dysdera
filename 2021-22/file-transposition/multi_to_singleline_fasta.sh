
#! /bin/bash

# Input the path to the fasta of interest.
fasta=$1

# Prints in stdout the single-lined-file.
  # The first sed command surrounds each header with '@'
  # The second one deletes all occurrences of '\n' char.
  # Next transliterate transforms @ to \n.
  # To cap it off, deletes the first line and adds a closing newline to the file.
    # sed '$a\' adds a newline at the end of the file if there is none.

sed '/^>/ {s/$/@/g; s/^/@/g}' $fasta | tr -d '\n' | tr '@' '\n' | sed '1d;$a\'

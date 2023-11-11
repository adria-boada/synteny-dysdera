#! /bin/bash

## PARSE ARGUMENTS ##

# Parse command-line arguments to define two FASTA files and which labels
# should be used for each FASTA.
while [ "$1" != "" ]; do
  case $1 in
    -t | --threads )
      shift
      THREADS=$1
      ;;
    --one )
      shift
      LABEL1="$(echo $1 | cut -d, -f2)"
      FILENAME1="$(echo $1 | cut -d, -f1)"
      coma1=$(echo $1 | grep , -o)
      ;;
    --two )
      shift
      LABEL2="$(echo $1 | cut -d, -f2)"
      FILENAME2="$(echo $1 | cut -d, -f1)"
      coma2=$(echo $1 | grep , -o)
      ;;
    * )
      echo "Please call the script as follows:"
      echo "${0##*/} --one FILEPATH,LABEL --two FILEPATH,LABEL [-t THREADS or --threads THREADS]"
      echo
      echo "The options '--one' and '--two' point towards two different genomes in FASTA format. These will be aligned with minimap2."
      echo "+ LABEL: A shorthand for the species name that will be used to label the output alignments."
      echo "+ FILEPATH: The path to the FASTA files."
      echo "+ THREADS [optional]: The amount of threads that will be allowed for minimap2 to use."
      echo "Consider that the computer must possess spare threads for I/O purposes;"
      echo "do not allow all threads to minimap2! By default, allows 24 (23 + 1 I/O)."
      exit 1
      ;;
  esac
  shift
done
# Make sure both --one and --two are well defined
if [ -z "$LABEL1" ] || [ -z "$LABEL2" ] || [ -z "$FILENAME1" ] || [ -z "$FILENAME2" ]; then
  echo "ERROR: One of the two genomes or labels was not correctly defined!"
  exit 1
elif [ "$coma1" != "," ] || [ "$coma2" != "," ] ; then
  echo "ERROR: One of the two input expressions (LABEL,FASTA) did not contain a comma!"
  exit 1
fi
# If the --threads lever has not been parameterised, set $THREADS to '23'.
if [ -z "$THREADS" ]; then
  THREADS=23
fi

## STARTUP PREPARATIONS ##

# source PATH to minimap2
source /users-d3/adria.boada/.bashrc
echo "minimap2 version number: $(minimap2 --version)"
echo "minimap2 path: $(which minimap2)"

echo -e "\nINFO: 6 mappings will be produced with varying -x asmY presets; asm20, asm10
and asm5. Both given genomes will be ran as query and target, which can produce
different results; mapping is not commutative. The option --cs (improve mapping
at the expense of computational resources) and -c (output PAF file instead of
SAM) have been enabled. It is possible to increase the mapping precision by
changing the default settings of 'asm20', 'asm10' or 'asm5' to better fit the data."

echo ==================================================
echo

## RUNNING MINIMAP2 ##

# Use the preset 'asm20' (up to 15% divergence)
echo "Computing $LABEL1 as TARGET against $LABEL2 as QUERY with settings asm20"
out="t${LABEL1}_q${LABEL2}_asm20.paf"
echo "Sending output to $out"
minimap2 -t $THREADS -cx asm20 --cs $FILENAME1 $FILENAME2 > $out
echo "Computing $LABEL2 as TARGET against $LABEL1 as QUERY with settings asm20"
out="t${LABEL2}_q${LABEL1}_asm20.paf"
echo "Sending output to $out"
minimap2 -t $THREADS -cx asm20 --cs $FILENAME2 $FILENAME1 > $out
echo ----------------------------------------

# Use the preset 'asm10' (up to ~5% divergence)
echo "Computing $LABEL1 as TARGET against $LABEL2 as QUERY with settings asm10"
out="t${LABEL1}_q${LABEL2}_asm10.paf"
echo "Sending output to $out"
minimap2 -t $THREADS -cx asm10 --cs $FILENAME1 $FILENAME2 > $out
echo "Computing $LABEL2 as TARGET against $LABEL1 as QUERY with settings asm10"
out="t${LABEL2}_q${LABEL1}_asm10.paf"
echo "Sending output to $out"
minimap2 -t $THREADS -cx asm10 --cs $FILENAME2 $FILENAME1 > $out
echo ----------------------------------------

# Use the preset 'asm5' (less than ~5% divergence)
echo "Computing $LABEL1 as TARGET against $LABEL2 as QUERY with settings asm5"
out="t${LABEL1}_q${LABEL2}_asm5.paf"
echo "Sending output to $out"
minimap2 -t $THREADS -cx asm5 --cs $FILENAME1 $FILENAME2 > $out
echo "Computing $LABEL2 as TARGET against $LABEL1 as QUERY with settings asm5"
out="t${LABEL2}_q${LABEL1}_asm5.paf"
echo "Sending output to $out"
minimap2 -t $THREADS -cx asm5 --cs $FILENAME2 $FILENAME1 > $out


#! /bin/bash

# PARSE ARGUMENTS
# Define FASTA files and their titles/labels.
while [ "$1" != "" ]; do
  case $1 in
    -t | --threads )
      shift
      THREADS=$1
      ;;
    --one )
      shift
      LABEL1="$(echo $1 | cut -d, -f1)"
      FILENAME1="$(echo $1 | cut -d, -f2)"
      coma1=$(echo $1 | grep , -o)
      ;;
    --two )
      shift
      LABEL2="$(echo $1 | cut -d, -f1)"
      FILENAME2="$(echo $1 | cut -d, -f2)"
      coma2=$(echo $1 | grep , -o)
      ;;
    * )
      echo "Please call the script as follows:"
      echo "${0##*/} --one LABEL,FASTA --two LABEL,FASTA [-t THREADS or --threads THREADS]"
      echo
      echo "The label will be used to create the output filenames."
      echo "'FASTA' refers to a file containing a genome sequence."
      echo "Default threads are 24 (23 plus one thread for I/O)."
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
  echo "ERROR: One of the two input expressions (LABEL,FASTA) did not contain a comma"
  exit 1
fi

# source PATH to minimap2
source /users-d3/adria.boada/.bashrc
echo "minimap2 version number: $(minimap2 --version)"
echo "minimap2 path: $(which minimap2)"

echo -e "\nINFO: 6 mappings will be produced with varying -x asmY presets; asm20, asm10
and asm5. Both given genomes will be ran as query and target, which can produce
different results; mapping is not commutative. The option --cs (improve mapping
at the expense of computational resources) and -c (output PAF file instead of
SAM) have been enabled."

echo ==================================================
echo

# asm20
echo "Computing Query.$LABEL1 vs Target.$LABEL2 with settings asm20"
out="q${LABEL1}_t${LABEL2}_asm20.paf"
echo "Sending output to $out"
minimap2 -t $THREADS -cx asm20 --cs $FILENAME1 $FILENAME2 > $out
echo "Computing Query.$LABEL2 vs Target.$LABEL1 with settings asm20"
out="q${LABEL2}_t${LABEL1}_asm20.paf"
echo "Sending output to $out"
minimap2 -t $THREADS -cx asm20 --cs $FILENAME2 $FILENAME1 > $out
echo ----------------------------------------

# asm10
echo "Computing Query.$LABEL1 vs Target.$LABEL2 with settings asm10"
out="q${LABEL1}_t${LABEL2}_asm10.paf"
echo "Sending output to $out"
minimap2 -t $THREADS -cx asm10 --cs $FILENAME1 $FILENAME2 > $out
echo "Computing Query.$LABEL2 vs Target.$LABEL1 with settings asm10"
out="q${LABEL2}_t${LABEL1}_asm10.paf"
echo "Sending output to $out"
minimap2 -t $THREADS -cx asm10 --cs $FILENAME2 $FILENAME1 > $out
echo ----------------------------------------

# asm5
echo "Computing Query.$LABEL1 vs Target.$LABEL2 with settings asm5"
out="q${LABEL1}_t${LABEL2}_asm5.paf"
echo "Sending output to $out"
minimap2 -t $THREADS -cx asm5 --cs $FILENAME1 $FILENAME2 > $out
echo "Computing Query.$LABEL2 vs Target.$LABEL1 with settings asm5"
out="q${LABEL2}_t${LABEL1}_asm5.paf"
echo "Sending output to $out"
minimap2 -t $THREADS -cx asm5 --cs $FILENAME2 $FILENAME1 > $out


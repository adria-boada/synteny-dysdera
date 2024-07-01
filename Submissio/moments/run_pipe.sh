#! /bin/bash -x
#$ -cwd
#$ -V

## usage: bash run_pipe.sh (without 'qsub')

# Command-line parameters.
QSUB_OPTIONFILE="$1"
PYT_PARAM_FILE="$2"
# Check command-line arguments were given.
if [ -z "$QSUB_OPTIONFILE" ] || [ -z "$PYT_PARAM_FILE" ] ; then
	echo "The arguments have not been specified"
	exit 1
fi

# Utilizem un 'environment' de Python amb els paquets necessaris preparats.
PYT_BIN="$HOME/anaconda3/envs/moments/bin/python"

i=0
cat $PYT_PARAM_FILE | while read line; do
	# Un nom de fitxer per a cada ordre qsub.
	STDOUT="$(date +%y%m%d-%H%M)-moments$i.log"
	qsub -@ "$QSUB_OPTIONFILE" -N "moments$i" -o "$STDOUT" -j yes -b y "$PYT_BIN" "$line"
	i=$(expr $i + 1)
done


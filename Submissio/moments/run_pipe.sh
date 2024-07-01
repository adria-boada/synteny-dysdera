#! /bin/bash -x
#$ -cwd
#$ -V

## usage: bash run_pipe.sh (without 'qsub')

source $HOME/../.bash_profile
# Activa el paquet de Python3 de 'moments'
conda activate moments
# Revisa la versió de Python3 i els paquets disponibles.
echo $(python3 --version)
echo $(python3 -c 'import moments')
# Podries fer un script que importi tots els paquets
# i n'asseguri el correcte funcionament.

QSUB_OPTIONFILE="$1"
PYT_PARAM_FILE="$2"

i=0
cat $PYT_PARAM_FILE | while read line; do
	# Un nom de fitxer per a cada ordre qsub.
	STDOUT="$(date +%y%m%d-%H%M)-moments.log"
	qsub -@ "$QSUB_OPTIONFILE" -N "moments$i" -o "$STDOUT" -j yes -b y python3 "$line"
	i=$(expr $i + 1)
done


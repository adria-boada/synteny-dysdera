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
STDOUT="$(date +%y%m%d-%H%M)-moments.log" 

cat $PYT_PARAM_FILE | while read line; do
	qsub -@ "$QSUB_OPTIONFILE" -N "moments" -o "$STDOUT" -j yes -b y python3 "$line"
done

#!/bin/bash

# Imprimeix els treballs enviats per qsub amb info util general.

# Espia als co-workers amb la opcio "-u" de qstat:
# Si el primer argument "$1" es un coworker, el busca.
# Per buscar tots els treballs de tots els usuaris, input="*" (incloent cometes)

if [ -z "$1" ]; then  # Si "$1" és zero (cadena buida, '').
	# Recull l'usuari que llança la ordre si el primer argument ($1) és buit. 
	user=$(whoami)
else
	user="$1"
fi

# Si no hi han treballs, no imprimeixis l'espai inicial
if [ -z "$(qstat -u "$user")" ]; then
	# si qstat es buit no facis res.
	true
else
	# Espai estetic i barra '===' inicial:
	#~ echo "" # [INHABILITAT]
	printf "%80s" | tr " " "="; echo ""
fi


# Per a tota la llista de treballs 
# (extreta de la ordre entre accents oberts, `qstat ... -f1`,
# que seria el mateix que invocar un subshell mitjançant '$()' ),
# extreu-ne la informacio rellevant:

for i in `qstat -u "$user" | tail -n+3 | cut -d " " -f1`; do

	# Imprimeix info general:
	qstat -u "$user" -j "$i" | grep "job_number"
	qstat -u "$user" -j "$i" | grep "job_name"
	qstat -u "$user" -j "$i" | grep "submission_time"

	# Atrapa l'estat que mostra qstat (running...).
	# [FALTA]

	# Atrapa la data de submissio del treball.
	SUBM_DATE=$(qstat -u "$user" -j "$i" | grep "submission_time" | cut -d" " -f13-18)

	# Transforma les dates a segons, per poder fer-ne la diferencia facilment:
	d1=$(date -d"$SUBM_DATE" +%s)
	d2=$(date +%s)

	# Diferencia entre data de subm. i data actual (en hores; dividit per 3600):
	echo "Hours since submission:     $(( ( $d2-$d1 )/3600 )) hours"

	# Atrapa'n el temps de CPU i memoria.
	qstat -u "$user" -j "$i" | grep "usage.*cpu.*vmem"

	# A quin node i quants processadors reservats
	# Sembla complicat pero el bloc nomes fa que quedi recte/justificat
	title_parallel=$(qstat -u "$user" -j "$i" | grep "parallel" | cut -d" " -f1-2)
	line_parallel=$(qstat -u "$user" -j "$i" | grep "parallel" | cut -d" " -f3-6)
	echo "$title_parallel      $line_parallel"

	# Arguments enviats
	qstat -u "$user" -j "$i" | grep "args"

	# Shell script usat
	qstat -u "$user" -j "$i" | grep "script_file:"

	# Barres entre treballs:
	printf "%80s" | tr " " "="; echo ""

done

echo # Línia buida final estètica

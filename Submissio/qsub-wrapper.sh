#! /bin/bash

# Parse arguments and options.
while [ "$1" != "" ]; do
	case $1 in
		# Number of threads requested.
		-t ) shift; ncpu="$1" ;;
		# Name given to the job submitted.
		-N ) shift; NAME="$1" ;;
		# Node {1, 8..13} where the job is sent.
		-n ) shift; NODE="$1" ;;
		# The rest of the args (script, program, options for them...)
		* )		ARGS+="$1 "
	esac
	shift
done
# Remove the ending whitespace "$1 "...
ARGS=${ARGS:0: -1}

# Quit and echo instructions if a missing arg. is found.
	# la variable `${0##*/}` crida el nom de l'script ($0=qsub_wrapper.sh)
	# sense afegir les barres i els directoris ('home\qsub_wrapper.sh').
for i in "$ncpu" "$NAME" "$NODE" ; do
	if [ -z "$i" ] ; then 
		echo "A variable has not been set properly"
		echo 'Call the script as follows:'
		echo "${0##*/} -N 'job-name' -n 'node[1; 8..13]' -t 'threads' 'SCRIPT-and-ARGS'"
		echo # move to a new line
		echo "Example with shellscript (keep in mind chmod u+x and shebang):"
		echo "${0##*/} -N prova -n 8 -t 2 script.sh all the other script args"
		echo "Example with binary (only to -n[ode] 1):"
		echo "${0##*/} -N find -n 1 -t 1 find . -size +1G -print"
		exit 1
	fi 
done

# Detect to which hostname was the query submitted to:
# Read single digits and translate to cluster addresses.
if [[ $NODE -eq 8 ]]; then
	Q_NODO="h0809.q"
	L_NODO="hostname=hercules08"
	THREADS="ompi24"
	INFO_NODE="Enviat a hercules08, a la cua ompi24 amb $ncpu CPUs"
elif [[ $NODE -eq 9 ]]; then
	Q_NODO="h0809.q"
	L_NODO="hostname=hercules09"
	THREADS="ompi24"
	INFO_NODE="Enviat a hercules09, a la cua ompi24 amb $ncpu CPUs"
elif [[ $NODE -eq 10 ]]; then
	Q_NODO="h10.q"
	L_NODO="hostname=hercules10"
	THREADS="ompi64h10"
	INFO_NODE="Enviat a hercules10, a la cua ompi64h10 amb $ncpu CPUs"
elif [[ $NODE -eq 11 ]]; then
	Q_NODO="h11.q"
	L_NODO="hostname=hercules11"
	THREADS="ompi64h11"
	INFO_NODE="Enviat a hercules11, a la cua ompi64h11 amb $ncpu CPUs"
elif [[ $NODE -eq 12 ]]; then
	Q_NODO="h12.q"
	L_NODO="hostname=hercules12"
	THREADS="ompi128h12"
	INFO_NODE="Enviat a hercules12, a la cua ompi128h12 amb $ncpu CPUs"
elif [[ $NODE -eq 13 ]]; then
	Q_NODO="h13.q"
	L_NODO="hostname=hercules13"
	THREADS="ompi255h13"
	INFO_NODE="Enviat a hercules13, a la cua ompi255h13 amb $ncpu CPUs"
elif [[ $NODE -eq 1 ]]; then
	Q_NODO="h0107.q"  # cas especial de la cua h0107...
	L_NODO="(auto)"   # El node tant fa, qui estigui lliure...
	THREADS="(auto)"  # Amb un thread, de sobres...
	ncpu="(1)"        # Amb un thread, de sobres...
	INFO_NODE="Enviat a la cua h0107.q"
else
	# Not a single case has been found to match...
	echo 'A number from 8 to 13 (representing a node) has not been given as the `-n` arg.'
	echo 'Exiting without launching a qsub'
	exit 1
fi

# Name for the stdout+stderr log-files, formatted as DATE-TITLE.log
# Prepends $(date) to 'job-name' ($NAME).
STDOUT="$(date +%y%m%d-%H%M)-${NAME}.log" 

# Demana si la ordre és correcta:
echo "La ordre que es llançarà és la següent (s'han eludit algunes opcions 'bàsiques' dins els punts suspensius):"
echo "qsub (...) -N $NAME -q $Q_NODO -l $L_NODO -pe $THREADS $ncpu -o $STDOUT $ARGS"
echo # move to a new line
read -p "Desitja continuar? [S/n]" -n 1 -r # La opció -n 1 deixa entrar tan sols un caràcter.
echo # move to a new line

# Si contesta positivament (S o Y), llança la ordre...
if [[ $REPLY =~ ^[YS]$ ]]; then
	# Omple output amb info.
	echo "==== Output pel treball '$NAME' ====" >> $STDOUT  # Nom, [id. de $JOB_ID]
	echo $(date) >> $STDOUT  # data completa amb `date`
	# Retalla 'hostname' de la variable $L_NODO; separa pel caràcter '=' en dues parts.
	echo "$INFO_NODE" >> $STDOUT  # Quin node, quantes cpus...
	echo >> $STDOUT  # Espai estètic
	echo "Ordre enviada:" >> $STDOUT
	echo "qsub -cwd -V -N $NAME -q $Q_NODO -l $L_NODO -pe $THREADS $ncpu -o $STDOUT -j yes $ARGS" >> $STDOUT
	echo "-- -- -- -- -- -- -- -- -- --" >> $STDOUT
	echo >> $STDOUT  # Espai estètic

	if [[ $NODE != 1 ]]; then  # Si NO envies a la cua h0107.q ...
		qsub -cwd -V -N "$NAME" -q "$Q_NODO" -l "$L_NODO" -pe "$THREADS" "$ncpu" -o "$STDOUT" -j yes $ARGS
	else  # Si envies a la cua h0107.q, comanda qsub diferent...
		qsub -cwd -V -N "$NAME" -q "$Q_NODO" -o "$STDOUT" -j yes $ARGS
	fi

	# La opció '-h' llança la ordre en espera... no comença fins a sotmetre `qalter -h U "$j_id" `.
	# Recollir job_id de l'últim treball enviat: `qstat | tail -n 1 | cut -f1 -d' '`
	#~JOB_ID=$(qstat | tail -n 1 | cut -f1 -d' ')

	# Ara, comença el treball
	#~qalter -h U "$JOB_ID"

# Si no esta satisfet amb la ordre i no prem 'S', no fa res.
else
	echo "S'ha anul·lat el llançament de qsub"
fi


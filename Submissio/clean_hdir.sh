#!/bin/bash

# Moure la resta d'outputs a la carpeta de "Job_sub/"

mv *.[oe][0-9][0-9][0-9][0-9][0-9][0-9][0-9]* ~/Job_sub/
mv *.p[oe][0-9][0-9][0-9][0-9][0-9][0-9][0-9]* ~/Job_sub/

# Elimina fitxers buits a dins de la carpeta especificada,
# que hauria de ser `Job_sub/`.
# Es pot canviar noms de carpetes per '.'.
# Buscar fitxers buits. L'opcio -delete elimina les troballes.

find Job_sub/ -maxdepth 1 -type f -empty -print -delete

#!/bin/bash

# Ordre associada a l'script que parseja un fitxer `.paf` : paf_parser.py
# Redirecciona output amb `>>` cap a un fitxer de sortida.

# List of paf-files or single paf-file: `ls paf-file1 paf-fileN >> list` or path/to/single/paf-file.
paf_file_or_list=$1


# En el cas de que 1 sigui un fitxer paf, executa segons el 1er paràgraf.
# Alternativament, interpreta que és una llista i corre el 2on paràgraf.

case "$1" in
  *.paf)
      # 1er arg és paf-file.
      i=$1
      echo "^file   ${i##*/}" ; paf_parser.py "$i"
      ;;

  *)
      # 1er arg és llista
      for i in $(cat "$1") ; do echo "^file   ${i##*/}" ; paf_parser.py "$i" ; done
      ;;

esac


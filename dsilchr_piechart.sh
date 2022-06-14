
#! bin/bash

repeat_file="/home/adria/Escritorio/Data/RepeatMasker/Dsil_V2.3.fasta.renamed2_RMclasses.out"
idx_file="/home/adria/Escritorio/Data/RepeatMasker/Dsil_V2.3.fasta.idx"

# per a tots els cromosomes principals:
for crm in DsilChr1 DsilChr2 DsilChr3 DsilChr4 DsilChr5 DsilChr6 DsilChru1 DsilChru2 DsilChrx Scaffold; do
	# create a table file
	python3 ~/Dipogits/synteny-dysdera/repeatmasker_parser3.py "$repeat_file" dsil "$crm" > tbl
	# show the table to the terminal (least its head:)
	head tbl
	# read the table file with piechart.py and create a piechart.
	python3 ~/Dipogits/synteny-dysdera/repeatmasker_piechart.py tbl -i "$idx_file" -c "$crm"
	# move the output
	mv pchart.png "${crm}_piechart.png"
done

# Rename scaffold so they do not overwrite between silvatica's and catalonica's.
mv Scaffold_piechart.png DcatScff_piechart.png

# remove tmp tble files...
rm tbl


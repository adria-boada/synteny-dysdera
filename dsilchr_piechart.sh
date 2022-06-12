
#! bin/bash

repeat_file="/home/adria/Dsil_V2.3_genome.fasta_renamed.out"
idx_file="/home/adria/Dsil_V2.3_genome.fasta.idx"

# per a tots els cromosomes principals:
for crm in Dsilchr1 Dsilchr2 Dsilchr3 Dsilchr4 Dsilchr5 Dsilchr6 Dsilchru1 Dsilchru2 Dsilchrx Scaffold; do
	# create a table file
	python3 ~/Dipogits/synteny-dysdera/repeatmasker_parser3.py "$repeat_file" dsil "$crm" > tbl
	# show the table to the terminal (least its head:)
	head tbl
	# read the table file with piechart.py and create a piechart.
	python3 ~/Dipogits/synteny-dysdera/repeatmasker_piechart.py tbl -i "$idx_file" -c "$crm"
	# move the output
	mv pchart.png "${crm}_piechart.png"
done

# remove tmp tble files...
rm tbl


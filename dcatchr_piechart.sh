
#! bin/bash

repeat_file="$1"
idx_file="$2"

# per a tots els cromosomes principals:
for crm in DcatChr1 DcatChr2 DcatChr3 DcatChr4 DcatChr5 DcatChrx Scaffold; do
	# create a table file
	python3 ~/Dipogits/synteny-dysdera/repeatmasker_parser3.py "$repeat_file" dcat35 "$crm" > tbl
	# show the table to the terminal (least its head:)
	head tbl
	# read the table file with piechart.py and create a piechart.
	python3 ~/Dipogits/synteny-dysdera/repeatmasker_piechart.py tbl -i "$idx_file" -c "$crm"
	# move the output
	mv pchart.png "${crm}_piechart.png"
done

# Rename scaffold so they do not overwrite between silvatica's and catalonica's.
mv Scaffold_piechart.png DsilScff_piechart.png

# remove tmp tble files...
rm tbl



#! bin/bash

repeat_file="~/"
idx_file=""

# per a tots els cromosomes principals:
for crm in dcatchr1 dcatchr2 dcatchr3 dcatchr4 dcatchr5 dcatchrx; do
	# create a table file
	python3 ~/Dipogits/synteny-dysdera/repeatmasker_parser3.py "$repeat_file" dsil "$crm" > tbl
	# show the table to the terminal (least its head:)
	head tbl
	# read the table file with piechart.py and create a piechart.
	python3 ~/Dipogits/synteny-dysdera/repeatmasker_piechart.py tbl -i "$idx_file" -c "$crm"
	# move the output
	mv piechart.png "${crm}_piechart.png"
done

# remove tmp tble files...
rm tbl


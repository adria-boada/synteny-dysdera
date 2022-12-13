#! /bin/bash

# El primer argument és el filename.stats.tmp
# que prové de juntar les seguents comandes:
#~ samtools stats BAM >> filename.stats.tmp
#~ samtools flagstats BAM >> filename.stats.tmp
fn=$1

# Filtra els stats generals dins el fitxer .stats.tmp
# i extreu-ne les dades d'interès.
qr=$(grep "^SN" $fn|cut -f2-|grep "raw total sequences")
echo "+ Number of raw reads queried: $(echo $qr|cut -f4 -d' ')"
prim=$(grep "primary mapped" $fn|cut -f1 -d' '|tr -d '(')
echo "+ Primary alignments: $prim"
nprim=$(grep "^SN" $fn|cut -f2-|grep "non-primary alignments")
echo "+ Non-primary alignments: $(echo $nprim|cut -f3 -d' ')"
supp=$(grep "^SN" $fn|cut -f2-|grep "supplementary alignments")
echo "+ Supplementary alignments: $(echo $supp|cut -f3 -d' ')"
readlen=$(grep "^SN" $fn|cut -f2-|grep "total length:")
readlen=$(echo $readlen|cut -f3 -d' ')
echo "+ Raw reads' cumulative length: $readlen"
nuc=$(grep "^SN" $fn|cut -f2-|grep "(cigar):")
nuc=$(echo $nuc|cut -f4 -d' ')
echo "+ Bases mapped (more accurate): $nuc"
err=$(grep "^SN" $fn|cut -f2-|grep "error rate")
echo "+ Error rate: $(echo $err|cut -f3 -d' ')"
insert=$(grep "^SN" $fn|cut -f2-|grep "insert size")
echo "+ Insert size and its standard deviation: $(echo $insert|cut -f4 -d' ') \$\pm\$ $(echo $insert|cut -f9 -d' ')"
perc=$(grep "primary mapped" $fn|cut -f6 -d' '|tr -d '(')
echo "+ Percentage of primary reads mapped: $perc"

# Seria interessant calcular "Expected(Coverage)" per a mides
# genòmiques de 3 Gb i 1.5 Gb.
# E(cov.) = Raw-read-cumulative-len / Expected-gensize
# Empraria la variable $readlen.
est30=$(python3 -c "print(round(($readlen/3000000000), 4))")
est15=$(python3 -c "print(round(($readlen/1500000000), 4))")
echo ''
echo "Expected coverages (valid if input is not filtered by regions, e.g. BUSCOs):"
echo ''
echo "+ If Gsize == 3.0 Gb, Exp(Cvg.) = $est30"
echo "+ If Gsize == 1.5 Gb, Exp(Cvg.) = $est15"


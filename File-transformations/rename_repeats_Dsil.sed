
# Afegeix signe d'interrogació al guió baix "_".
# Detecta expressions amb guió o sense guió baix.
# Acaba amb un espai (expressió \s). Torna a
# afegir aquest espai durant la substitució???

# Ordenats /sintènicament/ (1+4, 5+6+U1, etc)
# Substitute scaffold IDs with chromosome names in column 5:
s/Scaffold_\?14804_\?HRSCAF_\?18385\s/Dsilchr1 /
s/Scaffold_\?15272_\?HRSCAF_\?19649\s/Dsilchr4 /
s/Scaffold_\?3225_\?HRSCAF_\?3759\s/Dsilchr5 /
s/Scaffold_\?15361_\?HRSCAF_\?19822\s/Dsilchr6 /
s/Scaffold_\?13616_\?HRSCAF_\?15698\s/DsilchrU1 /
s/Scaffold_\?15005_\?HRSCAF_\?18897\s/Dsilchr3 /
s/Scaffold_\?14178_\?HRSCAF_\?16784\s/Dsilchr2 /
s/Scaffold_\?15362_\?HRSCAF_\?19823\s/DsilchrX /
s/Scaffold_\?15360_\?HRSCAF_\?19821\s/DsilchrU2_chrX /


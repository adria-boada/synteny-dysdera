
BEGIN{OFS="\t"}

# if there is a sixteenth field, replace it by "True"
# if there is not, create a sixteenth field containing "False"
{if ( $16 ) { $16="True" } else { $16="False" } }

# Afegeix signe d'interrogació al guió baix "_".
# Detecta expressions amb guió o sense guió baix.

# Ordenats /sintènicament/ (1+4, 5+6+U1, etc)
$5~/Scaffold_?14804_?HRSCAF_?18385/ {$5="Dsilchr1"}
$5~/Scaffold_?15272_?HRSCAF_?19649/ {$5="Dsilchr4"}
$5~/Scaffold_?3225_?HRSCAF_?3759/ {$5="Dsilchr5"}
$5~/Scaffold_?15361_?HRSCAF_?19822/ {$5="Dsilchr6"}
$5~/Scaffold_?13616_?HRSCAF_?15698/ {$5="DsilchrU1"}
$5~/Scaffold_?15005_?HRSCAF_?18897/ {$5="Dsilchr3"}
$5~/Scaffold_?14178_?HRSCAF_?16784/ {$5="Dsilchr2"}
$5~/Scaffold_?15362_?HRSCAF_?19823/ {$5="DsilchrX"}
$5~/Scaffold_?15360_?HRSCAF_?19821/ {$5="DsilchrU2_chrX"}

1


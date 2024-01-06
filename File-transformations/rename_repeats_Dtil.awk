
BEGIN{OFS="\t"}

# if there is a sixteenth field, replace it by "True"
# if there is not, create a sixteenth field containing "False"
{if ( $16 ) { $16="True" } else { $16="False" } }

# Afegeix signe d'interrogació al guió baix "_".
# Detecta expressions amb guió o sense guió baix.
$5~/Scaffold_?1/ {$5="DtilchrX"}
$5~/Scaffold_?2/ {$5="Dtilchr2"}
$5~/Scaffold_?3/ {$5="Dtilchr1"}
$5~/Scaffold_?4/ {$5="Dtilchr3"}
$5~/Scaffold_?5_?Vmt/ {$5="Dtilchr4"}
$5~/Scaffold_?6/ {$5="Dtilchr5"}
$5~/Scaffold_?7/ {$5="Dtilchr6"}

1


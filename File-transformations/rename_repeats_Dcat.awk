
BEGIN{OFS="\t"}

# if there is a sixteenth field, replace it by "True"
# if there is not, create a sixteenth field containing "False"
{if ( $16 ) { $16="True" } else { $16="False" } }

# Afegeix signe d'interrogació al guió baix "_".
# Detecta expressions amb guió o sense guió baix.
$5~/Scaffold_?1_?Vmt/ {$5="Dcatchr2"}
$5~/Scaffold_?2/ {$5="Dcatchr1"}
$5~/Scaffold_?3/ {$5="Dcatchr3"}
$5~/Scaffold_?4/ {$5="Dcatchr4"}
$5~/Scaffold_?5/ {$5="DcatchrX"}

1


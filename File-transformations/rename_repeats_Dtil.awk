
BEGIN{OFS="\t"}

# if there is a sixteenth field, replace it by "True"
# if there is not, create a sixteenth field containing "False"
{if ( $16 ) { $16="True" } else { $16="False" } }

$5=="Scaffold_1" {$5="DtilchrX"}
$5=="Scaffold_2" {$5="Dtilchr2"}
$5=="Scaffold_3" {$5="Dtilchr1"}
$5=="Scaffold_4" {$5="Dtilchr3"}
$5=="Scaffold_5_Vmt" {$5="Dtilchr4"}
$5=="Scaffold_6" {$5="Dtilchr5"}
$5=="Scaffold_7" {$5="Dtilchr6"}

1


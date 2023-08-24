
BEGIN{OFS="\t"}

{if ( $16 ) { $16="True" } else { $16="False" } }

$5=="Scaffold_1_Vmt" {$5="Dcatchr2"}
$5=="Scaffold_2" {$5="Dcatchr1"}
$5=="Scaffold_3" {$5="Dcatchr3"}
$5=="Scaffold_4" {$5="Dcatchr4"}
$5=="Scaffold_5" {$5="DcatchrX"}

1



BEGIN{OFS="\t"}

$1=="Scaffold_1" {$1="DtilchrX"}
$1=="Scaffold_2" {$1="Dtilchr2"}
$1=="Scaffold_3" {$1="Dtilchr1"}
$1=="Scaffold_4" {$1="Dtilchr3"}
$1=="Scaffold_5_Vmt" {$1="Dtilchr4"}
$1=="Scaffold_6" {$1="Dtilchr5"}
$1=="Scaffold_7" {$1="Dtilchr6"}

$6=="Scaffold_1_Vmt" {$6="Dcatchr2"}
$6=="Scaffold_2" {$6="Dcatchr1"}
$6=="Scaffold_3" {$6="Dcatchr3"}
$6=="Scaffold_4" {$6="Dcatchr4"}
$6=="Scaffold_5" {$6="DcatchrX"}

1


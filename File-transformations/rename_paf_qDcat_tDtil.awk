
BEGIN{OFS="\t"}

$6=="Scaffold_1" {$6="DtilchrX"}
$6=="Scaffold_2" {$6="Dtilchr2"}
$6=="Scaffold_3" {$6="Dtilchr1"}
$6=="Scaffold_4" {$6="Dtilchr3"}
$6=="Scaffold_5_Vmt" {$6="Dtilchr4"}
$6=="Scaffold_6" {$6="Dtilchr5"}
$6=="Scaffold_7" {$6="Dtilchr6"}

$1=="Scaffold_1_Vmt" {$1="Dcatchr2"}
$1=="Scaffold_2" {$1="Dcatchr1"}
$1=="Scaffold_3" {$1="Dcatchr3"}
$1=="Scaffold_4" {$1="Dcatchr4"}
$1=="Scaffold_5" {$1="DcatchrX"}

1



BEGIN{OFS="\t"}

$4=="Dcat" && $14=="Scaffold_1_Vmt" {$14="Dcatchr2"}
$4=="Dcat" && $14=="Scaffold_2" {$14="Dcatchr1"}
$4=="Dcat" && $14=="Scaffold_3" {$14="Dcatchr3"}
$4=="Dcat" && $14=="Scaffold_4" {$14="Dcatchr4"}
$4=="Dcat" && $14=="Scaffold_5" {$14="DcatchrX"}

$4=="Dtil" && $14=="Scaffold_1" {$14="DtilchrX"}
$4=="Dtil" && $14=="Scaffold_2" {$14="Dtilchr2"}
$4=="Dtil" && $14=="Scaffold_3" {$14="Dtilchr1"}
$4=="Dtil" && $14=="Scaffold_4" {$14="Dtilchr3"}
$4=="Dtil" && $14=="Scaffold_5_Vmt" {$14="Dtilchr4"}
$4=="Dtil" && $14=="Scaffold_6" {$14="Dtilchr5"}
$4=="Dtil" && $14=="Scaffold_7" {$14="Dtilchr6"}

1

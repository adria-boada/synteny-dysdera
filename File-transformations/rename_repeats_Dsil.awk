
BEGIN{OFS="\t"}

# if there is a sixteenth field, replace it by "True"
# if there is not, create a sixteenth field containing "False"
{if ( $16 ) { $16="True" } else { $16="False" } }

# ordenats /sintènicament/ (1+4, 5+6+U1, etc)
$5=="Scaffold_14804_HRSCAF_18385" {$5="Dsilchr1"}
$5=="Scaffold_15272_HRSCAF_19649" {$5="Dsilchr4"}
$5=="Scaffold_3225_HRSCAF_3759" {$5="Dsilchr5"}
$5=="Scaffold_15361_HRSCAF_19822" {$5="Dsilchr6"}
$5=="Scaffold_13616_HRSCAF_15698" {$5="DsilchrU1"}
$5=="Scaffold_15005_HRSCAF_18897" {$5="Dsilchr3"}
$5=="Scaffold_14178_HRSCAF_16784" {$5="Dsilchr2"}
$5=="Scaffold_15362_HRSCAF_19823" {$5="DsilchrX"}
$5=="Scaffold_15360_HRSCAF_19821" {$5="DsilchrU2_chrX"}

1


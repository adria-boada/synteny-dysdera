
# Afegeix signe d'interrogació al guió baix "_".
# Detecta expressions amb guió o sense guió baix.
# Acaba amb un espai (expressió \s). Torna a
# afegir aquest espai durant la substitució???

# Substitute scaffold IDs with chromosome names in column 5:
s/Scaffold_\?1_\?Vmt\s/Dcatchr2 /
s/Scaffold_\?2\s/Dcatchr1 /
s/Scaffold_\?3\s/Dcatchr3 /
s/Scaffold_\?4\s/Dcatchr4 /
s/Scaffold_\?5\s/DcatchrX /


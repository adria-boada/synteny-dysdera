
# Afegeix signe d'interrogació al guió baix "_".
# Detecta expressions amb guió o sense guió baix.
# Acaba amb un espai (expressió \s). Torna a
# afegir aquest espai durant la substitució???

# Substitute scaffold IDs with chromosome names in column 5:
s/Scaffold_\?1\s/DtilchrX /
s/Scaffold_\?2\s/Dtilchr2 /
s/Scaffold_\?3\s/Dtilchr1 /
s/Scaffold_\?4\s/Dtilchr3 /
s/Scaffold_\?5_\?Vmt\s/Dtilchr4 /
s/Scaffold_\?6\s/Dtilchr5 /
s/Scaffold_\?7\s/Dtilchr6 /


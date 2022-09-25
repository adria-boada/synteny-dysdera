#! bin/bash

# $1=taula de repeticions de silvatica, per transformar-la a un format comú,
# compartit amb D. catalonica.

sed -e '/Pararetrovirus/{s/$/@/}' -e '/Retrotransposon/{s/$/@/}' -e '/\*macro\*/{s/$/@/}' "$1" | tr '@' '\n'


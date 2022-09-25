#! bin/bash

# $1=taula de repeticions de catalonica, per transformar-la a un format comú,
# compartit amb D. silvatica.

sed -e '/L1-dep/{s/$/@/}' "$1" | tr '@' '\n'


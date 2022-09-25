#!/usr/bin/awk -f

# Script awk per trobar línies que tinguin una
# columna entre una parella de ranges.
(200 <= $2 && $2 <= 600) {print}


#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# paflines_into_fastaexcerpts.py
#
# 13 de juny 2022  <adria@molevol-OptiPlex-9020>

"""Takes in a few lines from a paf-file, tracks the coordinates
and chromosome name for the query and target on each line,
and tries to spot those fragments inside a pair of query and
target fasta files.
"""

import sys

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Si es crida com a script:
if __name__ == '__main__':
    ### VALORS D'ARGUMENTS A LA TERMINAL ###
    import argparse

#    parser = argparse.ArgumentParser(description='')
#    # file-name: positional arg.
#    parser.add_argument('filename', type=str, help='... file-name')
#    # integer argument
#    parser.add_argument('-a', '--numero_a', type=int, help='Paràmetre "a"')
#    # choices argument
#    parser.add_argument('-o', '--operacio', 
#            type=str, choices=['suma', 'resta', 'multiplicacio'],
#            default='suma', required=False,
#            help='')
#
#    args = parser.parse_args()
#    # call a value: args.operacio or args.filename.
    


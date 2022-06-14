#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# coordlines_into_fastaexcerpts.py
#
# 14 jun 2022  <adria@molevol-OptiPlex-9020>

"""Takes in a few lines from a coordinate file, tracks the regions and
chromosome name for each line, and tries to spot those fragments inside a fasta
file.

OUTPUTS TO STDOUT.
"""

import sys
import paflines_into_fastaexcerpts as fex

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Si es crida com a script:
if __name__ == '__main__':
    ### VALORS D'ARGUMENTS A LA TERMINAL ###
    import argparse

    parser = argparse.ArgumentParser(
            description='Coordinate extraction from FASTA')
    # file-name: positional arg.
    parser.add_argument('fasta', type=str, help='FASTA file-path')
    # file-name: positional arg.
    parser.add_argument('coordinates', type=str, help='coordinate file-path')

    args = parser.parse_args()
    # call a value: args.operacio or args.filename.
    

    # Acquire regions of interesting coordinates...
    # from the specified coordinate file.
    with open(args.coordinates) as cc:
        # init dict, which will save pairs of header patterns
        # and lists of ranges of seq inside that pattern.
        dc = {}
        # Jump the header of the table of coordinates:
        capçalera_coords = next(cc)

        for line in cc:
            if line == '\n':
                # Si arriba al final del fitxer, apaga màquines.
                break

            # the coord file should look like:
            pattern, begin, end, tag = line.strip('\n').split()
            # int-ize coordinates:
            begin, end = int(begin), int(end)
            # si no s'ha creat una entrada per aquest header-pattern:
            if not pattern in dc:
                # init key.
                dc[pattern] = []
            # add the list of coords.
            dc[pattern] += [ [begin, end, tag] ]


    # Acquire and print to stdout the sequences:
    fex.fasta_singleline_excerpts(args.fasta, dc)




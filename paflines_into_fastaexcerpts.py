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
import pafparser as pfp

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Si es crida com a script:
if __name__ == '__main__':
    ### VALORS D'ARGUMENTS A LA TERMINAL ###
    import argparse

    parser = argparse.ArgumentParser(description='Takes in a few lines of a paf-file, \
            tracks its coordinates and chromosome names for the target and query pairs \
            on each line, and tries to spot those fragments inside the corresponding \
            fasta files.')
    # file-name: positional arg.
    parser.add_argument('paflines', type=str, help='Path to a PAF-file with selected lines.')
    # 'optional' arguments
    parser.add_argument('query', type=str, 
            help='Query single-line-FASTA file-name. Scaffold names compatible with PAF-file.')
    parser.add_argument('target', type=str, 
            help='Target single-line-FASTA file-name. Scaffold names compatible with PAF-file.')
#    # choices argument
#    parser.add_argument('-o', '--operacio', 
#            type=str, choices=['suma', 'resta', 'multiplicacio'],
#            default='suma', required=False,
#            help='')

    args = parser.parse_args()
    # call a value: args.operacio or args.filename.
    
    paf_file = args.paflines
    fasta_target = args.query
    fasta_query = args.target

    pafnum=0
    with open(paf_file) as paf:
        for pafline in paf:
            # Create a pafnum to keep in mem the number of the
            # actual line, which will be appended to the outputted
            # resulting fasta (in order to pair sequences from query
            # and target of the same line). 
            pafnum+=1

            # Create a Mapping object, easier to handle.
            obj = pfp.Mapping(pafline)
            # begin coordinates for paf-ali:
            # begin for paf-files are zero-based, no need to modify.
            bt, bq = obj.tstart, obj.qstart
            # end coordinates for paf-ali:
            # add +1 because python slices exclude end.
            et, eq = obj.tend+1, obj.qend+1
            # name of query and target:
            nt, nq = obj.tname, obj.qname
            # Zip target and query fasta with their coordinates:
            # example:
            #~for coord, fasta in ((coord1, fasta1), (coord2, fasta2))
            for coord, fasta in (([bt, et, nt], fasta_query), ([bq, eq, nq], fasta_target)):
                # open each fasta file
                with open(fasta) as fn:
                    while True:
                        try:
                            # Extract a couple of contiguous lines:
                            # The first one should be header, the second the seq.
                            # does not .strip('\n') because the later print takes this into consideration
                            cap, seq = next(fn), next(fn)
                        # Si s'arriba al final del fitxer, next(cap) donarà
                        # l'error 'StopIteration' i es trencarà el while loop.
                        except StopIteration:
                            break
                        # Si la capçalera conté el patró d'interés, analitza:
                        if coord[2] in cap:
                            print(cap,'_',pafnum, sep='', end='')
                            print(seq[coord[0]:coord[1]], end='')



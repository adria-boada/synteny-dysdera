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


def unique_paftag_list(file: "File with alignment, PAF formatted."):
    """Inputs a PAF file and outputs the query and
    target scaffold names, in a dict, with unique entries.
    """
    out = {'q':[], 't':[]}

    with open(file) as fn:
        for line in fn:
            obj = pfp.Mapping(line)
            if not obj.qname[2:] in out['q']:
                out['q'] += [obj.qname[2:]]
            if not obj.tname[2:] in out['t']:
                out['t'] += [obj.tname[2:]]

    return out


def fasta_singleline_excerpts(filepath: "Path to a FASTA file", 
        dc: "Dictionary with chr keys and coord-list values"):
    """
    This function opens a FASTA and tries to find the regions
    denoted by the 'dict_coord' input variable. Once found,
    prints them in stdout.

    INPUT
    =====

    The first var is the path to a FASTA file from which we want
    to acquire specific regions.

    !!! MUST BE SINGLELINE FASTA !!!

    single-line FASTA format:
     * NO starting blank lines. First line is a header.
     * NO blank lines between sequences.
     * Headers followed by sequences. If a sequence or
     header is unpaired, it will screw up the rest of the
     analysis. 

    The second var ('dc') is a dictionary with header-patterns and
    coordinate-list pairs. The keys should be found inside the header, and the
    value for each key is a list of coordinates to extract from that
    header-seq.

    Example: Find the patterns 'dsilchr1' and 'dsilchr2' inside the header, and
    extract begin--end regions from the following sequence . If there is a tag,
    append it to the end of the copy-pasted header from the file.

    dc = {'dsilchr1': [[begin, end, tag], (...), [begin, end, tag]],
        'dsilchr2': [[begin, end, tag], (...), [begin, end, tag]] 
        (...): [(...)]
    }

    OUTPUT
    ======

    Prints to stdout pairs of header-sequence. 

    Use pipes and dashes to save into a new fasta file: '|', '>', '>>'.
    """

    # open the FASTA file:
    with open(filepath) as fasta:
        while True:
            try:
                # Extract a couple of contiguous lines:
                # The first one should be header, the second the seq.
                # does not .strip('\n') because the later print takes this into consideration
                cap, seq = next(fasta).strip('\n'), next(fasta).strip('\n')

                # Si la capçalera no conté '>', raise error i surt de la funció.
                # voldria dir que el FASTA no és single-line. 
                if not '>' in cap:
                    raise ValueError("S'ha trobat una capçalera dins el FASTA sense '>'.")

            # Si arriba al final del fitxer i no aconsegueix extreure la
            # pròxima línia, trenca el while loop.
            except StopIteration:
                break

            # Busca si algun patró cromosòmic (claus del diccionari) encaixa:
            for pat in dc.keys():
                # Si el patró el trobes dins la capçalera del FASTA:
                if pat.casefold() in cap.casefold():  # cerca case-insensitive. 
                    # Comença a extreure la llista de regions:
                    # Haurien de tenir el format [begin, end, tag]
                    for coord in dc[pat]:
                        # imprimeix el header amb l'identificador.
                        # Si les regions contenen etiquetes, append them:
                        try:
                            # si existeix coord[2] (etiqueta), append.
                            print(cap, coord[2], sep='')
                        except IndexError:
                            # Si no, sols imprimeix el header.
                            print(cap)
                        # Imprimeix entre begin i end del pafalign dict.
                        print(seq[ coord[0]:coord[1] ])


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
    parser.add_argument('paflines', type=str, 
            help='Path to a PAF-file with selected lines.')
    parser.add_argument('query', type=str, 
            help='Query single-line-FASTA file-name. Scaffold names compatible with PAF-file.')
    parser.add_argument('target', type=str, 
            help='Target single-line-FASTA file-name. Scaffold names compatible with PAF-file.')

    args = parser.parse_args()
    # call a value: args.operacio or args.filename.
    
    # Init a dict which will store paf-lines:
    pafaligns = {'query':{}, 'target':{}}
        
    # Init all the unique chromosomes:
    # for query:
    for crm in unique_paftag_list(args.paflines)['q']:
        pafaligns['query'][crm] = []
    # for target:
    for crm in unique_paftag_list(args.paflines)['t']:
        pafaligns['target'][crm] = []
        
    pafnum=0
    # Open the three files using a 'with' statement:
    with open(args.paflines) as paf:
        for line in paf:
            # Create a pafnum to keep in mem the number of the
            # actual line, which will be appended to the outputted
            # resulting fasta (in order to pair sequences from query
            # and target of the same line). 
            pafnum+=1

            # Create a Mapping object, easier to handle.
            obj = pfp.Mapping(line)

            # begin coordinates for paf-ali:
            # begin for paf-files are zero-based, no need to modify.
            bt, bq = obj.tstart, obj.qstart
            # end coordinates for paf-ali:
            # add +1 because python slices exclude end.
            et, eq = obj.tend+1, obj.qend+1
            # name of query and target:
            # s'han d'eliminar les etiquetes identificadores que posa la classe
            # 'Mapping()':  # total        etiqueta   nom
                            # T.DsilChrX --> 'T.' + 'DsilChrX'
            nt, nq = obj.tname[2:], obj.qname[2:]

            # computer, do the dance.
            # The same id. pafnum is added to both query and target from the same line.
            pafaligns['query'][nq] += [[bq, eq, pafnum]]
            pafaligns['target'][nt] += [[bt, et, pafnum]]

    # Busca a dins de cada FASTA les regions d'interès del PAF, 
    # tant per a target com per a query:
    for file, dict_chr_coord in ((args.query, pafaligns['query']), 
            (args.target, pafaligns['target'])):
        # The following function prints in stdout:
        fasta_singleline_excerpts(file, dict_chr_coord)





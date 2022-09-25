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
import repeatmasker_parser3 as rm3

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def unique_pafchr_list(file: "File with an alignment, formatted as a PAF."):
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


def repeatmasker_finder(filepath: "Path to a FASTA file",
        dc: "Dictionary with chr keys and coord-list values"):
    """
    """

    with open(filepath) as rm:
        while True:
            try:
                # Extract a single line of the rm file.
                line = next(rm).strip('\n')
                # Create a Repeat() obj.
                obj = rm3.Repeat(line)

            except StopIteration:
                break

            # Busca si algún patró cromosòmic encaixa dins la línia RM:
            for pat in dc.keys():
                # si rm.refid encaixa amb pafline.crm:
                if pat.casefold() in obj.refid.casefold():
                    # Comprova que la línia es troba dins el rang de
                    # coordenades.
                    for coord in dc[pat]:
                        if (coord[0] < obj.end) and (coord[1] > obj.begin):
                            if obj.orient=='C': obj.orient='-'
                            # If every condition is satisfied, print(line).
                            print(pat+coord[2],
                                  # Extract the relative beggining and end of
                                  # the repetition (relative to paf-ali).
                                  obj.begin-coord[0], obj.end-coord[0],
                                  # Strand, orientation:
                                  obj.orient,
                                  # paf begin and end
                                  coord[0], coord[1],
                                  # Get all the classification attributes.
                                  obj.Class, obj.Subclass, obj.Order, obj.Super_family,
                                  sep="\t")
                    # Si amb el primer chr match no aconsegueix entrar dins el
                    # rang proposat per les coords de paflines, trenca el loop:
                    break


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Si es crida com a script:
if __name__ == '__main__':
    ### VALORS D'ARGUMENTS A LA TERMINAL ###
    import argparse, pathlib

    parser = argparse.ArgumentParser(description='Takes in a few lines of a paf-file, \
            tracks its coordinates and chromosome names for the target and query pairs \
            on each line, and tries to spot those fragments inside the corresponding \
            fasta files.')
    # file-name: positional arg.
    parser.add_argument('paflines', type=pathlib.Path,
            help='Path to a PAF-file with selected lines.')
    parser.add_argument('--fasta', type=pathlib.Path,
            help='Optional Single-line-FASTA file-name. Scaffold names compatible with \
                        PAF-file. Both query and target FASTA files should be \
                        concatenated into a single one before submission.',
                        default=argparse.SUPPRESS)
    #~parser.add_argument('target', type=str,
    #~        help='Target single-line-FASTA file-name. Scaffold names compatible with PAF-file.')
    parser.add_argument('--rmfile', type=pathlib.Path, required=False,
            help='Optional RM file from which to also extract repeats inside \
                        the ranges of the given PAF file.',
                        # If there is no given rmfile, supress and do not do
                        # analysis. Go to the end of the file to see the if statement.
                        default=argparse.SUPPRESS)

    args = parser.parse_args()
    # call a value: args.operacio or args.filename.

    # Init a dict which will store paf-lines:
    pafaligns = {'query':{}, 'target':{}}

    # Init all the unique chromosomes:
    # for query:
    for crm in unique_pafchr_list(args.paflines)['q']:
        pafaligns['query'][crm] = []
    # for target:
    for crm in unique_pafchr_list(args.paflines)['t']:
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
            # However, we're after surrounding region of DNA, so substract a few
            # kilobases of margin to look around.
            look_around = 5000
            bt, bq = obj.tstart-look_around, obj.qstart-look_around
            # end coordinates for paf-ali:
            # add +1 because python slices exclude end.
            # As we also would like to look around, adds kbs:
            et, eq = obj.tend+1+look_around, obj.qend+1+look_around
            # name of query and target:
            # s'han d'eliminar les etiquetes identificadores que posa la classe:
            # 'Mapping() Class': # default        tag    crm.name
                                 # T.DsilChrX --> 'T.' + 'DsilChrX'
            nt, nq = obj.tname[2:], obj.qname[2:]

            # computer, do the dance.
            # The same id. pafnum is added to both query and target from the same line.
            # (tagged at the end of the header of each FASTA sequence).
            pafaligns['query'][nq] += [[bq, eq,
                                        '_paf'+str(pafnum)+'_'+str(bq)+'+'+str(eq-bq)]]
            pafaligns['target'][nt] += [[bt, et,
                                         '_paf'+str(pafnum)+'_'+str(bt)+'+'+str(et-bt)]]

    # Busca a dins de cada FASTA les regions d'interès del PAF,
    # tant per a target com per a query:
    if "fasta" in args:
        for dict_chr_coord in (pafaligns['query'],
                 pafaligns['target']):
            # The following function prints in stdout:
            fasta_singleline_excerpts(args.fasta, dict_chr_coord)

    if "rmfile" in args:
        # print a header for Rstudio.
        # Also, a useful separator...
        print('refid', # name of the chromosome/scaffold
              'relative_start', 'relative_end', # Repeat pos relative to
              # mapping.
              'strand', # + (normal) or C (compli)
              'paf_start', 'paf_end', # Position inside chromosome
              'class', 'subclass', 'order', 'superfam', # classification.
              sep="\t")

        for dict_chr_coord in (pafaligns['query'],
                                 pafaligns['target']):
            repeatmasker_finder(args.rmfile, dict_chr_coord)



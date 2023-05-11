#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# broccoli_to_mcl.py
#
# 11 May 2023  <adria@molevol-OptiPlex-9020>

"""
Turn the output table of Broccoli (modified by laboratuvar arkadaşlarım) into
an MCL-like output: a txt file where each line is a family, containing all gene
IDs for that family.
"""

import sys
import pandas as pd

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class Broccoli:

    def __init__(self, file:"Path to the table that requires modification",
                 species=['Dcat', 'Dtil']):
        """
        For sp in species (Dcat, Dtil) and orthogroups in dataframe (OG1, OG2,
        ...OGN) subset the original table (inputted file) and print a row with
        all gene IDs in each family (unique combinations of species+OG).
        """
        # try to open the file (make sure provided path is correct)
        try:
            with open(file) as ft:
                pass
        except:
            sys.exit('ERROR: The file-path does not exist?')

        # read and prepare dataframe from tabulated file
        df = pd.read_table(file)
        families_out = {}
        for sp in species:
            print(f'Status: Families for species {sp}')
            for og in df['OGid'].unique():
                df_subsetted = df.query(f'Species == "{sp}" & OGid == "{og}"')
                families_out[og] = []
                # print the genes belonging to the subsetted
                # family (sp, og) in a tabulated format
                for gene in df_subsetted['GeneID']:
                    families_out[og] += [gene]
            # create a list of sublists
            # each sublist is a family (group of geneIDs)
            l = list(families_out.values())
            # sort by length of family (amount of geneIDs)
            l.sort(key=lambda x: len(x), reverse=True)
            for fam in l:
                for gene in fam[:-1]:
                    print(gene, end='\t')
                print(fam[-1])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Si es crida com a script:
if __name__ == '__main__':

    ### PARSE ARGUMENTS ###
    import argparse

    parser = argparse.ArgumentParser(description='')
    # file-name: positional arg.
    parser.add_argument('file', type=str, help='Path to "table" file')
    # integer argument
#    parser.add_argument('-a', '--numero_a', type=int, help='Paràmetre "a"')
    # choices argument
#    parser.add_argument('-o', '--operacio', 
#            type=str, choices=['suma', 'resta', 'multiplicacio'],
#            default='suma', required=False,
#            help='')

    args = parser.parse_args()
    # call a value: args.operacio or args.filename.
    Broccoli(args.file, species=['Dtil'])



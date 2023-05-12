#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# broccoli_to_mcl.py
#
# 11 May 2023  <adria@molevol-OptiPlex-9020>

"""
Turn the output table of Broccoli (modified by laboratuvar arkadaşlarım) into
an MCL-like output: a txt file where each line is an orthogroup, containing all gene
IDs for that orthogroup.

-- Example output --
gene1   gene2   gene3
gene4
gene5
--------------------
The first 3 genes comprise OG1, and lines 2 and 3 would contain orphan genes.
Genes should be separated by \t.
"""

import sys
import pandas as pd

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class Broccoli:

    def __init__(self, file:"File-path to the table that requires modification"):
        """
        1. Check if input `file` exists and `sys.exit()` if it does not.
        2. Create a pandas.DataFrame() from the input `file`

        -- TO FURTHER FILTER INPUT TABLE... --

        If you desire to filter the dataframe by types of orthogroups ("broccoli
        buenos", filtering by quality, filtering by RNAseq status...) it can be
        done right after reading from file with:
        ` df.query(f'Column == "value" | Column == "{variable}"') `

        f'': formatted-string, can introduce variables to strings within
        braces '{}'
        To combine logical expressions, use '|' (OR) and '&' (AND).
        """

        # try to open the file (make sure provided path is correct)
        try:
            with open(file) as ft:
                pass
        except:
            sys.exit('ERROR: The inputted file-path does not exist?')
        # read and prepare dataframe from tabulated file
        self.df_input_table = pd.read_table(file)
        # init df_dups; if it is later needed, it will be computed
        # and stored in this variable (overwritting the now empty DataFrame)
        self.df_dups = pd.DataFrame()


    def paranome_OG_dict(self, species,
            out: 'Outfile-path where to write results'=1):
        """
        Write the 'paranome.tsv.mcl' (one orthogroup per line, each line the set of
        genes for that OG separated by \t) onto a file.

        For `sp` in species list (Dcat, Dtil) and orthogroups in dataframe (OG1, OG2,
        ...OGN) subset the original table (inputted file) and print a row with
        all gene IDs in each family (unique combinations of species+OG).
        """

        # create outfile pathname string if it was not provided (default 1)
        if out==1:
            out = f'broquil_{species}_paranome.tsv.mcl'
        # check that the outfile does not exist (no file will be overwritten):
        try:
            with open(out) as o:
                sys.exit(f'ERROR: The given outfile path {out} '+
                     'already exists and would be overriden!')
        except:
            pass
        # shorten var calls
        df = self.df_input_table
        # store genes in dict. before writing to file:
        families_out = {}
        print(f'Status: Computing `paranome` for species {species}')
        # collect unique OGids present in 'species':
        for og in df.query(f'Species == "{species}"')['OGid'].unique():
            # create a list in which all genes from og orthogroup
            # will be stored
            families_out[og] = df.query(f'Species == "{species}" & '+
                                        f'OGid == "{og}"')['GeneID'].to_list()
        # write the dictionary to file:
        self.write_OG_dict(families_out, out)

        return (families_out, out)

    def one_v_one_orthologs_dict(self, sp_pair,
            out: 'Outfile-path where to write results'=1):
        """
        Write all N to N orthologs (where N != 0) for a pair of species.
        Just like if it were an MCL file (rows of geneIDs for each OG).
        """
        # create outfile string if it was not provided (default 1)
        if out==1:
            out = f'broquil_{sp_pair[0]}_{sp_pair[1]}_one_to_one.tsv.mcl'
        # check that the outfile does not exist and will not be overwritten:
        try:
            with open(out) as o:
                sys.exit(f'ERROR: The given outfile path {out} '+
                     'already exists and would be overriden!')
        except:
            pass
        # shorten var calls
        df = self.df_input_table
        # verify that the given sp_pair exists; species name correctly
        # spelled...
        for sp in sp_pair:
            # if `sp` not in dataframe:
            if not df[df['Species']==sp]['Species'].any():
                sys.exit("ERROR: A species given to the function "+
                         "`write_orthonome_ovo_mcl` "+
                         "has been mispelled (not found in dataframe)")
        # do we need to compute `tmp_df_amounts()` (yes if it is empty)?
        if self.df_dups.empty:
            self.df_dups = self.tmp_df_amounts()
        df_one_to_one = self.df_dups.query(f'{sp_pair[0]} == 1 & '+
                                           f'{sp_pair[1]} == 1')
        # store genes in dict. before writing to file:
        families_out = {}
        print(f'Status: Computing `orthonome` for species pair {sp_pair}')
        # collect unique OGids present in 'species':
        for og in df_one_to_one.index:
            # create a list in which all genes from og orthogroup
            # will be stored
            # The query could be sorted by species (so one always precedes)
            families_out[og] = df.query(f'(Species == "{sp_pair[0]}"   | '+
                        f'Species == "{sp_pair[1]}")   & '+
                        f'OGid == "{og}"')['GeneID'].to_list()
        # write the dictionary to file:
        self.write_OG_dict(families_out, out)

        return (families_out, out)

    def write_OG_dict(self, dinput, outname):
        """
        Write the given `dict` to `outname` filepath in a correctly tabulated
        manner.
        """

        # create a list of sublists with the generated dicts
        # each sublist in the list is a family (group of geneIDs)
        l = list(dinput.values())
        # sort by length of family (amount of geneIDs)
        l.sort(key=lambda x: len(x), reverse=True)
        # write to outfile:
        with open(outname, 'w') as o:
            for fam in l:
                for gene in fam[:-1]:
                    # sep each gene with \t
                    o.write(gene + '\t')
                # end line with \n
                o.write(fam[-1] + '\n')

        return None

    def tmp_df_amounts(self,
                       OGtype_species_order: "Order of species in OGtype col"=(
                           'Dcat', 'Dtil', 'Dsil', 'Dver', 'Dban', 'Dgom')):
        """
        Create a tmp. subset of a given dataframe. Returns it.

        One column with OGids. One column per given species, where the amount of
        dups for that sp-OG is stored.

        It is a very slow function (queries a massive dataframe, many times)
        """

        # shorten var calls
        df = self.df_input_table
        # initialise output dataframe
        df_dups = pd.DataFrame(0,
                               # one col per species
                               columns=OGtype_species_order,
                               # one row per OG
                               index=df['OGid'].unique())
        print(df_dups) #DEBUG
        print(df_dups.columns) #DEBUG
        for i in df_dups.index:
           for sp in df_dups.columns:
                df_dups.loc[i, sp] = len(df.query(f'OGid == "{i}" & Species == "{sp}"'))
        print(df_dups) #DEBUG

        return df_dups


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
    broc = Broccoli(args.file)
    broc.one_v_one_orthologs_dict(['Dcat', 'Dtil'])
    broc.paranome_OG_dict('Dcat')
    broc.paranome_OG_dict('Dtil')



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
import matplotlib.pyplot as plt, numpy as np

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

        To subset all OGids within a list...
        df[df['OGid'].isin(list_OGids)]
        """

        # try to open the file (make sure provided path is correct)
        try:
            with open(file) as ft:
                pass
        except:
            sys.exit('ERROR: The inputted file-path does not exist?')
        # read and prepare dataframe from tabulated file
        self.df_input_table = pd.read_table(file,
            dtype={ # specify dtypes of a few columns (faster file reading)
            'OGid': str,
            'Method': str,
            'Species': str,
            'BRAKER_status':str,
            'ChimericGeneStatus': bool,
            'RNAseq_status': str, 'RNAseq_OGstatus': str,
            'GO_OGstatus': str, 'GO_status': str,
            'GeneStart': int,
            'GeneEnd': int,
            })
        # init df_amounts; if it is later needed, it will be computed
        # and stored in this variable (overwritting the now empty DataFrame)
        self.df_amounts = pd.DataFrame()


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
                     'already exists and would be overwritten!')
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
        df_one_to_one = self.df_dups().query(f'{sp_pair[0]} == 1 & '+
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

    def df_dups(self,
                       OGtype_species_order: "Order of species in OGtype col"=(
                           'Dcat', 'Dtil', 'Dsil', 'Dver', 'Dban', 'Dgom')):
        """
        Create a tmp. subset of a given dataframe. Returns it.

        One column with OGids. One column per given species, where the amount of
        dups for that sp-OG is stored.

        It is a very slow function (queries a massive dataframe, many times)
        """

        # do we need to compute `df_dups()` (yes if self.df_dups is empty)
        if not self.df_amounts.empty:
            return self.df_amounts
        # shorten var calls
        df = self.df_input_table
        # initialise output dataframe
        df_dups = pd.DataFrame(0,
                               # one col per species
                               columns=OGtype_species_order,
                               # one row per OG
                               index=df['OGid'].unique())
        for i in df_dups.index:
           for sp in df_dups.columns:
                df_dups.loc[i, sp] = len(df.query(f'OGid == "{i}" & Species == "{sp}"'))
        self.df_amounts = df_dups

        return df_dups

    def count_minor_major_scaffold_occupancy(self, species,
            df: "It is possible to supply a filtered dataframe"= pd.DataFrame(),
            up_to: "Up to how many OGtypes included?"= 4):
        """
        Count the orthogroups present in only minor, only major, and mix of
        major-minor scaffolds. Group minor-OGs, major-OGs, and their
        intersection.
        """
        if len(df)==0:
            df = self.df_input_table
        if len(species)!=2:
            sys.exit('ERROR: supplied a list of species longer than 2')
        # df storing the recounting of OGs per OGtype...
        df_count_ogtypes = pd.DataFrame(
            columns=[species[0]+str(i) for i in range(1, up_to)]+[species[0]+f'>={up_to}'],
            index=[species[1]+str(i) for i in range(1, up_to)]+[species[1]+f'>={up_to}']
        )
        # compute the inner 1 to up_to-1 matrix
        for i in range(0, up_to):
            for j in range(0, up_to):
                OGid_interested_list = self.df_dups().query(
                    f'{species[0]} == {i} & {species[1]} == {j}'
                ).index.to_list()
                # subset DF for OGids with i and j amount of orthologs (OGtypes)
                d = df.query(f'OGid in {OGid_interested_list}')
                ##print(d[['OGid', 'Species', 'ScaffoldCategory', 'OGtype']])##DEBUG
                df_count_ogtypes.loc[species[1]+str(j), species[0]+str(i)] = (
                len(d['OGid'].unique()))

        # compute both i=(1,3) and j>=up_to for both species
        for i in range(0, up_to):
            OGid_interested_list = self.df_dups().query(
                f'{species[0]} == {i} & {species[1]} >= {up_to}'
            ).index.to_list()
            d = df.query(f'OGid in {OGid_interested_list}')
            ##print(d[['OGid', 'Species', 'ScaffoldCategory', 'OGtype']])##DEBUG
            df_count_ogtypes.loc[species[1]+'>='+str(up_to), species[0]+str(i)] = (
            len(d['OGid'].unique()))
            OGid_interested_list = self.df_dups().query(
                f'{species[1]} == {i} & {species[0]} >= {up_to}'
            ).index.to_list()
            d = df.query(f'OGid in {OGid_interested_list}')
            ##print(d[['OGid', 'Species', 'ScaffoldCategory', 'OGtype']])##DEBUG
            df_count_ogtypes.loc[species[1]+str(i), species[0]+'>='+str(up_to)] = (
            len(d['OGid'].unique()))

        # compute both species >=up_to
        OGid_interested_list = self.df_dups().query(
            f'{species[1]} >= {up_to} & {species[0]} >= {up_to}'
        ).index.to_list()
        d = df.query(f'OGid in {OGid_interested_list}')
        ##print(d[['OGid', 'Species', 'ScaffoldCategory', 'OGtype']])##DEBUG
        df_count_ogtypes.loc[species[1]+'>='+str(up_to), species[0]+'>='+str(up_to)] = (
        len(d['OGid'].unique()))

        return df_count_ogtypes

    def heatmap(self, df,
                filename = 'heatmap_brocc_py.png'):
        """
        Create a heatmap from the pairwise connections between chr
        """
        plt.pcolormesh(df)
        plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
        plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns,
                   # rotate 45 degrees xticks
                   rotation=45)
        plt.savefig(filename)


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
    # Create MCL-like files
    ##broc.one_v_one_orthologs_dict(['Dcat', 'Dtil'])
    ##broc.paranome_OG_dict('Dcat')
    ##broc.paranome_OG_dict('Dtil')

    # let us filter the input dataframe by OGs present in both Dcat and Dtil:
    # (remove present only in one of these two species)
    OGid_interested_list = broc.df_dups().query(
        'Dcat > 0 & Dtil > 0').index.to_list()
    df_dcat_dtil = broc.df_input_table.query(f'OGid in {OGid_interested_list}')
    # remove other species, which we are not interested in...
    df_dcat_dtil = df_dcat_dtil.query('Species in ["Dcat", "Dtil"]')
    # create a new column, OGtype per species
    # instead of 1.2.3.4.5.6 etc. Dcat
    # only '1'
    df_dcat_dtil['OGtype_perspecies']=0
    for row in broc.df_dups().loc[:, ['Dcat', 'Dtil']].iterrows():
        og = row[0]
        row = row[1]
        for species in row.index:
            df_dcat_dtil.loc[(df_dcat_dtil['Species'] == species) &
                             (df_dcat_dtil['OGid'] == og),
                             'OGtype_perspecies'] = row[species]
    pd.set_option('display.min_rows', 10)
    pd.set_option('display.max_columns', 100)
    print(df_dcat_dtil.loc[:, ['OGtype', 'Species',
                               'OGtype_perspecies']].drop_duplicates())#DEBUG
    # write df to tsv
    df_dcat_dtil.to_csv('broquil_og_inboth_dcat_dtil.tsv', sep='\t',
        na_rep='NA') # don't remove missing values (GO_OGstatus)

    # filter df by OGids with all chr
    minor_dfs = []
    major_dfs = []
    mix_dfs = []
    unique_OGids = broc.df_input_table['OGid'].unique()
    for og in unique_OGids:
        filtered_df = broc.df_input_table.query(f'OGid == "{og}" &'+
                                    'Species in ["Dtil", "Dcat"]')
        if all(filtered_df['ScaffoldCategory']=='Chr'):
            major_dfs.append(filtered_df)
        elif all(filtered_df['ScaffoldCategory']=='Scaffold'):
            minor_dfs.append(filtered_df)
        else:
            mix_dfs.append(filtered_df)
    df_mix = pd.concat(mix_dfs).query('Method == "BROC"')
    df_major = pd.concat(major_dfs).query('Method == "BROC"')
    df_minor = pd.concat(minor_dfs).query('Method == "BROC"')
    df_total = broc.df_input_table.query('Method == "BROC"')
    print("Counts for the complete dataframe")
    df_count_total = broc.count_minor_major_scaffold_occupancy(['Dcat', 'Dtil'],
                                            up_to=6, df=df_total)
    print(df_count_total)
    print("Counts for major scaffolds only")
    df_count_major = broc.count_minor_major_scaffold_occupancy(['Dcat', 'Dtil'],
                                            up_to=6, df=df_major)
    print(df_count_major)
    print("Counts for minor scaffolds only")
    df_count_minor = broc.count_minor_major_scaffold_occupancy(['Dcat', 'Dtil'],
                                            up_to=6, df=df_minor)
    print(df_count_minor)
    print("Counts for OGs mixed between minor/major scaffolds")
    df_count_mix = broc.count_minor_major_scaffold_occupancy(['Dcat', 'Dtil'],
                                            up_to=6, df=df_mix)
    print(df_count_mix)
    # s'haurien de passar el dtype de cada cel·la a integer pq numpy pugui
    # interpretar-ho i fer-ne els heatmaps... df = df.astype(int)
    broc.heatmap(df_count_mix.astype(int), 'heatmap_brocc_mix.png')
    broc.heatmap(df_count_major.astype(int), 'heatmap_brocc_major.png')
    broc.heatmap(df_count_minor.astype(int), 'heatmap_brocc_minor.png')
    broc.heatmap(df_count_total.astype(int), 'heatmap_brocc_total.png')

###    # let us filter the input dataframe by OGs present in both Dcat and Dtil:
###    # (remove present only in one of these two species)
###    OGid_interested_list = broc.df_dups().query(
###        'Dcat > 0 & Dtil > 0').index.to_list()
###    df_dcat_dtil = broc.df_input_table.query(f'OGid in {OGid_interested_list}')
###    # remove other species, which we are not interested in...
###    df_dcat_dtil = df_dcat_dtil.query('Species in ["Dcat", "Dtil"]')
###    # create a new column, OGtype per species
###    # instead of 1.2.3.4.5.6 etc. Dcat
###    # only '1'
###    df_dcat_dtil['OGtype_perspecies']=0
###    for row in broc.df_dups().loc[:, ['Dcat', 'Dtil']].iterrows():
###        og = row[0]
###        row = row[1]
###        for species in row.index:
###            df_dcat_dtil.loc[(df_dcat_dtil['Species'] == species) &
###                             (df_dcat_dtil['OGid'] == og),
###                             'OGtype_perspecies'] = row[species]
###    pd.set_option('display.min_rows', 10)
###    pd.set_option('display.max_columns', 100)
###    print(df_dcat_dtil.loc[:, ['OGtype', 'Species',
###                               'OGtype_perspecies']].drop_duplicates())#DEBUG
###    print(df_dcat_dtil)
###    # write df to tsv
###    df_dcat_dtil.to_csv('broquil_og_inboth_dcat_dtil.tsv', sep='\t')


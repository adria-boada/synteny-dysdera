##! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# repeatmasker_parser4.py
#
# 10 Jun 2023  <adria@molevol-OptiPlex-9020>

"""

Let's parse from default output.
All changes ought to be done by pandas, temporally...
The input RM.out was modified as little as possible.

However, some changes are inevitable. With `awk`:
awk 'BEGIN{OFS="\t"} ; {if ( $16 ) { $16="True" } else { $16="False" } } 1' RM.out
There are rows with an asterisk as 16th field, and rows with only 15 fields.
Make sure all rows are of the same length (16 fields for all)

"""

import sys
import pandas as pd

# some csv are too big to read; they must be split beforehand
# (failed attempt)
import os, csv # reading and splitting csv files

# plotting
import seaborn as sns, matplotlib.pyplot as plt

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def print_green(text_string):
    """
    Print text to terminal with green colour
    """
    print('\033[0;32m'+
          text_string +
          '\033[0;0m')

    return None

def split_csv_file(input_file, output_prefix, lines_per_file):
    """
    Read a csv and split it in tmp files
    """
    with open(input_file, 'r') as file:
        reader = csv.reader(file)
        # pass the first 3 rows (unused multirow header)
        [next(reader) for i in range(3)]
        file_number = 1
        line_count = 0

        out_file_names = []
        while True:
            output_file = f"{output_prefix}_{file_number}.tmp"
            out_file_names.append(output_file)
            with open(output_file, 'w', newline='') as output:
                writer = csv.writer(output)

                for line in reader:
                    writer.writerow(line)
                    line_count += 1
                    if line_count >= lines_per_file:
                        break

                if line_count < lines_per_file:
                    break # end of read file

                line_count = 0
                file_number += 1

    return out_file_names

class Repeat:

    def __init__(self,
        file_list: "List of file-paths to a table.out given by RepeatMasker",
        species_names: "List of species names of each file-path, "+
            "in the same order (in order to pair file and sp.name)",
        gensizes_dict: "Dict with gensize value for species in species_names"):
        """
        """
        # try to open the file (make sure provided path is correct)
        for file in file_list:
            try:
                with open(file) as ft:
                    pass
            except:
                sys.exit('ERROR: The inputted file-path does not exist?')

        dfs = []
        for i in range(len(file_list)):
            # check filesize and split if it were too big
            # (always split to tmp file to shorten coding time)
            file = file_list[i]
            species = species_names[i]
            file_size = os.path.getsize(file)
            print_green("STATUS: Size of `" + str(file) + "` is " +
                        str(file_size) + " (bytes?)")
            tmp_file_list = split_csv_file(file, 'fragment_rmparser', 1000000)

            for tfl in tmp_file_list:
                # read and prepare dataframe from tabulated file
                df = pd.read_table(tfl,
                header=None,
                sep='\s+',
                names=[
                'score_SW',
                'perc_divg',
                'perc_del',
                'perc_ins',
                'sequid',
                'begin',
                'end',
                'left',
                'orient',
                'name',
                'default_repclass',
                'begin_match',
                'end_match',
                'left_match',
                'ID',
                'overlapping'],
                dtype={
                    'score_SW': int,        # Smith-Waterman score
                    'perc_divg': float,     # % of substit. bps in matching region
                    'perc_del': float,      # % of delet. bps in matching region
                    'perc_ins': float,      # % of insert. bps in matching region
                    'sequid': str,          # queried sequid name
                    'begin': int,           # match begin in sequid
                    'end': int,             # match end in sequid
                    'left': str,            # left to end of sequid
                    'orient': str,          # orientation/strand
                    'name': str,            # internal RM repeat name
                    'default_repclass': str,# default repeat class
                    'begin_match': str,     # begin in database match
                    'end_match': str,       # end in database match
                    'left_match': str,      # left to end of database match
                    'ID': int,              # identifier (not always unique)
                    'overlapping': bool},   # overlapping with a higher quality R.E.
                )
                # add more columns
                df['Species'] = species
                df['perc_indel'] = df['perc_del'] + df['perc_ins']
                df['replen'] = abs(df['end'] - df['begin']) +1
                dfs.append(df)

            # once read, remove all tmp files
            if any(".tmp" in s for s in tmp_file_list):
                print_green('STATUS: Removing temp. files from tmp_file_list')
                print(tmp_file_list)
                [os.remove(tfl) for tfl in tmp_file_list]

        # print('DEBUG: Preconcat')
        self.df_input_table = pd.concat(dfs)
        # print('DEBUG: Postconcat')

        # shorten calls to dataframe
        df = self.df_input_table
        # print(df) # DEBUG

        # Long paragraphs to classify "default_repclass" into
        # four new columns, easier to read and employ

        # UNKNOWN
        df.loc[
            df["default_repclass"]=='__unknown',
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Unclassified', 'NA', 'NA', 'NA')

        # "ENZYMATIC" RNAs
        df.loc[
            df["default_repclass"]=='tRNA',
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Other', 'tRNA', 'NA', 'NA')

        df.loc[
            df["default_repclass"]=='rRNA',
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Other', 'rRNA', 'NA', 'NA')

        df.loc[
            df["default_repclass"]=='srpRNA',
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Other', 'srpRNA', 'NA', 'NA')

        df.loc[
            df["default_repclass"]=='snRNA',
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Other', 'snRNA', 'NA', 'NA')

        df.loc[
            df["default_repclass"]=='scRNA',
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Other', 'scRNA', 'NA', 'NA')

        # GENERIC REPETITIVE ELEMENTS
        df.loc[
            df["default_repclass"]=='Simple_repeat',
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Tandem_repeat', 'Simple_repeat', 'NA', 'NA')

        df.loc[
            df["default_repclass"]=='Satellite',
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Tandem_repeat', '*Satellite*', 'NA', 'NA')

        df.loc[
            df["default_repclass"]=='Satellite/W-chromosome',
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Tandem_repeat', 'Satellite', '*W-chromosome*', 'NA')

        df.loc[
            df["default_repclass"]=='Satellite/Y-chromosome',
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Tandem_repeat', 'Satellite', '*Y-chromosome*', 'NA')

        df.loc[
            df["default_repclass"].str.contains('Satellite/acro'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Tandem_repeat', 'Satellite', '*acromeric*', 'NA')

        df.loc[
            df["default_repclass"].str.contains('Satellite/centr'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Tandem_repeat', 'Satellite', '*centromeric*', 'NA')

        df.loc[
            df["default_repclass"].str.contains('Satellite/macro'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Tandem_repeat', 'Satellite', '*macro*', 'NA')

        df.loc[
            df["default_repclass"].str.contains('Low_complexity'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Other', 'Low_complexity', 'NA', 'NA')

        df.loc[
            df["default_repclass"].str.contains('Other'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('NA', 'NA', 'NA', 'NA')


        # GENERIC RETROTRANSPOSONS
        df.loc[
            df["default_repclass"].str.contains('Retroposon') |
            df["default_repclass"].str.contains('__ClassI'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'NA', 'NA')

        df.loc[
            df["default_repclass"].str.contains('LINE'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'LINE', 'NA')

        df.loc[
            (df["default_repclass"].str.contains('_LTR')) |
            (df["default_repclass"] == 'LTR'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'LTR', 'NA')

        df.loc[
            df["default_repclass"].str.contains('SINE'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'SINE', 'NA')

        df.loc[
            df["default_repclass"].str.contains('Retroposon') &
            (df["default_repclass"].str.contains('SVA') |
            df["default_repclass"].str.contains('L1-dep')),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'NA', 'L1-dependent')

        df.loc[
            df["default_repclass"].str.contains('Retroposon/RTE-derived'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'NA', 'RTE-derived')

        # LINEs
        df.loc[
            (df["default_repclass"].str.contains('LINE') &
            df["default_repclass"].str.contains('Penelope')) |
            (df["default_repclass"].str.contains('ClassI') &
            df["default_repclass"].str.contains('PLE')),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'PLE', 'Penelope')

        df.loc[
            df["default_repclass"].str.contains('LINE/CR1') |
            df["default_repclass"].str.contains('LINE/L2') |
            df["default_repclass"].str.contains('LINE/Rex-Babar'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'LINE', 'CR1')

        df.loc[
            df["default_repclass"].str.contains('LINE/CRE'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'LINE', 'CRE')

        df.loc[
            df["default_repclass"].str.contains('LINE/Deceiver'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'LINE', 'Deceiver')

        df.loc[
            df["default_repclass"].str.contains('LINE/Dong-R4'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'LINE', 'Dong-R4')

        df.loc[
            df["default_repclass"].str.contains('LINE/I') |
            df["default_repclass"].str.contains('nLTR_LINE_I'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'LINE', 'I')

        df.loc[
            df["default_repclass"].str.contains('LINE/I-Jockey'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'LINE', 'Jockey')

        df.loc[
            df["default_repclass"].str.contains('LINE/L1') |
            df["default_repclass"].str.contains('nLTR_LINE_L1'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'LINE', 'L1')

        df.loc[
            df["default_repclass"].str.contains('LINE/Proto1'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'LINE', 'Proto1')

        df.loc[
            df["default_repclass"].str.contains('LINE/Proto2'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'LINE', 'Proto2')

        df.loc[
            df["default_repclass"].str.contains('LINE/R1') |
            df["default_repclass"].str.contains('LINE/Tad1'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'LINE', 'R1')

        df.loc[
            df["default_repclass"].str.contains('LINE/R2') |
            df["default_repclass"].str.contains('nLTR_LINE_R2'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'LINE', 'R2')

        df.loc[
            df["default_repclass"].str.contains('LINE/RTE') |
            df["default_repclass"].str.contains('nLTR_LINE_RTE'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'LINE', 'RTE')

        # LTRs
        df.loc[
            df["default_repclass"].str.contains('LTR') &
            df["default_repclass"].str.contains('DIRS'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'DIRS', 'DIRS')

        df.loc[
            df["default_repclass"].str.contains('LTR') &
            df["default_repclass"].str.contains('Ngaro'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'DIRS', 'Ngaro')

        df.loc[
            df["default_repclass"].str.contains('LTR/Viper'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'DIRS', 'Viper')

        df.loc[
            df["default_repclass"].str.contains('LTR/Pao') |
            df["default_repclass"].str.contains('LTR_BEL'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'LTR', 'Bel-Pao')

        df.loc[
            df["default_repclass"].str.contains('LTR') &
            df["default_repclass"].str.contains('Copia'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'LTR', 'Copia')

        df.loc[
            df["default_repclass"].str.contains('LTR') &
            df["default_repclass"].str.contains('ERV'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'LTR', 'ERV')

        df.loc[
            df["default_repclass"].str.contains('LTR') &
            df["default_repclass"].str.contains('Gypsy'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'LTR', 'Gypsy')

        df.loc[
            df["default_repclass"].str.contains('LTR') &
            df["default_repclass"].str.contains('ERV'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'LTR', 'ERV')

        df.loc[
            df["default_repclass"].str.contains('LTR') &
            df["default_repclass"].str.contains('Caulimovirus'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'LTR', 'Pararetrovirus')

        # SINEs
        df.loc[
            df["default_repclass"].str.contains('SINE') &
            df["default_repclass"].str.contains('5S'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'SINE', '5S')

        df.loc[
            df["default_repclass"].str.contains('SINE') &
            (df["default_repclass"].str.contains('7SL') |
            df["default_repclass"].str.contains('ALU', case=False) |
            df["default_repclass"].str.contains('B2') |
            df["default_repclass"].str.contains('B4')),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'SINE', '7SL')

        df.loc[
            df["default_repclass"].str.contains('Retroposon/L1-derived'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'SINE', 'L1-derived')

        df.loc[
            df["default_repclass"].str.contains('Retroposon/L2-derived'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'SINE', 'L2-derived')

        df.loc[
            df["default_repclass"].str.contains('SINE') &
            df["default_repclass"].str.contains('tRNA'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'SINE', 'tRNA')

        df.loc[
            df["default_repclass"].str.contains('Retroposon/R4-derived'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Retrotransposon', 'NA', 'SINE', 'R4-derived')

        # GENERIC DNAs
        df.loc[
            (df["default_repclass"] == 'DNA') |
            (df["default_repclass"] == '__ClassII'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('DNA', 'NA', 'NA', 'NA')

        # DNA I
        df.loc[
            df["default_repclass"].str.contains('DNA') &
            df["default_repclass"].str.contains('Crypton'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('DNA', '1', 'Crypton', 'Crypton')

        df.loc[
            df["default_repclass"].str.contains('DNA') &
            df["default_repclass"].str.contains('Zator'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('DNA', '1', 'NA', 'Zator')

        df.loc[
            df["default_repclass"].str.contains('ClassII_MITE'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('DNA', '1', 'NA', 'MITE')

        df.loc[
            df["default_repclass"].str.contains('ClassII_nMITE'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('DNA', '1', 'NA', 'nMITE')

        df.loc[
            df["default_repclass"].str.contains('DNA') &
            df["default_repclass"].str.contains('Academ'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('DNA', '1', 'TIR', 'Academ')

        df.loc[
            (df["default_repclass"].str.contains('DNA') &
            df["default_repclass"].str.contains('CMC')) |
            df["default_repclass"].str.contains('ClassII_DNA_CACTA'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('DNA', '1', 'TIR', 'CACTA')

        df.loc[
            (df["default_repclass"].str.contains('DNA') &
            df["default_repclass"].str.contains('Dada')),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('DNA', '1', 'TIR', 'Dada')

        df.loc[
            (df["default_repclass"].str.contains('DNA') &
            df["default_repclass"].str.contains('Kolobok')),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('DNA', '1', 'TIR', 'Kolobok')

        df.loc[
            (df["default_repclass"].str.contains('DNA') &
            df["default_repclass"].str.contains('Ginger')),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('DNA', '1', 'TIR', 'Ginger')

        df.loc[
            (df["default_repclass"].str.contains('DNA') &
            df["default_repclass"].str.contains('MULE')),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('DNA', '1', 'TIR', 'MULE')

        df.loc[
            (df["default_repclass"].str.contains('DNA') &
            df["default_repclass"].str.contains('Merlin')),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('DNA', '1', 'TIR', 'Merlin')

        df.loc[
            (df["default_repclass"].str.contains('DNA/P')),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('DNA', '1', 'TIR', 'P')

        df.loc[
            (df["default_repclass"].str.contains('DNA/PIF')) |
            (df["default_repclass"].str.contains('ClassII_DNA_Harbinger')),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('DNA', '1', 'TIR', 'PIF-Harbinger')

        df.loc[
            (df["default_repclass"].str.contains('PiggyBac')),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('DNA', '1', 'TIR', 'PiggyBac')

        df.loc[
            (df["default_repclass"].str.contains('DNA/Sola')),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('DNA', '1', 'TIR', 'Sola')

        df.loc[
            (df["default_repclass"].str.contains('DNA/Sola')),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('DNA', '1', 'TIR', 'Sola')

        df.loc[
            (df["default_repclass"].str.contains('DNA/TcMar') |
            df["default_repclass"].str.contains('DNA_TcMar')),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('DNA', '1', 'TIR', 'TcMar')

        df.loc[
            (df["default_repclass"].str.contains('DNA') &
            df["default_repclass"].str.contains('hAT')),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('DNA', '1', 'TIR', 'hAT')

        df.loc[
            (df["default_repclass"].str.contains('DNA_Mutator')),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('DNA', '1', 'TIR', 'Mutator')

        df.loc[
            df["default_repclass"].str.contains('Maverick'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('DNA', '1', 'Maverick', 'Maverick')

        df.loc[
            df["default_repclass"].str.contains('IS3EU'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('DNA', '1', 'NA', 'IS3EU')

        df.loc[
            df["default_repclass"].str.contains('Novosib'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('DNA', '1', 'NA', 'Novosib')

        df.loc[
            df["default_repclass"].str.contains('Zisupton'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('DNA', '1', 'NA', 'Zisupton')

        # DNA II
        df.loc[
            (df["default_repclass"].str.contains('Helitron')),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('DNA', '2', 'Helitron', 'Helitron')

        # MITE or nMITE?
        df.loc[
            df["default_repclass"].str.contains('_MITE'),
            ['mite']
        ] = (True)

        df.loc[
            df["default_repclass"].str.contains('_nMITE'),
            ['mite']
        ] = (False)

        # check dtypes have recorded correctly
        print_green(f'STATUS: Dataframe dtypes:')
        self.df_input_table.info()

        # create a complementary dataframe which summarizes `df_input_table`
        # computes sum of replens (cumulative length of repeat class)
        # mean of replen (average length of repeat class)
        # count of replen (number of elements of repeat class)
        g = df.groupby(['Species', 'class', 'subclass', 'order',
                        'superfam'])['replen'].agg([
                        'sum', 'mean', 'count',]).reset_index()
        g = g.rename(columns={
            'sum': 'cumulative_bp_len',
            'mean': 'avgele_bp_len',
            'count': 'num_elements'})

        # compute proportion of genome occupancy for all classes,
        # for a given species
        for sp in gensizes_dict.keys():
            g.loc[g['Species'] == sp, 'genome_occupancy'] = (
                g.loc[g['Species'] == sp, 'cumulative_bp_len'] /
                int(gensizes_dict[sp])) * 100 # percent

        # compute proportion of repetitive fraction for all classes,
        # given a species
        for sp in g['Species'].unique():
            g.loc[g['Species'] == sp, 'relative_repeat_fraction'] = (
                g.loc[g['Species'] == sp, 'cumulative_bp_len'] /
                g.loc[g['Species'] == sp, 'cumulative_bp_len'].sum()) * 100

        # summarize `df_groupby_summary` (bin superfamilies into classes)
        # reduces information even further than `g`/`df_groupby_summary`
        df_compressed_summary = df.groupby(['Species', 'class'])['replen'
                                              ].agg(['sum','mean','count']
                                                  ).reset_index()
        total = df.groupby(['Species'])['replen'
                                              ].agg(['sum','mean','count']
                                                  ).reset_index()
        total['class'] = 'Total'
        df_compressed_summary = pd.concat([df_compressed_summary, total])
        df_compressed_summary = df_compressed_summary.rename(columns={
            'sum': 'cumulative_bp_len',
            'mean': 'avgele_bp_len',
            'count': 'num_elements'})

        # Compute proportions for the last df, too...
        for sp in gensizes_dict.keys():
            df_compressed_summary.loc[
                df_compressed_summary['Species'] == sp, 'genome_occupancy'] = (
                df_compressed_summary.loc[
                    df_compressed_summary['Species'] == sp, 'cumulative_bp_len'] /
                int(gensizes_dict[sp])
                ) * 100 # percent
        for sp in df_compressed_summary['Species'].unique():
            # total value of repeated bps:
            tot_bp = int(df_compressed_summary.loc[
                (df_compressed_summary['Species'] == sp) &
                (df_compressed_summary['class'] == 'Total'),
                'cumulative_bp_len'])
            df_compressed_summary.loc[
                df_compressed_summary['Species'] == sp, 'relative_repeat_fraction'] = (
                    df_compressed_summary.loc[
                        df_compressed_summary['Species'] == sp, 'cumulative_bp_len'] /
                    tot_bp) * 100

        # hauria de fer una funció per agrupar, però qui té temps per funcionar?
        # Broaden `df_compressed_summary` scope to introduce content per chr.

        # subset all repeats in chrs (and minor sequids)
        df = self.df_input_table.loc[:]
        # rename all scaffolds, so they share a single name; simply, "Scaffold"
        # (instead of scaffold_1, scaffold_2, etc.)
        df.loc[df['sequid'].str.contains('Scaffold'), 'sequid'] = "Scaffold"
        # groupby species, sequid and class...
        df_summ_per_seq = self.df_input_table.groupby(['Species', 'sequid', 'class'])['replen'
                                              ].agg(['sum','mean','count']
                                                  ).reset_index()
        total = df.groupby(['Species', 'sequid'])['replen'
                                              ].agg(['sum','mean','count']
                                                  ).reset_index()
        total['class'] = 'Total'
        df_summ_per_seq = pd.concat([df_summ_per_seq, total])
        df_summ_per_seq = df_summ_per_seq.rename(columns={
            'sum': 'cumulative_bp_len',
            'mean': 'avgele_bp_len',
            'count': 'num_elements'})

        # Compute proportions for the last df, too...
        for sp in gensizes_dict.keys():
            msp = df_summ_per_seq['Species']==sp
            for sequid in df_summ_per_seq.loc[msp,
                                              'sequid'].unique():
                msq = df_summ_per_seq['sequid']==sequid
                df_summ_per_seq.loc[(msp) & (msq), 'genome_occupancy'] = (
                    df_summ_per_seq.loc[(msp) & (msq), 'cumulative_bp_len'] /
                int(gensizes_dict[sp])
                ) * 100 # percent
        for sp in gensizes_dict.keys():
            msp = df_summ_per_seq['Species']==sp
            for sequid in df_summ_per_seq.loc[msp,
                                              'sequid'].unique():
                msq = df_summ_per_seq['sequid']==sequid
                tot_bp = int(df_summ_per_seq.loc[
                    (msq) & (msp) &
                    (df_summ_per_seq['class'] == 'Total'),
                    'cumulative_bp_len'])
                df_summ_per_seq.loc[
                    (msq) & (msp),
                    'relative_repeat_fraction'] = (
                df_summ_per_seq.loc[
                    (msq) & (msp),
                    'cumulative_bp_len'] / tot_bp) * 100

        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_colwidth', None)
        print_green("STATUS: Summary of repeatmasker's dataframe.out:")
        # print df in markdown format (requires tabulate package)
        print(g.round(decimals=2).to_markdown(None))
        self.df_groupby_summary = g
        self.df_compressed_summary = df_compressed_summary.reset_index(drop=True)
        self.df_summ_per_seq = df_summ_per_seq.reset_index(drop=True)

    # this is a long __init__ function...

    def check_classification(self):
        """
        Checks how the default classes have been re-classified
        """
        df = self.df_input_table # check on created self.dataframe

        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_colwidth', None)
        print_green(f'STATUS: Reclassification of dataframe:')
        # show the new classification
        # (drop rows sharing the same repeat class)
        df_new_classes = df.drop_duplicates(
            subset=['default_repclass', 'class',
                    'subclass', 'order', 'superfam'])
        # sort dataframe by classes (then subclass, etc.)
        df_new_classes = df_new_classes.sort_values(
            by=['class', 'subclass', 'order', 'superfam'])
        # reset index
        df_new_classes = df_new_classes.reset_index()
        # only print columns with repeat class names
        df_new_classes = df_new_classes.loc[:, ['class',
                                     'subclass', 'order',
                                     'superfam', 'mite',
                                     'default_repclass']]
        # print df in markdown format (requires tabulate package)
        print(df_new_classes.round(decimals=2).to_markdown(None))
        # return to default pandas options
        pd.set_option('display.max_rows', 60)
        pd.set_option('display.max_colwidth', 80)

        return df_new_classes

    def boxplot_entrance(self):
        """
        Draw boxplots of dataframe by
            + X variable (e.g. SW score, genome content)
            + Y category (e.g. repeat class)
            + Z hue (e.g. species)
        """
        df = self.df_input_table
        fig_titles = {
            'score_SW': 'Smith-Waterman score',
            'perc_divg': 'Divergence from consensus sequence',
            'perc_indel': 'Indels from consensus sequence',
            'replen': 'Repetitive element lengths',
        }
        fig_xlabel = {
            'score_SW': 'score',
            'perc_divg': '% divergence',
            'perc_indel': '% indels (?)',
            'replen': 'basepairs',
        }
        for xvar in ['score_SW', 'perc_divg', 'perc_indel',
                     'replen']:
            sns.boxplot(data=df, x=xvar, y='class',
                        hue='Species', showmeans='True',
                        meanprops=dict(color='green',
                                       marker='*',
                                       markeredgecolor='black',
                                       markersize='15'))
            plt.title(fig_titles[xvar])
            plt.subplots_adjust(left=0.25, bottom=0.1, right=0.95, top=0.91)
            plt.xlabel(fig_xlabel[xvar])
            plt.savefig(str(xvar)+'_matplotlib.png', dpi=300)
            plt.close('all')

        return None

    def boxplot_genome_content(self):
        """
        Draw boxplots of genome content occupied
        by each class in each species.
        """
        g  = self.df_groupby_summary
        fig_titles = {
            'avgele_bp_len': 'Average repeat class length',
            'num_elements': 'Number of elements of repeat class',
            'genome_occupancy': 'Genome occupancy for repeat class'
        }
        fig_xlabel = {
            'avgele_bp_len': 'basepairs',
            'num_elements': '#amount eles',
            'genome_occupancy': '% occupancy'
        }
        for xvar in ['avgele_bp_len', 'num_elements', 'genome_occupancy']:
            sns.boxplot(data=g, x=xvar, y='class',
                        hue='Species', showmeans='True',
                        meanprops=dict(color='green',
                                       marker='*',
                                       markeredgecolor='black',
                                       markersize='15'))
            # lets try to label outlier subclasses in a class
            xvar_q1 = g.groupby(['class', 'Species']
                            )[xvar].quantile(0.25)
            xvar_q3 = g.groupby(['class', 'Species']
                            )[xvar].quantile(0.75)
            outlier_top_lim = xvar_q3 + 1 * (xvar_q3 - xvar_q1)
            outlier_bot_lim = xvar_q1 - 1 * (xvar_q3 - xvar_q1)

            categorical_pos = {
                'Dcat*Low_complexity*': 0,
                'Dtil*Low_complexity*': 0,
                'DcatDNA': 1,
                'DtilDNA': 1,
                'DcatRetrotransposon': 2,
                'DtilRetrotransposon': 2,
                'DcatTandem_repeat': 3,
                'DtilTandem_repeat': 3,
                'DcatUnclassified': 4,
                'DtilUnclassified': 4,
                'DcatrRNA': 5,
                'DtilrRNA': 5,
                'DcattRNA': 6,
                'DtiltRNA': 6,
            }
            # labeling values not bonded by limits
            for row in g.to_dict(orient='index').items():
                index = row[0]
                rclass = row[1]['class']
                sp = row[1]['Species']
                pos = categorical_pos[sp+rclass]
                val = row[1][xvar]
                if (val > outlier_top_lim[rclass, sp]) or (
                    val < outlier_bot_lim[rclass, sp]):
                    plt.text(val, pos, f' {index}',
                             ha='left', va='center', fontsize=4)

            plt.title(fig_titles[xvar])
            plt.subplots_adjust(left=0.25, bottom=0.1, right=0.95, top=0.91)
            plt.xlabel(fig_xlabel[xvar])
            plt.savefig(str(xvar)+'_matplotlib.png', dpi=300)
            plt.close('all')

        return None


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Si es crida com a script:
if __name__ == '__main__':

    ### PARSE ARGUMENTS ###
    if len(sys.argv) < 3:
        sys.exit(f"Lacking RepeatMasker's output\n{sys.argv[0]} rm_sp1.out "+
                 "species1_name rm_sp2.out species2_name ...")

    species = []
    files = []
    for i in range(1, len(sys.argv), 2):
        files.append(sys.argv[i])
        species.append(sys.argv[i+1])

    # send genome length in bases (compute % genome occupancy)
    repeats = Repeat(files, species,
                     gensizes_dict=dict(
                         Dcat=3279770110,
                         Dtil=1550000000))
    # check wheather classification is done correctly
    # store reclassification in a TSV file
    repeats.check_classification(
        ).round(decimals=2).to_csv('repeat_reclassification.tsv',
        sep='\t', na_rep='NA',)
    # moreover, store summary df as TSV
    repeats.df_groupby_summary.round(decimals=2).to_csv(
        'repeat_table_superfamilies.tsv', sep='\t', na_rep='NA')
    # store dfs that further summarize information;
    # more concisely:
    repeats.df_compressed_summary.round(decimals=2).to_csv(
        'repeat_table_classes.tsv', sep='\t', na_rep='NA')
    # per sequid/chromosomes:
    repeats.df_summ_per_seq.round(decimals=2).to_csv(
        'repeat_table_per_sequid.tsv', sep='\t', na_rep='NA')
    # create boxplots of whole df and summary df
    repeats.boxplot_entrance()        # whole
    repeats.boxplot_genome_content()  # summary


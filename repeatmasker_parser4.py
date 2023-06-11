#! /usr/bin/env python3
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
import os # reading
import csv # splitting

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

    def __init__(self, file:"File-path to a table.out given by RepeatMasker"):
        """
        """
        # try to open the file (make sure provided path is correct)
        try:
            with open(file) as ft:
                pass
        except:
            sys.exit('ERROR: The inputted file-path does not exist?')

        # check filesize and split if it were too big
        file_size = os.path.getsize(file)
        print_green("STATUS: CSV file-size of " + str(file_size) + " (bytes?)")
        if file_size > 500000000: # file-size bigger than 500 Mb?
            # split before loading
            print_green(f'STATUS: File-size exceeding 5 Mb, memory compromised. '+
                        'File will be split in tmp files and read '+
                        'individually.')

            tmp_file_list = split_csv_file(file, 'fragment_rmparser', 1000000)

        else:
            tmp_file_list = [file]

        print('DEBUG: split file')

        dfs = []
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
            dfs.append(df)

        print('DEBUG: Preconcat')
        self.df_input_table = pd.concat(dfs)
        print('DEBUG: Postconcat')
        # once read, remove all tmp files
        print_green('STATUS: Removing files from tmp_file_list')
        print(tmp_file_list)
        [os.remove(tfl) for tfl in tmp_file_list]

        print_green(f'STATUS: The file {file} '+
                    'has been read as the following dataframe:')
        print(self.df_input_table)

        # shorten calls to dataframe
        df = self.df_input_table

        # Long paragraphs to classify "default_repclass" into
        # four new columns, easier to read and employ

        # UNKNOWN
        df.loc[
            df["default_repclass"]=='__unknown',
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Unknown', 'NA', 'NA', 'NA')

        # "ENZYMATIC" RNAs
        df.loc[
            df["default_repclass"]=='tRNA',
            ['class', 'subclass', 'order', 'superfam']
        ] = ('tRNA', 'NA', 'NA', 'NA')

        df.loc[
            df["default_repclass"]=='rRNA',
            ['class', 'subclass', 'order', 'superfam']
        ] = ('rRNA', 'NA', 'NA', 'NA')

        df.loc[
            df["default_repclass"]=='srpRNA',
            ['class', 'subclass', 'order', 'superfam']
        ] = ('srpRNA', 'NA', 'NA', 'NA')

        df.loc[
            df["default_repclass"]=='snRNA',
            ['class', 'subclass', 'order', 'superfam']
        ] = ('snRNA', 'NA', 'NA', 'NA')

        df.loc[
            df["default_repclass"]=='scRNA',
            ['class', 'subclass', 'order', 'superfam']
        ] = ('scRNA', 'NA', 'NA', 'NA')

        # GENERIC REPETITIVE ELEMENTS
        df.loc[
            df["default_repclass"]=='Simple_repeat',
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Tandem_repeat', 'Simple_repeat', 'NA', 'NA')

        df.loc[
            df["default_repclass"]=='Satellite',
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Tandem_repeat', 'Satellite', 'NA', 'NA')

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
        ] = ('Low_complexity', 'Low_complexity', 'NA', 'NA')

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

        # check transformation
        pd.set_option('display.max_rows', None)
        print_green(f'STATUS: New columns added:')
        print(df.drop_duplicates(subset=['default_repclass',
                                         'class', 'subclass',
                                         'order',
                                         'superfam']).loc[:, [
                                     'class', 'subclass', 'order',
                                     'superfam', 'mite',
                                     'default_repclass']].sort_values(
                                        by=['class', 'subclass', 'order',
                                            'superfam']))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Si es crida com a script:
if __name__ == '__main__':

    ### PARSE ARGUMENTS ###
    import argparse

    parser = argparse.ArgumentParser(description='')
    # file-name: positional arg.
    parser.add_argument('filename', type=str, help='Path to ... file-name')
    # integer argument
#    parser.add_argument('-a', '--numero_a', type=int, help='Paràmetre "a"')
    # choices argument
#    parser.add_argument('-o', '--operacio', 
#            type=str, choices=['suma', 'resta', 'multiplicacio'],
#            default='suma', required=False,
#            help='')

    args = parser.parse_args()
    # call a value: args.operacio or args.filename.

    repeats = Repeat(args.filename)

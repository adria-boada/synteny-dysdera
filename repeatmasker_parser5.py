#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# repeatmasker_parser5.py
#
# 22 Jun 2023  <adria@molevol-OptiPlex-9020>

"""
"""

import sys

# data (frames) handling
import pandas as pd
import os # compute file-size

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

"""
Funció per transformar raw_file.out a modified_file.out

Parse 16th row of 'overlapping' repeats; it can either be an asterisk or a
blank field ("*" or ""). Blank fields are not loaded by pandas.read_table().
Turn asterisk into "True" and blank into "False" so every row has a 16th field.

with open(inputfile) as file:
  [file.readline() for i in range(3)]
  for line in file:
    l = line.split()
    x = l[len(l) - 1]
    if x=="*":
      sixtenth = "True"
    else:
      sixtenth = "False"
      ...
"""

def remove_overlapping(intervals):
    """
    Parse intervals of a given Species-sequid and remove overlapping.

    INPUT

    * intervals: A list with sublists. Each sublist is a couple of coordinates;
    begin and end (integers). Add score for intervals as the third item of
    sublists. Include a tuple with (class, subclass, order, superfam)

    [[begin_1, end_1, score_1, class_tuple_1], ...,
    [begin_N, end_N, score_N, class_tuple_N]]

    OUTPUT

    Sum of lengths, each one computed as `(end-begin)+1`. Removes all
    overlapping segments (counts each nucleotide a single time).
    """
    # Init dataframe which will be returned
    # compute unique repeat types
    unique_class_tuples = list(set(tuple(x[3]) for x in intervals))
    unique_class_tuples.sort()
    # covert into pd.MultiIndex()
    idx = pd.MultiIndex.from_tuples(unique_class_tuples, names=[
            'class', 'subclass', 'order', 'superfam'])
    # create df
    df_resulting = pd.DataFrame(index=idx)
    df_resulting['algor_bp'] = 0

    # Create a point-list from an interval-list
    # from [begin, end], [begin, end], etc.
    # to [coord, begin, 0], [coord, end, 0], [coord, begin, 1], etc.
    point_list = []
    for i in range(0, len(intervals)):
        # Left/starting points labeled '0'
        point_list += [ [intervals[i][0], 0, i] ]
        # Right/ending points labeled '1'
        point_list += [ [intervals[i][1], 1, i] ]
    # sort point list by position (by first value in sublists)
    point_list.sort()

    # Init algorithm variables
    currentBest = -1
    open_intervals = []

    # for each point in point list:
    for i in range(0, len(point_list)):
        # DEBUG #
##        print(df_resulting)
##        print('currentBest:', currentBest)
##        print('iterate:', point_list[i])
##        print('open_intervals:', open_intervals)
        # DEBUG #

        # If the loop landed on a left point (0)
        # opens an interval
        if point_list[i][1] == 0:

            # if there is no other interval opened:
            if currentBest == -1:
                # enters interval 'i'
                currentBest = point_list[i][2]
                currentBegin = int(point_list[i][0])
                currentScore = int(intervals[currentBest][2])

            # else, there already was an open interval:
            # (and it is of lower score than currentBest)
            elif currentScore > intervals[point_list[i][2]][2]:
                # do not close currentBest, because the next interval (repeat)
                # is of lower score; however, append to open intervals list
                iden = point_list[i][2]
                open_intervals.append(
                    [int(intervals[iden][2]),     # score
                     iden])                       # ID
                open_intervals.sort(reverse=True) # sort by score

            # else, currentScore should be higher
            else:
                # compute length up until now (begin2 -begin1)
                interval_length = point_list[i][0] -currentBegin
                # add this length to the resulting pandas df
                reptype = intervals[currentBest][3]
                df_resulting.loc[reptype, "algor_bp"] += interval_length
                # append index of past Best to open_intervals
                iden = currentBest
                open_intervals.append(
                    [int(intervals[iden][2]),     # score
                     iden])                       # ID
                open_intervals.sort(reverse=True) # sort by score
                # and change currentBest
                currentBest = point_list[i][2]
                currentBegin = point_list[i][0]
                currentScore = intervals[currentBest][2]

        # If the loop landed on a right point (1)
        # closes an interval
        # moreover, it has to be the right point of the
        # currently open interval 'i'
        elif point_list[i][2] == currentBest:
            # DEBUG
##            print('point_list[i][2] (', point_list[i], ') == currentBest')
            # compute length up until now (end -begin +1)
            interval_length = point_list[i][0] -currentBegin +1
            # add this length to the resulting pandas df
            reptype = intervals[currentBest][3]
            df_resulting.loc[reptype, "algor_bp"] += interval_length
            # Close this interval and open the next best one in open_intervals
            if len(open_intervals) == 0:
                currentBest = -1
            else:
                # remove second best from open_intervals and use it as
                # currentBest
                c = open_intervals.pop(0)
                currentBest = c[1]
                currentScore = c[0]
                currentBegin = point_list[i][0] +1

        # Otherwise, it is closing an interval listed in `open_intervals`
        else:
            # remove the interval from `open_intervals`
            iden = point_list[i][2]
            open_intervals.remove(
                [int(intervals[iden][2]), # score
                iden])                    # ID

    return df_resulting

class Repeat:

    def __init__(self,
                 file_list: "List of file-paths to `modified_table.out` given "+
                 "by Repeatmasker",
                 seqsizes_dict: "Dictionary with species included, their "+
                 "original sequid name, their refined sequid name, sequid "+
                 "lengths, etc. See example below.",
                 ):
        """
        INPUT:

        file_list: A list of 'refined' or 'modified' files. Can be of any length.
        Should be congruent with `seqsizes_dict` species order.
        ["modified_table_sp1.out", modified_table_sp2.out", etc.]

        seqsizes_dict: A very informative dictionary.
        `chr1`, `chr2`, etc. should be a regex which matches the sequid(s) of
        interest (and nothing else!).
        Use '|' for OR matching: `string1|string2|string3`, `Scaff|ctg`...
        Use '&' for AND matching: `Sca&ffolds`...

        dict('species1': {'chr1': length.int,
                          'chr2': length.int,
                            ...
                          'minor_scaffolds': sum_lengths.int
                         },
             'species2': {...},
             ...)

        OUTPUT:

        Dataframes ready to analyze.
        """

        # might be useful further on
        self.seqsizes_dict = seqsizes_dict

        # complete genome sizes
        gensizes_dict = dict()
        for species in seqsizes_dict.keys():
            gensizes_dict[species] = 0
            for sequid in seqsizes_dict[species].keys():
                gensizes_dict[species] += seqsizes_dict[species][sequid]
        self.gensizes_dict = gensizes_dict

        # make sure path is correct
        for file in file_list:
            try:
                with open(file) as ft:
                    pass
            except:
                sys.exit(f'ERROR: Cannot read the provided path: {file}')

        # read and load input TSV files as dataframes
        dfs = []
        print_green("STATUS: FILESIZES")
        for i in range(len(file_list)):
            # check file-size
            file = file_list[i]
            species = list(seqsizes_dict.keys())[i]
            file_size = os.path.getsize(file)
            print("File `" + str(file) + "` is " +
                        str(file_size) + " (bytes?)")

            # read tabulated file into df with pandas
            df = pd.read_table(file,
            header=None,
            skiprows=3,    # skip first 3 header rows
            sep='\s+',     # fields separated by '\s+' regex
            names=[        # use these column names...
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

        print()
        # concat all given files from different species
        df_concat = pd.concat(dfs)

        # count rows matching a sequid in seqsizes_dict
        # and rows not matching a single sequid
        # also, display what is the key matching (print unique matching fields)
        print_green(
            "STATUS: REVISING REGEX MATCHES BETWEEN PROVIDED DICT AND TSV")
        for species in seqsizes_dict.keys():
            df_species = df_concat.loc[df_concat["Species"] == species]
            df_species_totnrows = df_species.shape[0]
            df_species_matchrows = 0
            for sequid_regex in seqsizes_dict[species].keys():
                df_species_and_regex = df_species.loc[
                    df_species["sequid"].str.contains(sequid_regex)]
                print("--", sequid_regex, "--")
                print(" + Nrows:", df_species_and_regex.shape[0])
                df_species_matchrows += df_species_and_regex.shape[0]
                print(" + Regex matching count:",
                      len(list(df_species_and_regex["sequid"].unique())),
                      "unique sequids")
                print(" + Regex matching list:",
                      list(df_species_and_regex["sequid"].unique()))
                # add new col to keep track of these 'sequid_types'
                df_concat.loc[
                    (df_concat["Species"] == species) &
                    (df_concat["sequid"].str.contains(sequid_regex)),
                    "sequid_type"] = sequid_regex
            print("-----------------------------")
            print(species, "summary: matched", df_species_matchrows,
                  "out of", df_species_totnrows, "rows (",
                  round((df_species_matchrows/df_species_totnrows) *100, ndigits=3),
                  "%)\n")

        # Long paragraphs to classify "default_repclass" into
        # four new columns, easier to read and employ
        df = df_concat

        # UNKNOWN
        df.loc[
            df["default_repclass"]=='__unknown',
            ['class', 'subclass', 'order', 'superfam']
        ] = ('Unknown', 'NA', 'NA', 'NA')

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

        # check whether dtypes have been recorded properly
        print_green(
            "STATUS: INFO ABOUT LOADED DATAFRAME")
        df_concat.info()
        df_concat.head()

        # Summarize dataframe into absolute bp count per pairs of Species,
        # sequid-regex for each type of repeat.
        df_naive = df.groupby(
            ["Species", "sequid_type", "class", "subclass", "order",
                    "superfam"])["replen"].agg(["count",   # sum lens of repeats
                                                "mean",  # avg len of repeat
                                                "sum"] # amount of repeat
                                        )
        # rename columns
        df_naive = df_naive.rename(columns={
            'sum': 'naive_bpsum',
            'mean': 'naive_avglen',
            'count': 'naive_numele'})

        # remove 'Overlapping' tagged repeats by RepeatMasker from the sum
        df = df_concat.loc[df_concat["overlapping"] == False]
        df_tagged = df.groupby(
            ["Species", "sequid_type", "class", "subclass", "order",
                    "superfam"])["replen"].agg(["count",   # sum lens of repeats
                                                "mean",  # avg len of repeat
                                                "sum"] # amount of repeat
                                        )
        # rename columns
        df_tagged = df_tagged.rename(columns={
            'sum': 'tagged_bpsum',
            'mean': 'tagged_avglen',
            'count': 'tagged_numele'})

        # merge naive, tagged dataframes
        df_absolute_summary = df_naive.join([df_tagged]).reset_index()

        # Remove overlapping regions algorithmically
        # the function `remove_overlapping` outputs a df with a "algor_bp"
        # column
        df_absolute_summary["algor_bp"] = 0

        pd.options.display.max_rows = None # debug
        print(df_absolute_summary.head())  # debug

        for species in df_concat["Species"].unique():
            sequid_list = df_concat.loc[
                df_concat["Species"] == species,
                "sequid"].unique()
            for seq in sequid_list:
                d = df_concat.loc[
                    (df_concat["Species"]==species) &
                    (df_concat["sequid"]==seq), ]
                seqtype = d["sequid_type"].unique()[0]
                tuple_repclass = tuple(zip(
                    list(d["class"]),
                    list(d["subclass"]),
                    list(d["order"]),
                    list(d["superfam"]) ))
                intervals = list(zip(
                    list(d["begin"]),
                    list(d["end"]),
                    list(d["score_SW"]),
                    tuple_repclass
                ))
                df_return = remove_overlapping(intervals).reset_index()
                df_return.insert(0, "Species", [species]*df_return.shape[0])
                df_return.insert(1, "sequid_type", [seqtype]*df_return.shape[0])
                df_absolute_summary = pd.merge(right=df_absolute_summary, left=df_return,
                                               # preserve right dataframe's keys
                                               how="right",
                                               on=["Species", "sequid_type",
                                                   "class", "subclass", "order",
                                                   "superfam"])
                print(df_absolute_summary[["algor_bp_x", "algor_bp_y"]].info()) # debug
                df_absolute_summary["algor_bp_x"] = (
                df_absolute_summary["algor_bp_x"].fillna(0))
                df_absolute_summary["algor_bp_y"] = (
                    df_absolute_summary["algor_bp_x"] +
                    df_absolute_summary["algor_bp_y"])
                df_absolute_summary.drop(columns=["algor_bp_x"], inplace=True)
                df_absolute_summary.rename(columns={"algor_bp_y": "algor_bp"},
                                           inplace=True)
                print(df_absolute_summary[["algor_bp"]].info()) # debug

        # compute non-repeat size per sequid_type
        for species in df_absolute_summary["Species"].unique():
            for sequid in df_absolute_summary.loc[
                    df_absolute_summary["Species"] == species,
                    "sequid_type"].unique():
                # compute fraction of repetitive repeats
                repetitive_bps = {
                    'naive': int(df_absolute_summary.loc[
                    (df_absolute_summary["Species"] == species) &
                    (df_absolute_summary["sequid_type"] == sequid),
                    "naive_bpsum"].sum()),
                    'tagged': int(df_absolute_summary.loc[
                    (df_absolute_summary["Species"] == species) &
                    (df_absolute_summary["sequid_type"] == sequid),
                    "tagged_bpsum"].sum()),
                }
                # compute non-repeat fraction
                # add a row with non-repeat and repeat total bps
                new_rows = {
                    "Species": [species] *2,
                    "sequid_type": [sequid] *2,
                    "class": ["Nonrepetitive_fraction",
                              "Repetitive_fraction"],
                    "subclass": ["NA"] *2,
                    "order": ["NA"] *2,
                    "superfam": ["NA"] *2,
                    "naive_numele": ["NA"] *2,
                    "naive_avglen": ["NA"] *2,
                    "naive_bpsum": [
                        seqsizes_dict[species][sequid] - repetitive_bps["naive"],
                        repetitive_bps["naive"] ],
                    "tagged_numele": ["NA"] *2,
                    "tagged_avglen": ["NA"] *2,
                    "tagged_bpsum": [
                        seqsizes_dict[species][sequid] - repetitive_bps["tagged"],
                        repetitive_bps["tagged"]],
                }
                df_absolute_summary = pd.concat([df_absolute_summary,
                                        pd.DataFrame(new_rows)])

        # avgele length might not be interpreted as float dtype

        # compute non-overlapping length
        # removing all repeats with overlapping overshoots the correct bp
        # the value is ~twice as close, although smaller instead of larger
        # compute all true. then add overlapping between False and True?

        print(df_absolute_summary.to_markdown()) # debug
        df_absolute_summary.round(decimals=2).to_csv(
            'rpm5_prova.tsv', sep="\t", na_rep="NA")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Si es crida com a script:
if __name__ == '__main__':

    ### PARSE ARGUMENTS ###
    if len(sys.argv) < 3:
        sys.exit(f"Lacking RepeatMasker's output\n{sys.argv[0]} rm_sp1.out "+
                 "species1_name rm_sp2.out species2_name ...")

    files = []
    for i in range(1, len(sys.argv)):
        files.append(sys.argv[i])

    repeats = Repeat(files,
    # send genome length in bases (to compute % genome occupancy)
                     seqsizes_dict=dict(
                        Dcat={
                            "Dcatchr1": 827927493,
                            "Dcatchr2": 604492468,
                            "Dcatchr3": 431906129,
                            "Dcatchr4": 421735357,
                            "DcatchrX": 662930524,
                            "Scaffold|ctg": 330778139},
                        Dtil={
                            "Dtilchr1": 210826211,
                            "Dtilchr2": 210901293,
                            "Dtilchr3": 205634564,
                            "Dtilchr4": 158554970,
                            "Dtilchr5": 152324531,
                            "Dtilchr6": 125744329,
                            "DtilchrX": 434793700,
                            "Scaffold": 50989031}
                        ))



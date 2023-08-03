#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# repeatmasker_parser5.py
#
# 22 Jun 2023  <adria@molevol-OptiPlex-9020>

"""
INPUT:
======

 ## RepeatMasker.out modified.

Output given by RepeatMasker with a few modifications...

awk 'BEGIN{OFS="\t"} ; {if ( $16 ) { $16="True" } else { $16="False" } } 1' RM.out
There are rows with an asterisk as 16th field, and rows with only 15 fields.
Make sure all rows are of the same length (16 fields for all)

Moreover, make sure to correctly label minor scaffolds from chromosomes in column 5

Lastly, there can be few asterisks '*' in `position (left)`. These obstruct the
loading/reading of the file because that column is otherwise an integer.
Substitute them by zeroes.

 ## Index file/Sequid lengths

TSV file with pairs of "sequid" and "sequid_length"
"Sequids" are regexes to match in the RepeatMasker file.
ie. 'chr1' will match 'Species_chr1', and 'Scaffold' will match 'Scaffold_1527'

```idx
Species       whatever
Chr1          3600
Chr2          5000
...           ...
ChrN          100
Scaffolds     $(sum_len_scaffolds)
```

WARNING:
========

RepeatMasker output is massive. The file can weigh Gbytes. Pandas.DataFrame is
loaded into memory. If the file is too big, the script might be killed. My
solution for the whilebeing was using computers with a lot of RAM. Another
solution could be to employ the `dask` module for `pandas`.
"""

import sys
import math
# to read input's file size and creating folders
import os

# data (frames) handling
import pandas as pd, numpy as np
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

        # make sure path to provided files exists
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
            "STATUS: REVISING REGEX MATCHES BETWEEN PROVIDED CHR.INDEX AND RM.OUT")
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
#                print(" + Regex matching list:",
#                      list(df_species_and_regex["sequid"].unique()))
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

        # reclassify repeat types into 4 new neat columns
        df_concat = self.reclassification(df_concat)

        # Summarize dataframe into absolute bp count per pairs of Species,
        # sequid-regex for each type of repeat.
        df_naive = df_concat.groupby(
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
        # the function `remove_overlapping` outputs a df with a new column
        # called "algor_bpsum"
        df_absolute_summary["algor_bpsum"] = 0
        print_green("STATUS: COMPUTING BP ACCOUNTING FOR OVERLAPPING")
        # subset the dataframe for each sequid in each species
        for species in df_concat["Species"].unique():
            sequid_list = df_concat.loc[
                df_concat["Species"] == species,
                "sequid"].unique()
            for seq in sequid_list:
                # subset the dataframe for each sequid in each species
                d = df_concat.loc[
                    (df_concat["Species"]==species) &
                    (df_concat["sequid"]==seq), ]
                # instead of grouping by the loop var `seq`, group by
                # sequid_type because then ALL scaffolds will be together
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
                df_return = self.remove_overlapping(intervals).reset_index()
                df_return.insert(0, "Species", [species]*df_return.shape[0])
                df_return.insert(1, "sequid_type", [seqtype]*df_return.shape[0])
                df_absolute_summary = pd.merge(right=df_absolute_summary, left=df_return,
                                               # preserve right dataframe's keys
                                               how="right",
                                               on=["Species", "sequid_type",
                                                   "class", "subclass", "order",
                                                   "superfam"])
#                print(df_absolute_summary[["algor_bpsum_x", "algor_bpsum_y"]].info()) # debug
                df_absolute_summary["algor_bpsum_x"] = (
                df_absolute_summary["algor_bpsum_x"].fillna(0))
                df_absolute_summary["algor_bpsum_y"] = (
                    df_absolute_summary["algor_bpsum_x"] +
                    df_absolute_summary["algor_bpsum_y"])
                df_absolute_summary.drop(columns=["algor_bpsum_x"], inplace=True)
                df_absolute_summary.rename(columns={
                    "algor_bpsum_y": "algor_bpsum"}, inplace=True)
#                print(df_absolute_summary[["algor_bpsum"]].info()) # debug

        # print the sum of algorithmically overlapping bps
        self.how_much_overlapping(df_absolute_summary)

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
                    'algor': int(df_absolute_summary.loc[
                    (df_absolute_summary["Species"] == species) &
                    (df_absolute_summary["sequid_type"] == sequid),
                    "algor_bpsum"].sum()),
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
                    "algor_bpsum": [
                        seqsizes_dict[species][sequid] - repetitive_bps["algor"],
                        repetitive_bps["algor"]],
                }
                df_absolute_summary = pd.concat([df_absolute_summary,
                                        pd.DataFrame(new_rows)])

        # check whether dtypes have been recorded properly
        df_absolute_summary[["naive_numele", "naive_avglen",
                             "tagged_numele", "tagged_avglen"]] = (
        df_absolute_summary[["naive_numele", "naive_avglen",
                             "tagged_numele", "tagged_avglen"]].apply(
                                 pd.to_numeric, errors="coerce"))
        print_green(
            "\nSTATUS: INFO ABOUT INPUT/RECEIVED DATAFRAME")
        df_concat.info()
        df_absolute_summary.info()
        self.df_complete = df_absolute_summary
        self.df_main = df_concat

        # compute num_eles for 'Repetitive_fraction' (=total by sequid) class
        for species in df_absolute_summary["Species"].unique():
            sequid_list = df_absolute_summary.loc[
                df_absolute_summary["Species"]==species,
                "sequid_type"].unique()
            for seq in sequid_list:
                msp = df_absolute_summary["Species"]==species
                msq = df_absolute_summary["sequid_type"]==seq
                df_absolute_summary.loc[
                    (df_absolute_summary["class"]=="Repetitive_fraction") &
                    (msp) & (msq),
                "naive_numele"] = df_absolute_summary.loc[
                    (msp)&(msq), "naive_numele"].sum()
                df_absolute_summary.loc[
                    (df_absolute_summary["class"]=="Repetitive_fraction") &
                    (msp) & (msq),
                "tagged_numele"] = df_absolute_summary.loc[
                    (msp)&(msq), "tagged_numele"].sum()

        (self.df_groupby_seq_and_species,
         self.df_groupby_species,
         self.df_groupby_reptype_and_species) = (
            df_complete_summary_grouping(self.df_complete))

        # Write newly created dataframes to tabular files

        # store the classification of REs in case it has to be checked...
        self.check_classification(
            ).round(decimals=2).to_csv('repeat_reclassification.tsv',
            sep='\t', na_rep='NA',)

        # "refined" RM out
        self.df_main.round(decimals=2).to_csv(
            "repeat_df_main.tsv", sep="\t",
            na_rep="NA", index=False)

        # df with a complete summary of the raw RM output
        self.df_complete.round(decimals=2).to_csv(
            "repeat_df_complete.tsv", sep="\t",
            na_rep="NA", index=False)

        # Groupings of sequence/species
        self.df_groupby_seq_and_species.round(decimals=2).to_csv(
            "repeat_df_gby_seq_species.tsv", sep="\t",
            na_rep="NA", index=False)
        self.df_groupby_species.round(decimals=2).to_csv(
            "repeat_df_gby_species.tsv", sep="\t",
            na_rep="NA", index=False)
        self.df_groupby_reptype_and_species.round(decimals=2).to_csv(
            "repeat_df_gby_reptype_species.tsv", sep="\t",
            na_rep="NA", index=False)

        # another long __init__ function
        return None

    def how_much_overlapping(self,
        df: "input `df_absolute_summary` from __init__"):
        """
        Compute sum of bps for each method (naive_bpsum, tagged_bpsum and
        algor_bpsum) and percentages relative to naive (how much overlapping of
        repeats do we encounter)

        Print it all to terminal
        """

        for species in df["Species"].unique():
            d = df.loc[df["Species"] == species]
            naive_bpsum = d["naive_bpsum"].sum()
            print()
            print(" ----- How much overlapping of REs is present in", species,
                  "-----")
            print("+ Naive sum of basepairs:", naive_bpsum, "/ 100 %")
            print("+ Bps removing overlapping REs:",
                d["tagged_bpsum"].sum(), "/",
                round((d["tagged_bpsum"].sum()/naive_bpsum)*100, 1),
                "%")
            print("+ Bps algorithmically accounting for overlaps:",
                d["algor_bpsum"].sum(), "/",
                round((d["algor_bpsum"].sum()/naive_bpsum)*100, 1),
                "%")

        return None

    def reclassification(self, df: "input df_concat from __init__"):
        """
        Long function of many paragraphs, each reclassifying the default repeat
        type into four new columns, easily read and employed.

        Attention!

        Some classes are reclassified with 'Y==X' (the field Y exactly matches
        X) but others are reclassified with 'Y.str.contains(X)' (the field Y
        contains X string)
        """

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

        df.loc[
            df["default_repclass"].str.contains('Maverick'),
            ['class', 'subclass', 'order', 'superfam']
        ] = ('DNA', '2', 'Maverick', 'Maverick')

        # MITE or nMITE?
        df.loc[
            df["default_repclass"].str.contains('_MITE'),
            ['mite']
        ] = (True)

        df.loc[
            df["default_repclass"].str.contains('_nMITE'),
            ['mite']
        ] = (False)

        return df

    def remove_overlapping(self,
        intervals):
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
        # convert into pd.MultiIndex()
        idx = pd.MultiIndex.from_tuples(unique_class_tuples, names=[
                'class', 'subclass', 'order', 'superfam'])
        # create df
        df_resulting = pd.DataFrame(index=idx)
        df_resulting['algor_bpsum'] = 0

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
                    df_resulting.loc[reptype, "algor_bpsum"] += interval_length
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
                df_resulting.loc[reptype, "algor_bpsum"] += interval_length
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

    def check_classification(self):
        """
        Checks how the default classes have been re-classified
        """
        df = self.df_main

        print_green("\nSTATUS: CHECKING RE-CLASSIFICATION OF REPEAT 'TYPES'")
        # show the new classification
        # (drop rows sharing the same repeat class/subclass/etc.)
        df_new_classes = df.drop_duplicates(
            subset=['default_repclass', 'class',
                    'subclass', 'order', 'superfam'])
        # sort dataframe by classes (then subclass, etc.)
        df_new_classes = df_new_classes.sort_values(
            by=['class', 'subclass', 'order', 'superfam'])
        # reset index
        df_new_classes = df_new_classes.reset_index()
        # only print columns with repeat class names
        df_new_classes = df_new_classes[['class',
                                     'subclass', 'order',
                                     'superfam', 'mite',
                                     'default_repclass']]
        # print df in markdown format (requires tabulate package)
        print(df_new_classes.round(decimals=2).to_markdown(None))

        return df_new_classes

def df_complete_summary_grouping(df):
    """
    Group complete summary by certain columns/variables

    Reduces information (e.g. from rep-length by chromosome
    to whole species rep-length)
    """

    # compute groupby 'sequid' and 'species' values
    df_groupby_seq_and_species = df.groupby([
        "Species", "sequid_type", "class"])[
            ["naive_numele", "naive_bpsum", "tagged_numele",
             "tagged_bpsum", "algor_bpsum"]].agg([
                "sum"]).reset_index()
    # remove MultiIndex column names
    df_groupby_seq_and_species.columns = (
        df_groupby_seq_and_species.columns.droplevel(1))

    # compute groupby 'Species-only' values
    df_groupby_species = df.groupby([
        "Species", "class"])[
            ["naive_numele", "naive_bpsum", "tagged_numele",
             "tagged_bpsum", "algor_bpsum"]].agg([
                "sum"]).reset_index()
    # remove MultiIndex column names
    df_groupby_species.columns = (
        df_groupby_species.columns.droplevel(1))

    # compute groupby 'repeat-type' and 'species' values
    df_groupby_reptype_and_species = df.groupby([
        "Species", "class", "subclass", "order", "superfam"])[
            ["naive_numele", "naive_bpsum", "tagged_numele",
             "tagged_bpsum", "algor_bpsum"]].agg([
                "sum"]).reset_index()
    # remove MultiIndex column names
    df_groupby_reptype_and_species.columns = (
        df_groupby_reptype_and_species.columns.droplevel(1))

    return (df_groupby_seq_and_species,
            df_groupby_species,
            df_groupby_reptype_and_species)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class Plotting:

    def __init__(
          self, tsv_main,
          tsv_complete_summary,
          seqsizes_dict):
        """
        Read TSVs from `Repeat` class
        """
        # list of file-paths
        file_list = [tsv_main, tsv_complete_summary]

        # might be useful further on
        self.seqsizes_dict = seqsizes_dict

        # complete genome sizes
        gensizes_dict = dict()
        for species in seqsizes_dict.keys():
            gensizes_dict[species] = 0
            for sequid in seqsizes_dict[species].keys():
                gensizes_dict[species] += seqsizes_dict[species][sequid]
        self.gensizes_dict = gensizes_dict

        # make sure path to provided files exists
        for file in file_list:
            try:
                with open(file) as ft:
                    pass
            except:
                sys.exit(f'ERROR: Cannot read the provided path: {file}')

        # read provided TSVs
        self.df_main = pd.read_table(tsv_main).fillna("NA")
        df = pd.read_table(tsv_complete_summary)

        # fill NaN with either zeroes or "NA" depending on column
        df.loc[:, "Species":"superfam"] = (
            df.loc[:, "Species":"superfam"].fillna("NA") )
        df.loc[:, "naive_numele":"algor_bpsum"] = (
            df.loc[:, "naive_numele":"algor_bpsum"].fillna(0) )
        self.df_complete_summary = df.copy()

        # create a dict to subset df by kinds of chromosome
        self.regex_to_genome_subset = {  # pair regex with the name of the genome subset
            "ChrX": r"chrX", "Autosomes": r"chr\d", "Genome": r".*"}
            # contains `chr[[digit]]`, `chrX` or `anything`

        # create a new column, where each class has a given colour
        # For available shades, take a look at:
        # https://matplotlib.org/stable/gallery/color/named_colors.html
        self.colours_by_class = {
            "DNA": "mediumslateblue",
            "Other": "orange",
            "Retrotransposon": "salmon",
            "Unknown": "olive",
            "Tandem_repeat": "yellow",
            "Nonrepetitive_fraction": "plum"
        }
        self.species_palette = {
            "Dcat": "#eb8f46",
            "Dtil": "#30acac",
            "Dsil": "mediumseagreen",
        }

        return None

    def histogram_divergence(
        self, df, write_to_filepath, title,
        density: bool,
        hue: "colname to hue by (e.g. Species)" =False):
        """
        Draw a distribution of divergence from main df or a subset of main df

        density: Density distribution (True) or absolute distribution (False)

        df colname with divergence values: 'perc_divg'
        df colname with colours for each bin: 'colours_class'
        `write_to_filepath`: save PNG figure to this path
        """
        # specify amount of bins (as many as max value?)
        amount_bins = math.floor(df["perc_divg"].max()) + 1
        if hue and density:
            sns.histplot(data=df, x='perc_divg',
                         # relative percentage normalized by "Species"
                         # instead of the whole dataset.
                         common_norm=False, stat="percent",
                         # amount of bins and draw KDE
                         bins=int(amount_bins), kde=True,
                         # dodge bars of multiple species
                         multiple="dodge", hue=hue,
                         palette=self.species_palette )
            plt.ylabel("Density")
        elif hue:
            sns.histplot(data=df, x='perc_divg',
                         # amount of bins and draw KDE
                         bins=int(amount_bins), kde=True,
                         # dodge bars of multiple species
                         multiple="dodge", hue=hue,
                         palette=self.species_palette )
            plt.ylabel("Count")
        else:
            ax = sns.histplot(data=df, x='perc_divg',
                              bins=int(amount_bins))

        plt.title(title, fontsize=12)
        plt.suptitle("Divergence from consensus sequence", fontsize=12)
        plt.xlabel("% Divergence")
        plt.savefig(write_to_filepath, dpi=300)
        plt.close('all')

        return None

    def histogram_both_species_and_reptype_pairs(
        self, folder=None, relative=True,
        filter_overlap_label=False,
        all_species=True):
        """
        Create histograms overlapping on the same plot/axis.

        If relative, compare the shape of divergence distributions between species.
        Else not relative, compare absolute counts in divergence bins.

        folder: String with folder name. By default does not create a folder.

        filter_overlap_label: Filter repeats labelled as overlapping. Remove
            them from the analysis.

        relative: Whether the plot should have a relative or absolute axis per
            species. Relative axis are "density-like" (density of each
            divergence bin per species) while absolute axis display the count of
            REs per bin. Relative axis allows comparison of distribution shapes,
            while absolute allows the accumulation and excesses of REs.
        """
        d1 = self.df_main.copy()
        unique_class = list(set(list(zip(list(d1["class"])))))
        # unique subclasses, orders and superfams from DNA and RT only, please
        d2 = d1.loc[d1["class"].isin(["DNA", "Retrotransposon"])]
        unique_subclass = list(set(list(zip(
            list(d2["class"]), list(d2["subclass"])))))
        unique_order = list(set(list(zip(
            list(d2["class"]), list(d2["subclass"]), list(d2["order"])))))
        unique_superfam = list(set(list(zip(list(d2["class"]),
            list(d2["subclass"]), list(d2["order"]), list(d2["superfam"])))))

        reptype_levels = ["class", "subclass", "order", "superfam"]

        if filter_overlap_label:
            # remove overlapping repeats from dataframe
            d1 = d1.loc[d1["overlapping"] == False]

        # for each reptype in the concatenated list with all different levels...
        for reptype in unique_class+unique_subclass+unique_order+unique_superfam:
            # filter dataframe and create histogram from it
            d2 = d1.copy()
            for i in range(0, len(reptype)):
                mask = d2[reptype_levels[i]] == reptype[i]
                d2 = d2.loc[mask]
            if all_species:
                if folder:
                    if not os.path.exists(folder):
                        os.makedirs(folder)
                    write_to_filepath = (folder+"/divhist_AllSp_"+
                        "_".join(reptype)+'.png')
                else:
                    write_to_filepath = ("divhist_AllSp_"+
                        "_".join(reptype)+'.png')

                self.histogram_divergence(
                    d2, write_to_filepath,
                    title="; ".join(reptype) + "; [filter "+
                        str(filter_overlap_label)+"]",
                    hue="Species",
                    density=relative)  # default True

            else: # not all species (every species individually)
                for species in d2["Species"].unique():
                    mask = d2["Species"] == species
                    d3 = d2.loc[mask]
                    if folder:
                        if not os.path.exists(folder):
                            os.makedirs(folder)
                        write_to_filepath = (folder+"/divhist_"+
                            species+"_"+"_".join(reptype)+'.png')
                    else:
                        write_to_filepath = ("divhist_"+
                            species+"_"+"_".join(reptype)+'.png')

                    self.histogram_divergence(
                        d3, write_to_filepath,
                        title="; ".join(reptype) + "; [filter "+
                            str(filter_overlap_label)+"]",
                        hue=False,
                        density=False)  # default True

        return None

    def boxplots_df_main(self, folder=None):
        """
        Draw boxplots of data from input dataframe (self.df_main) by
          + X variable (e.g. SW score, repeat length...)
          + Y category (e.g. repeat classes)
          + Z hue (e.g. species, sequid...)
        """
        df = self.df_main

        # pair x-variables with figure labels and titles
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
                        # determine species hue
                        palette=self.species_palette,
                        meanprops=dict(color='green',
                                       marker='*',
                                       markeredgecolor='black',
                                       markersize='7'))
            plt.title(fig_titles[xvar])
            plt.subplots_adjust(left=0.25, bottom=0.1, right=0.95, top=0.91)
            plt.xlabel(fig_xlabel[xvar])
            if folder:
                if not os.path.exists(folder):
                    os.makedirs(folder)
                plt.savefig(folder+'/matplotlib_boxplots_'+str(xvar)+'.png', dpi=300)
            else:
                plt.savefig('matplotlib_boxplots_'+str(xvar)+'.png', dpi=300)
            plt.close('all')

        return None

    def histogram_stacked_absolute(self):
        """
        INPUT
        =====
        `self.df_complete_summary`: Complete summary df
        `self.colours_by_class`: Colours applied to each RE class

        OUTPUT
        ======
        PNG file with the plotted histogram

        More info/ references:
        https://www.pythoncharts.com/matplotlib/stacked-bar-charts-labels/
        """
        # create a df copy of `df_complete_summary`
        df = self.df_complete_summary.loc[:]
        # create a dict to subset df by kinds of chromosome
        rgensub = self.regex_to_genome_subset
        # remove repetitive fraction (sum of all repeats)
        df = df.loc[(df["class"]=="Repetitive_fraction")==False]
        # we are interested in sorting classes manually into the following order
        ordering_classes = {
            "DNA": 1, "Retrotransposon": 2, "Tandem_repeat": 3,
            "Other": 4, "Unknown": 5, "Nonrepetitive_fraction": 6}
        # we are also interested in secondarily sorting rows by 'algor_bpsum'
        # (sum of basepairs for that entry/row in the table)
        df = df.sort_values(["class", "algor_bpsum"], key=lambda x: x.replace(ordering_classes),
                            ascending=[True, True])
        # kind of sort = mergesort: to preserve relative order of previous sort
        # reference:
        # https://stackoverflow.com/questions/47440621/numpy-argsort-behavior-for-equal-numbers
        df = df.sort_values(["Species", "sequid_type"], kind="mergesort").reset_index(drop=True)

        # create a dict to subset df by kinds of chromosome
        regex_to_gensubs = {  # pair regex with the name of the genome subset
            "Allosome": r"chrX", "Autosomes": r"chr\d", "Genome": r".*"}
            # contains `chr[[digit]]`, `chrX` or `anything`

        # remove edgecolor from bars
        fig, ax = plt.subplots(figsize=(8,14))
        # `zorder` sends the grid to background
        # it has to be also specified for `ax.bar()`
        ax.grid(axis='y', color="darkslategrey", linestyle="dashed",
                linewidth=2, zorder=0)

        bottom = np.zeros(3 * len(df["Species"].unique()))
        for rep_class in df["class"].unique():
            d1 = df.loc[df["class"] == rep_class]
            bp_values = dict()
            for rep_superfam in d1["superfam"].unique():
                d2 = d1.loc[d1["superfam"] == rep_superfam]
                for species in df["Species"].unique():
                    d3 = d2.loc[d2["Species"] == species]
                    for gsub_key in rgensub.keys():
                        # subset df where sequid type contains a regex
                        d4 = d3.loc[d3["sequid_type"].str.contains(
                            rgensub[gsub_key])]
                        xaxis = str(species) + "_" + str(gsub_key)
                        # sum of basepairs for the triple-subsetted df
                        bp_values[xaxis] = d4["algor_bpsum"].sum()

                # create the segment bar
                print("Creating segment bar for", # debug
                      rep_class, rep_superfam, list(bp_values.values())) # debug
                ax.bar(x=list(bp_values.keys()), height=list(bp_values.values()),
                       bottom=bottom, label=rep_class,
                       linewidth=0.4, edgecolor='black',
                       color=self.colours_by_class[rep_class],
                       zorder=3)
                # sum to bottom
                bottom += np.array(list(bp_values.values()))

        ax.set_title("Absolute occupancy (in Gb) by species and genome subset")
        plt.suptitle("Each black division segregates classes by superfamilies")
        handles, labels = plt.gca().get_legend_handles_labels()
        labels, ids = np.unique(labels, return_index=True)
        handles = [handles[i] for i in ids]
        ax.legend(handles, labels, loc='best')
        #plt.subplots_adjust(left=0.3, bottom=0.18, right=0.59, top=0.95)
        plt.ylabel("Gb", fontsize=12)
        plt.xticks(rotation=90)
        #plt.xlabel("")
        plt.savefig('matplotlib_histo_stacked_absolute.png', dpi=300)

        return None

    def col_to_categorical_palette(self, categories: "list of strings",
        shade: "what kind of shade will be used (red, blue, salmon...)"):
        """
        For available shades, take a look at:
        https://matplotlib.org/stable/gallery/color/named_colors.html

        + categories
        A list of strings, each string symbolizing a category that
        will be paired to a colour of the seaborn palette.

        Returns dict with keys being categories and values colours.
        """
        # if there is a single categorical entry do not
        # return a palette, but a single colour for all rows
        if len(categories) == 1:
            # split and remove 'dark:' or 'light:' seaborn prefix
            dict_cols = {categories[0]: shade}

        else:
            shade = "light:" + shade # use the light shade
            list_cols = sns.color_palette(shade, len(categories))
            list_cols.reverse() # send the "dark:" or "light:" colours to the end
            dict_cols = dict(zip(
                categories,  # keys
                list_cols))  # values

        return dict_cols

    def pie_chart(self, df: "unfiltered or filtered `df_complete_summary`",
        grouptype: "reptype colname which will groupby() the df",
        threshold: "% by which to cut abundant from scarce" =5,
        folder: "Folder to save the figure in" =None):
        """
        INPUT
        =====

        + df
        `df_complete_summary`. Could be filtered (only contains "DNA" class,
        only contains "LINEs", etc.

        + grouptype
        Groupby() REs by this categorical level ("class", "subclass", etc.)

        + threshold
        Pool REs below the specified percentage into unspecified `Others`

        + folder
        If specified, checks whether a folder with this name exists; if not, it
        creates it; afterwards, it stores PNGs in that folder.

        OUTPUT
        ======

        PNG plot/figure/piechart on the working directory/in the specified
        folder (new folder would be created if there didnt exist one)
        """
        # Create the filename where fig will be written to
        if folder:
            if not os.path.exists(folder):
                os.makedirs(folder)
            fig_filename = (str(folder) + "/matplotlib_piecharts_by" +
                            str(grouptype) + ".png")
        else:
            fig_filename = "matplotlib_piecharts_by"+str(grouptype)+".png"
        # create a dict to subset df by kinds of chromosome
        rgensub = self.regex_to_genome_subset
        # create a dict to find following classification item
        # (e.g. subclass for class, superfam for order, order for subclass)
        below_type = {"class": "subclass", "subclass": "order",
                      "order": "superfam"}
        # `class_col` will be used to colour subclass, order and superfams by
        # the class colour in `self.colours_by_class`
        class_col = str(self.colours_by_class[
            list(df["class"].unique())[0]])

        Species = list(df["Species"].unique())
        Gensub  = list(rgensub.keys())
        list_dfs_to_pie = list()
        # remove repetitive fraction overall total/sum
        df = df.loc[df["class"] != "Repetitive_fraction"]

        # create a list with species * gensub dataframes for the pies
        for i in range(0, len(Species)):
            for j in range(0, len(Gensub)):
                sp = Species[i]
                gs_key   = Gensub[j]
                gs_regex = rgensub[gs_key]
                # subset main df for each of the i*j pies
                d1 = df.loc[ (df["Species"]==sp) &
                       (df["sequid_type"].str.contains(gs_regex))]
                d1 = d1.groupby(grouptype)["algor_bpsum"].agg(["sum", "count"])
                d1 = d1.sort_values(["sum"], ascending=False).reset_index()
                # compute relative occupancy of the pie
                d1["percent"] = round( (d1["sum"] /d1["sum"].sum()) *100, 1)
                d1["megabases_span"]=round(d1["sum"]/10**6, 1)
                subplot_title = str(sp) +"_"+ str(gs_key)
                list_dfs_to_pie += [[d1, i, j, subplot_title]]

        # checking which types are in all species below threshold of occupancy
        # first, get genomic dfs for all species with [::3] slice
        ##print([pie for pie in list_dfs_to_pie[2::3]]) # debug
        df_check_threshold = pd.concat([pie[0] for pie in list_dfs_to_pie[2::3]])
        # remove rows above percent (ie occupancy) threshold
        df_check_threshold = df_check_threshold.loc[
            df_check_threshold["percent"] < int(threshold)]
        # make sure each reptype is below threshold in all species
        for reptype in df_check_threshold[grouptype].unique():
            df_check = df_check_threshold.loc[
                df_check_threshold[grouptype] == reptype]
            # if n_rows != n_species
            if df_check.shape[0] != len(Species):
                # remove these rows
                df_check_threshold = df_check_threshold.loc[
                    df_check_threshold[grouptype] != reptype]

        # if there is more than 1 reptype that fulfills previous requirements:
        # group all reptypes into a "bin" ("Group_Others")
        binning = list(df_check_threshold[grouptype].unique())

        # subsetting dfs by `binning` and plotting them
        fig, axs = plt.subplots(nrows=3, ncols=len(Species),
                             figsize=(18,12))
        for pie in list_dfs_to_pie:
            # unpack pie list
            df, i, j, subplot_title = pie
            df = df.copy() # otherwise warns that it refers to df in the tuple
            df_below_thresh = df.loc[df[grouptype].isin(binning)]
            # remove below threshold reptypes from pie df
            df = df.loc[ (df[grouptype].isin(binning)) == False]
            # create the palette on first sweep of the loop
            if i==0 and j==0:
                dict_cols = self.col_to_categorical_palette(
                    categories=list(df[grouptype].unique()),
                    shade=class_col)
            # create a palette of colour for REs above threshold;
            # below threshold will be coloured all grey/white...
            if grouptype == "class":
                df["colours"] = "k"
                for key, val in self.colours_by_class.items():
                    df.loc[df["class"] == key, "colours"] = val
            else:
                for reptype, col in dict_cols.items():
                    # we cannot use the simpler method `column.replace(dict)`
                    # because some values within the dictionary are RGB values
                    # stored in tuples and numpy does not allow {str: tuple}
                    mask = df[grouptype] == reptype
                    df.loc[mask, "colours"] = pd.Series([
                        col for x in range(df.loc[mask].shape[0])],
                        index=df.loc[mask].index)

            # add new row summing every RE below threshold
            new_row = {
                grouptype: ["Remaining"],
                "colours": ["lightgrey"],
                "sum": [df_below_thresh["sum"].sum().round(1)],
                "count": [df_below_thresh["count"].sum()],
                "percent": [df_below_thresh["percent"].sum().round(1)],
                "megabases_span": [df_below_thresh["megabases_span"].sum(
                    ).round(1)],
            }
            df = pd.concat([df, pd.DataFrame(new_row)])

            # append count of types and sum of bps to the label of the pie
            df["label"] = (
                df[grouptype].astype(str) +"_"+
                df["count"].astype(str) +"_"+
                df["percent"].astype(str) +"_"+
                df["megabases_span"].astype(str) +"Mb" )

            axs[j, i].pie(df["percent"], labels=df["label"],
                            colors=df["colours"])
            axs[j, i].set_title(subplot_title)

        fig.suptitle("By "+str(grouptype))
        plt.savefig(fig_filename, dpi=300)

        return None

    def tabulate_divergence(self,
        unique_reptypes: "list of the allowed combinations of repeat-types"):
        """
        Create a table that groups by unique repeat types. Each value in the
        list is an ordered (from class, subclass, order to superfam) sublist of
        combinations of repeat types, e.g.

        [["DNA", "1", "TIR", "hAT"], ["DNA", "1", "TIR", "TcMar"], etc.]

        For each combination, locate in main dataframe and compute relative (%)
        genome occupancy and absolute (Mb) genome occupancy. Then, compute mean,
        median and standard deviation for the column of repeat divergence.

        Giving a list of unique combinations of repeat types circumvents the
        problem of duplicated names, i.e. "NA" from both "DNA" and
        "Retrotransposons" being considered a single group (and pooling values).
        Instead pass as a variable each type you would like to consider:

        [["DNA", "NA"], ["RT", "NA"]]

        output df. columns:

        Species | reptype | overlap (True, False) | measures (median, %, Mb, etc.)
        """
        # unique_reptypes must be sorted as follows:
        reptype_levels = ["class", "subclass", "order", "superfam"]
        # shorten call to main df
        df = self.df_main.copy()
        df_summ = self.df_complete_summary.copy()
        # make sure all items on the `unique_reptypes` list are of the same len.
        # create a list with the len of each unique_reptypes,
        # then set (remove duplicates) and make sure that its len is 1
        # (only one unique length value for all sublists)
        if len(list(set([len(x) for x in unique_reptypes]))) != 1:
            print("The given unique_reptypes might be of different length")
            return None

        # init. dict. where data will be appended
        # and later on will be the seed for a pd.DataFrame
        new_df = {"Species": [], "overlap": [], "perc_genom": [],
            "abs_occup_Mb": [], "median_divg": [], "mean_divg": [],
                  "stdev_divg": [], "max_divg": []}
        # create a column for each included repeat type (maybe only class, maybe
        # class and subclass, etc)
        for i in range(0, len(unique_reptypes[0])):
            new_df[reptype_levels[i]] = []
        for reptype in unique_reptypes:
            # filter df and df_summ with all levels of reptype
            d1 = df.copy()
            d1_summ = df_summ.copy()
            for i in range(0, len(reptype)):
                mask = d1[reptype_levels[i]] == reptype[i]
                d1 = d1.loc[mask]
                mask = d1_summ[reptype_levels[i]] == reptype[i]
                d1_summ = d1_summ.loc[mask]
            for species in d1["Species"].unique():
                for overlap in [True, False]:
                    # mask df_summ only by species
                    # (overlapping will be filtered in if-branches,
                    # either selecting tagged or algor bp)
                    mask = d1_summ["Species"] == species
                    d2_summ = d1_summ.loc[mask]

                    # compute a few last variables before creating rows
                    if overlap:
                        # because here we remove overlapping, use
                        # the basepairs from "tagged" column
                        repeat_bp = d2_summ["tagged_bpsum"].sum()
                        # mask df by species and overlapping
                        mask = (d1["Species"] == species) & (
                            d1["overlapping"] == False)
                        d2 = d1.loc[mask]
                        overlap = "Excluded"
                    else:
                        # because here we dont remove overlap,
                        # use basepairs from "algor" column
                        repeat_bp = d2_summ["algor_bpsum"].sum()
                        # mask df only by species (include
                        # both overlapping True and False)
                        mask = (d1["Species"] == species)
                        d2 = d1.loc[mask]
                        overlap = "Included"
                    # compute relative occupancy, in %
                    gensize = self.gensizes_dict[species]
                    repeat_perc = round((repeat_bp /gensize)*100, 2)
                    # `repeat_bp` from basepairs to Mbps (division by 10**6)
                    repeat_bp = round(repeat_bp/(10**6), 1)

                    # store values `reptype`, `species` and `overlap`
                    # in a new row
                    new_df["Species"].append(species)
                    # pair first value in list with "class", etc.
                    for i in range(0, len(reptype)):
                        new_df[reptype_levels[i]].append(reptype[i])
                    new_df["overlap"].append(overlap)
                    new_df["abs_occup_Mb"].append(repeat_bp)
                    new_df["perc_genom"].append(repeat_perc)
                    new_df["median_divg"].append(
                                d2["perc_divg"].median())
                    new_df["mean_divg"].append(
                                d2["perc_divg"].mean())
                    new_df["stdev_divg"].append(
                                d2["perc_divg"].std())
                    new_df["max_divg"].append(
                                d2["perc_divg"].max())

        new_df = pd.DataFrame(new_df)
        # sort `new_df` so the same reptype is paired across species
        # (to facilitate comparisons across species)
        sorting_list = []
        for i in range(0, len(unique_reptypes[0])):
            sorting_list.append(reptype_levels[i])
        sorting_list.append("Species")
        sorting_list.append("overlap")
        new_df = new_df.sort_values(sorting_list
                            ).reset_index(drop=True)
        # rearrange cols (send repeat type columns to beginning)
        cols = new_df.columns.tolist()
        # send the first seven columns behind the last column
        send_to_end = [cols.pop(i) for i in [0]*8]
        cols = cols + send_to_end
        new_df = new_df[cols]

        return new_df.round(2)

    def sliding_windows(
          self,
          # input these from 'self.seqsizes_dict'
          # (single entry or looping all entries)
          sequid_name: str, sequid_len: int, species: str,
          fig_title: str,
          # it is possible to input a pre-filtered dataframe by reptypes
          df: "main dataframe, filtered or unfiltered",
          folder: str = None,
          window_size: int = 100000, jump_size: int = 0,):
        """
        INPUT
        ========

        `sequid_name`, `sequid_len` and `species` are obtained from
        `seqsizes_dict`. Input a single entry or loop over all dict. values.

        Does not accept scaffold amalgams.

        + sequid_name: Name of the sequid to label the plot.
        + sequid_len: Length (in basepairs)  of the sequid.
        + species: Name of the species, to label the plot.
        + df: Main dataframe from the __init__ function. Can be filtered by repeat
            types (class == DNA, superfam == hAT, etc) to compute sliding
            windows for a subset of repeat groups.
        + folder: In which folder should the plot be written in.
        + window_size: Size of the sliding windows. Defaults to 100 kbps.
        + jump_size: Gap between windows that will be skipped. Defaults to zero.

        OUTPUT
        =========

        Creates a figure and writes it to a PNG.
        """
        # filter main dataframe by sequid
        mask = (df["Species"] == species) & (df["sequid"] == sequid_name)
        df = df.loc[mask]
        # init loop variables
        window_begin = 0
        genome_values = {
            # window coordinates and RE count
            "window": [], "ele_count": [],
            # median and st. dev. repeat length
            "replen_med": [], "replen_std": [], "replen_mean": [],
            # median and st. dev. divergence from consensus sequence
            "divg_med": [], "divg_std": [], "divg_mean": []}

        # slide and compute windows across the genome
        while window_begin < sequid_len:
            window_end = window_begin + window_size
            # find all REs that start inside the interval
            # (value of "begin" column is inside [begin, end] interval)
            mask = (window_begin < df["begin"]) & (df["begin"] < window_end)
            # filter out from dataframe
            d1 = df.loc[mask]
            # compute the center coordinate of the window
            window_center = (window_end + window_begin) / 2
            # store window values per loop cycle
            genome_values["window"].append(window_center)
            # shape[0] returns the number of rows in df
            genome_values["ele_count"].append(d1.shape[0])
            if d1.empty:
                genome_values["replen_med"].append(0)
                genome_values["replen_mean"].append(0)
                genome_values["replen_std"].append(0)
                genome_values["divg_med"].append(0)
                genome_values["divg_mean"].append(0)
                genome_values["divg_std"].append(0)
            else:
                genome_values["replen_med"].append(d1["replen"].median())
                genome_values["replen_mean"].append(d1["replen"].mean())
                genome_values["divg_med"].append(d1["perc_divg"].median())
                genome_values["divg_mean"].append(d1["perc_divg"].mean())
                # windows with only one repeat cannot have std()
                if d1.shape[0] == 1:
                    genome_values["divg_std"].append(0)
                    genome_values["replen_std"].append(0)
                else:
                    genome_values["divg_std"].append(d1["perc_divg"].std())
                    genome_values["replen_std"].append(d1["replen"].std())
            # slide the window forward for the next cycle
            window_begin += window_size + jump_size

        # create 3x1 subplots that share X (window coordinates)
        fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, sharex=True)
        # plot() method creates the center line (median)
        # while fill_between() method fills in the range of deviation
        # "ele_count" does not have deviation
        ax1.plot(genome_values["window"], genome_values["ele_count"],
                 label="Count of elements", lw=1, color="C0")
        ax1.set_title("Count of elements")

        ax2.plot(genome_values["window"], genome_values["replen_med"],
                 label="Median element length", lw=1, color="C1")
        zip_mean_and_std = list(zip(genome_values["replen_mean"],
                                    genome_values["replen_std"]))
        lower = [t[0] - t[1] for t in zip_mean_and_std]
        upper = [t[0] + t[1] for t in zip_mean_and_std]
        ax2.fill_between(genome_values["window"], lower, upper, facecolor="C1",
                         alpha=0.5, label="1 sigma range")
        ax2.set_title("Median element length")

        ax3.plot(genome_values["window"], genome_values["divg_med"],
                 label="Median divergence", lw=1, color="C2")
        zip_mean_and_std = list(zip(genome_values["divg_mean"],
                                    genome_values["divg_std"]))
        lower = [t[0] - t[1] for t in zip_mean_and_std]
        upper = [t[0] + t[1] for t in zip_mean_and_std]
        ax3.fill_between(genome_values["window"], lower, upper, facecolor="C2",
                         alpha=0.5, label="1 sigma range")
        ax3.set_title("Median divergence")

        plt.xlabel("Chromosome coordinates of " + sequid_name)
        fig.suptitle(fig_title + "; W.size=" + str(window_size), fontsize="medium")
        # create the folder if it is requested and doesnt exist.
        if folder:
            if not os.path.exists(folder):
                os.makedirs(folder)
            plt.savefig(folder+"/sliding_"+sequid_name+".png", dpi=300)
        else:
            plt.savefig("sliding_"+sequid_name+".png", dpi=300)
        # Figures are retained unless explicitly closed and may consume too much
        # memory. Close them after being written to disk.
        plt.close()

        return None


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Si es crida com a script:
if __name__ == '__main__':

    ### PARSE ARGUMENTS ###
    if len(sys.argv) < 3:
        sys.exit(f"Lacking RepeatMasker's output\n{sys.argv[0]} [type:r|p]\n"+
                 f"For r[epeat-parsing] option:\n{sys.argv[0]} r sp1_RM.out "+
                 "sp1_sequidsizes.idx sp2_RM.out sp2_sequidsizes.idx ..."+
                 "\n"+
                 f"For p[lotting] option:\n{sys.argv[0]} p main.tsv summary.tsv "
                 "sp1_sequidsizes.idx sp2_sequidsizes.idx ...\n\n"+
                 "The argv. `sequidsizes.idx` is an index portraying the "+
                 "length of each sequid; TSV with two columns, one sequid "+
                 "name and the other sequid length.")

    if sys.argv[1]=='r':
        print("Selected 'r' option (repeat-parsing)\n")
        files = []
        seqsizes_dict = dict()

        for i in range(2, len(sys.argv), 2):
            # store a list of file-paths as strings
            files.append(sys.argv[i])
            # from idx tabular file to python dict...
            df = pd.read_table(
                sys.argv[i+1],
                sep='\s+', index_col=0)
            seqsizes_dict[df.index.name]=df.iloc[:,0].to_dict()

        repeats = Repeat(files,
        # send genome length in bases (to compute % genome occupancy)
                    seqsizes_dict=seqsizes_dict)

    elif sys.argv[1]=='p':
        print("Selected 'p' option (plotting)\n")
        # load TSV path created by previous Repeat class...
        tsv_main_path = sys.argv[2]
        tsv_complete_summary_path = sys.argv[3]

        seqsizes_dict = dict()
        for i in range(4, len(sys.argv)):
            # from idx tabular file to python dict...
            df = pd.read_table(
                sys.argv[i],
                sep='\s+', index_col=0)
            seqsizes_dict[df.index.name]=df.iloc[:,0].to_dict()

        plots = Plotting(
            tsv_main_path,
            tsv_complete_summary_path,
            seqsizes_dict)

        # plots from main df (highly computationally consuming)
        for species in plots.seqsizes_dict.keys():
            for sequid in plots.seqsizes_dict[species].keys():
                sl = plots.seqsizes_dict[species][sequid]
                df = plots.df_main.copy()
                plots.sliding_windows(sequid_name=sequid,
                                      sequid_len=sl,
                                      species=species,
                                      df=df, fig_title="All repeats",
                                      folder="SW_All_Repeats")

                # Repeat above sliding windows but for specific RE types
                d1 = df.loc[df["class"] == "DNA"]
                plots.sliding_windows(sequid_name=sequid,
                                      sequid_len=sl,
                                      species=species,
                                      df=d1, fig_title="Filtered by DNA",
                                      folder="SW_class_DNA")
                d1 = df.loc[df["class"] == "Retrotransposon"]
                plots.sliding_windows(sequid_name=sequid,
                                      sequid_len=sl,
                                      species=species,
                                      df=d1, fig_title="Filtered by RT",
                                      folder="SW_class_RT")
                d1 = df.loc[df["superfam"] == "hAT"]
                plots.sliding_windows(sequid_name=sequid,
                                      sequid_len=sl,
                                      species=species,
                                      df=d1, fig_title="Filtered by DNA-hAT",
                                      folder="SW_class_DNA_superfam_hAT")
                d1 = df.loc[df["superfam"] == "TcMar"]
                plots.sliding_windows(sequid_name=sequid,
                                      sequid_len=sl,
                                      species=species,
                                      df=d1, fig_title="Filtered by DNA-TcMar",
                                      folder="SW_class_DNA_superfam_TcMar")
                d1 = df.loc[df["superfam"] == "Gypsy"]
                plots.sliding_windows(sequid_name=sequid,
                                      sequid_len=sl,
                                      species=species,
                                      df=d1, fig_title="Filtered by RT-Gypsy",
                                      folder="SW_class_RT_superfam_Gypsy")
                d1 = df.loc[df["order"] == "LINE"]
                plots.sliding_windows(sequid_name=sequid,
                                      sequid_len=sl,
                                      species=species,
                                      df=d1, fig_title="Filtered by RT-LINE",
                                      folder="SW_class_RT_order_LINE")

        plots.histogram_both_species_and_reptype_pairs(
            folder="Divg_relative_density_comparison", relative=True,
            all_species=True, filter_overlap_label=False)
        plots.histogram_both_species_and_reptype_pairs(
            folder="Divg_relative_density_comparison_remove_overlap",
            relative=True, all_species=True, filter_overlap_label=True)
        plots.histogram_both_species_and_reptype_pairs(
            folder="Divg_absolute_counts_perSp", relative=False,
            all_species=True, filter_overlap_label=False)
        plots.histogram_both_species_and_reptype_pairs(
            folder="Divg_individual_species",
            all_species=False, filter_overlap_label=False)
        plots.boxplots_df_main(folder="Boxplots")

        # tabulate main df divergence values by repeat type levels
        # (highly computationally consuming)
        df = plots.df_main.copy()
        unique_reptypes = list(set(list(zip(list(df["class"])))))
        plots.tabulate_divergence(unique_reptypes=unique_reptypes).to_csv(
            "divg_table_byclass.tsv", sep="\t", na_rep="NA", index=False)
        unique_reptypes = list(set(list(zip(
            list(df["class"]), list(df["subclass"])))))
        plots.tabulate_divergence(unique_reptypes=unique_reptypes).to_csv(
            "divg_table_bysubclass.tsv", sep="\t", na_rep="NA", index=False)
        unique_reptypes = list(set(list(zip(
            list(df["class"]), list(df["subclass"]), list(df["order"])))))
        plots.tabulate_divergence(unique_reptypes=unique_reptypes).to_csv(
            "divg_table_byorder.tsv", sep="\t", na_rep="NA", index=False)
        unique_reptypes = list(set(list(zip(list(df["class"]),
            list(df["subclass"]), list(df["order"]), list(df["superfam"])))))
        plots.tabulate_divergence(unique_reptypes=unique_reptypes).to_csv(
            "divg_table_bysuperfam.tsv", sep="\t", na_rep="NA", index=False)

        # plots from summary df
        plots.histogram_stacked_absolute()
        plots.pie_chart(df=plots.df_complete_summary,
                        grouptype="class")
        for reptype in ["subclass", "order", "superfam"]:
            df = plots.df_complete_summary.copy()
            df_dna = df.loc[df["class"] == "DNA"]
            plots.pie_chart(df=df_dna, folder="DNA",
                            grouptype=reptype)
            df_retro = df.loc[df["class"] == "Retrotransposon"]
            plots.pie_chart(df=df_retro, folder="Retrotransposon",
                            grouptype=reptype)

    else:
        print("Send either 'r' or 'p' as first argv "+
              "to specify in which mode the script should run")



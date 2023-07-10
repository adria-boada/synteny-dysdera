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
"""

import sys

# data (frames) handling
import pandas as pd, numpy as np
import os # compute file-size

# plotting
import seaborn as sns, matplotlib.pyplot as plt
import matplotlib as mpl

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
        # covert into pd.MultiIndex()
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

    # compute groupby 'Species-only' values
    df_groupby_species = df.groupby([
        "Species", "class"])[
            ["naive_bpsum", "tagged_bpsum", "algor_bpsum"]].agg([
                "count", "sum"]).reset_index()

    # compute groupby 'repeat-type' and 'species' values
    df_groupby_reptype_and_species = df.groupby([
        "Species", "class", "subclass", "order", "superfam"])[
            ["naive_bpsum", "tagged_bpsum", "algor_bpsum"]].agg([
                "count", "sum"]).reset_index()

    return (df_groupby_seq_and_species,
            df_groupby_species,
            df_groupby_reptype_and_species)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class Plotting:

    def __init__(self,
                 tsv_main,
                 tsv_complete_summary,
                 seqsizes_dict,):
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
        self.df_complete_summary = pd.read_table(tsv_complete_summary).fillna("NA")

        # create a new column, where each class has a given colour
        self.colours_by_class = {
            "DNA": "lightblue",
            "Other": "orange",
            "Retrotransposon": "red",
            "Unknown": "lime",
            "Tandem_repeat": "yellow",
            "Nonrepetitive_fraction": "dimgray"
        }

        return None

    def histogram_divergence(self, df, write_to_filepath,
                             title="Divergence from consensus sequence",
                             hue: "colname to hue by (e.g. Species)" =False):
        """
        Draw a distribution of divergence from main df or a subset of main df

        df colname with divergence values: 'perc_divg'
        df colname with colours for each bin: 'colours_class'
        `write_to_filepath`: save PNG figure to this path
        """
        if hue:
            ax = sns.histplot(data=df, x='perc_divg',
                              hue=hue, multiple="dodge",
                        palette={"Dcat": "#eb8f46",
                                 "Dtil": "#30acac"})
            plt.title(title, fontsize=12)
            plt.xlabel("% divergence")
            plt.savefig(write_to_filepath, dpi=300)
            plt.close('all')

        else:
            ax = sns.histplot(data=df, x='perc_divg')
            plt.title(title, fontsize=12)
            plt.xlabel("% divergence")
            plt.savefig(write_to_filepath, dpi=300)
            plt.close('all')

        return None

    def histogram_both_species_and_reptype_pairs(self, folder=None,
                                            filter_overlap_label=False):
        """
        """
        d1 = self.df_main
        if filter_overlap_label:
            d1 = d1.loc[d1["overlapping"]]

        for rep_class in d1["class"].unique():
            d2 = d1.loc[d1["class"] == rep_class]
            if folder:
                if not os.path.exists(folder):
                    os.makedirs(folder)
                write_to_filepath = (folder+"/matplotlib_divhist_"+
                    str(rep_class)+'.png')
            else:
                write_to_filepath = ('matplotlib_divhist_'+
                    str(rep_class)+'.png')
            self.histogram_divergence(
                d2, write_to_filepath,
                title=str(rep_class)+" divergence from consensus sequence",
                hue="Species")

            # compute subclass, order and superfam, but only
            # for both 'DNA' and 'Retrotransposon' classes
            if rep_class in ["DNA", "Retrotransposon"]:
                # by subclass
                for rep_type in d2["subclass"].unique():
                    d3 = d2.loc[d2["subclass"] == rep_type]
                    if folder:
                        if not os.path.exists(folder):
                            os.makedirs(folder)
                        write_to_filepath = (folder+'/matplotlib_divhist_'+
                            str(rep_class)+"_"+str(rep_type)+".png")
                    else:
                        write_to_filepath = ('matplotlib_divhist_'+
                            str(rep_class)+"_"+str(rep_type)+".png")
                    self.histogram_divergence(
                        d3, write_to_filepath,
                        title=str(rep_class)+"-"+str(rep_type)+
                            " divergence from consensus sequence",
                        hue="Species")

                # by order
                for rep_type in d2["order"].unique():
                    d3 = d2.loc[d2["order"] == rep_type]
                    if folder:
                        if not os.path.exists(folder):
                            os.makedirs(folder)
                        write_to_filepath = (folder+'/matplotlib_divhist_'+
                            str(rep_class)+"_"+str(rep_type)+".png")
                    else:
                        write_to_filepath = ('matplotlib_divhist_'+
                            str(rep_class)+"_"+str(rep_type)+".png")
                    self.histogram_divergence(
                        d3, write_to_filepath,
                        title=str(rep_class)+"-"+str(rep_type)+
                            " divergence from consensus sequence",
                        hue="Species")

                # by superfam
                for rep_type in d2["superfam"].unique():
                    d3 = d2.loc[d2["superfam"] == rep_type]
                    if folder:
                        if not os.path.exists(folder):
                            os.makedirs(folder)
                        write_to_filepath = (folder+'/matplotlib_divhist_'+
                            str(rep_class)+"_"+str(rep_type)+".png")
                    else:
                        write_to_filepath = ('matplotlib_divhist_'+
                            str(rep_class)+"_"+str(rep_type)+".png")
                    self.histogram_divergence(
                        d3, write_to_filepath,
                        title=str(rep_class)+"-"+str(rep_type)+
                            " divergence from consensus sequence",
                        hue="Species")

        return None

    def histogram_for_species_and_reptype_pairs(self, folder=None):
        """
        """
        df = self.df_main

        for species in df["Species"].unique():
            # subset by species
            d1 = df.loc[df["Species"] == species]
            for rep_class in d1["class"].unique():
                d2 = d1.loc[d1["class"] == rep_class]
                if folder:
                    if not os.path.exists(folder):
                        os.makedirs(folder)
                    write_to_filepath = (folder+'/matplotlib_divhist_'+
                        str(species)+'_'+str(rep_class)+'.png')
                else:
                    write_to_filepath = ('matplotlib_divhist_'+
                        str(species)+'_'+str(rep_class)+'.png')
                self.histogram_divergence(
                    d2, write_to_filepath,
                    title=str(rep_class)+" divergence from consensus sequence")

                # compute subclass, order and superfam, but only
                # for both 'DNA' and 'Retrotransposon' classes
                if rep_class in ["DNA", "Retrotransposon"]:
                    # by subclass
                    for rep_type in d2["subclass"].unique():
                        d3 = d2.loc[d2["subclass"] == rep_type]
                        if folder:
                            if not os.path.exists(folder):
                                os.makedirs(folder)
                            write_to_filepath = (folder+'/matplotlib_divhist_'+
                                str(species)+"_"+str(rep_class)+"_"+str(rep_type)+".png")
                        else:
                            write_to_filepath = ('matplotlib_divhist_'+
                                str(species)+"_"+str(rep_class)+"_"+str(rep_type)+".png")
                        self.histogram_divergence(
                            d3, write_to_filepath,
                            title=str(rep_class)+"-"+str(rep_type)+
                                " divergence from consensus sequence")

                    # by order
                    for rep_type in d2["order"].unique():
                        d3 = d2.loc[d2["order"] == rep_type]
                        if folder:
                            if not os.path.exists(folder):
                                os.makedirs(folder)
                            write_to_filepath = (folder+'/matplotlib_divhist_'+
                                str(species)+"_"+str(rep_class)+"_"+str(rep_type)+".png")
                        else:
                            write_to_filepath = ('matplotlib_divhist_'+
                                str(species)+"_"+str(rep_class)+"_"+str(rep_type)+".png")
                        self.histogram_divergence(
                            d3, write_to_filepath,
                            title=str(rep_class)+"-"+str(rep_type)+
                                " divergence from consensus sequence")

                    # by superfam
                    for rep_type in d2["superfam"].unique():
                        d3 = d2.loc[d2["superfam"] == rep_type]
                        if folder:
                            if not os.path.exists(folder):
                                os.makedirs(folder)
                            write_to_filepath = (folder+'/matplotlib_divhist_'+
                                str(species)+"_"+str(rep_class)+"_"+str(rep_type)+".png")
                        else:
                            write_to_filepath = ('matplotlib_divhist_'+
                                str(species)+"_"+str(rep_class)+"_"+str(rep_type)+".png")
                        self.histogram_divergence(
                            d3, write_to_filepath,
                            title=str(rep_class)+"-"+str(rep_type)+
                                " divergence from consensus sequence")

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
                        palette={"Dcat": "#eb8f46",
                                 "Dtil": "#30acac"},
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
        # remove repetitive fraction (sum of all repeats)
        df = df.loc[(df["class"]=="Repetitive_fraction")==False]
        # we are interested in sorting classes manually into the following order
        ordering_classes = {
            "DNA": 1, "Retrotransposon": 2, "Tandem_repeat": 3,
            "Other": 4, "Unknown": 5, "Nonrepetitive_fraction": 6}
        # we are also interested in secondarily sorting rows by 'algor_bpsum'
        # (sum of basepairs for that entry/row in the table)
        df = df.sort_values(["algor_bpsum"], ascending=True) # start by bp-sort
        df = df.sort_values(["class"], key=lambda x: x.map(ordering_classes))
        df = df.sort_values(["Species", "sequid_type"]).reset_index()
        print(df.to_markdown())

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

        bottom = np.zeros(6)#debug (change 6 to len(species+genome subset))
        for rep_class in df["class"].unique():
            d1 = df.loc[df["class"] == rep_class]
            bp_values = dict()
            for rep_superfam in d1["superfam"].unique():
                d2 = d1.loc[d1["superfam"] == rep_superfam]
                for species in df["Species"].unique():
                    d3 = d2.loc[d2["Species"] == species]
                    for gsub_key in regex_to_gensubs.keys():
                        # subset df where sequid type contains a regex
                        d4 = d3.loc[d3["sequid_type"].str.contains(
                            regex_to_gensubs[gsub_key])]
                        xaxis = str(species) + "_" + str(gsub_key)
                        # sum of basepairs for the triple-subsetted df
                        bp_values[xaxis] = d4["algor_bpsum"].sum()

                # create the segment bar
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

def histogram_stacked_absolute(df):
    """
    INPUT
    =====
    pandas.DataFrame from the function `compress_super_fams()`

    OUTPUT
    ======
    PNG file `filename` with the plotted histogram

    See more:
    https://www.pythoncharts.com/matplotlib/stacked-bar-charts-labels/
    """
    df["xval"] = df["Species"]+"_"+df["genome_subset"]
    df["reptype"] = df["class"]+"_"+df["superfam"]
    #first colour for nonrepetitive-fraction should be gray
    # the other colors we do not care as much
    needed_colours = len(df["reptype"].unique())
    colors = sns.color_palette("hls", needed_colours)
    colors[0] = "dimgrey"

    # remove edgecolor from bars
    plt.rcParams["patch.edgecolor"] = "none"
    # set the size of the plot (wider than tall)
    plt.figure(figsize=(8, 14))

    ax = sns.histplot(data=df, x='xval', hue="reptype", weights="sum",
                multiple="stack", shrink=0.9, palette=colors)
    plt.grid(axis='y')
    # move legend outside of the plotting box
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1,1))
    ax.tick_params(axis='y', labelsize=20, rotation=90)
    ax.tick_params(axis='x', labelsize=15)
    plt.title("Distribution of repetitive and nonrepetitive fraction", y=1.02,
              fontsize=20)
    # adjust margins of figure, so legend, axis, etc.
    # has enough space to be drawn
    plt.subplots_adjust(left=0.3, bottom=0.18, right=0.59, top=0.95)
    plt.ylabel("Gb", fontsize=18)
    plt.xticks(rotation=90)
    plt.xlabel("")
    plt.savefig('histo_stacked_absolute_matplotlib.png', dpi=300)
    plt.close('all')

    return None

def add_coloured_palette(df,
    colname: "string which names a categorical column"):
    """
    """
    df_return = df.copy()
    df_return["colours"]="k"
    unique_categories = df[colname].unique()
    list_cols = sns.color_palette("hls", len(unique_categories))
    for uniqcell, col in zip(unique_categories, list_cols):
        mask = df_return[colname]==uniqcell
        df_return.loc[mask, "colours"] = pd.Series([
            col for x in range(df_return.loc[mask].shape[0])],
            index=df_return.loc[mask].index)

    return df_return

def pies_relative(df: "`df_complete` TSV",
        percent_threshold: "Group rows below this percentage"=5):
    """
    """

    # label "autosomes", "allosome" and "genome" via regexes
    # of the column "sequid_type"
    regex_to_gensubs = {  # regex to genome subsetting
        "Autosomes": r"chr\d", "Allosome": r"chrX", "Genome": r".*"}
    # fill NaN values with an interpretable "NA" string
    df = df.fillna("NA")

    # PIES OF CLASSES #
    fig, axs = plt.subplots(nrows=3, ncols=len(df["Species"].unique()),
                             figsize=(18,12))
    # prepare the same colours for each given class across species
    # colname: colours; palette from seaborn.colors_palette("hls")
    d1 = add_coloured_palette(df, "class")
    # each species is a column, and each subset of genome a row
    i=0
    for species in d1["Species"].unique():
        d = d1.loc[d1["Species"]==species]
        # remove repetitive fraction sum from rows
        d = d.loc[(d["class"]=="Repetitive_fraction")==False]
        print("\n# "+species, "All classes")
        j=0
        for key_subset, regex in regex_to_gensubs.items():
            dd = d.loc[d["sequid_type"].str.contains(regex)]
            dd = dd.groupby(["class", "colours"])[
                "algor_bpsum"].agg(["sum", "count"])
            dd = dd.sort_values(["sum"], ascending=False).reset_index()
            dd["percent"] = (dd["sum"]/dd["sum"].sum())*100
            print(dd.to_markdown(index=False))
            # if there are 2 or more "percent" values
            # below the threshold, group them:
            pthresh_mask = dd["percent"] < percent_threshold
            if dd.loc[pthresh_mask].shape[0]>=2:
                new_row={
                    "class": ["Group_Others"],
                    "colours": ["grey"],
                    "sum": [dd.loc[pthresh_mask, "sum"].sum()],
                    "count": [dd.loc[pthresh_mask, "count"].sum()],
                    "percent": [dd.loc[pthresh_mask, "percent"].sum()],
                }
                # remove these rows below threshold
                dd = dd.loc[pthresh_mask==False]
                # add the new grouped row
                dd = pd.concat([dd, pd.DataFrame(new_row)])
            dd["megabases_span"]=round(dd["sum"]/10**6, 1)
            # append count of types and sum of bps to the label of the pie
            dd["class"] = (
            dd["class"]+"_"+dd["count"].astype(str)+"_"+dd["megabases_span"].astype(str)+"Mb")
            axs[j, i].pie(dd["percent"], labels=dd["class"],
                          autopct='%.2f%%', colors=dd["colours"])
            axs[j, i].set_title(species+"_"+str(key_subset))
            j+=1
        i+=1
    fig.suptitle("By class")
    plt.savefig('pies_byall_class_matplotlib.png', dpi=300)
    plt.close('all')

    # PIES OF SUBCLASSES #
    # do not iterate all unique classes, only ones of interest...
    for repeat_class in ["DNA", "Retrotransposon"]:
        fig, axs = plt.subplots(nrows=3, ncols=len(df["Species"].unique()),
                             figsize=(18,12))
        # subset df by `repeat_class` and create `colours`
        # column for each subclass
        d1 = df.loc[df["class"]==repeat_class]
        d1 = add_coloured_palette(d1, "subclass")
        i=0
        for species in d1["Species"].unique():
            # remove repetitive fraction sum from rows
            d2 = d1.loc[((d1["class"]=="Repetitive_fraction")==False) &
                      (d1["Species"]==species)]
            j=0
            print("\n# "+species, repeat_class)
            for key_subset, regex in regex_to_gensubs.items():
                d3 = d2.loc[d2["sequid_type"].str.contains(regex)]
                d3 = d3.groupby(["subclass", "colours"])[
                    "algor_bpsum"].agg(["sum", "count"])
                d3 = d3.sort_values(["sum"], ascending=False).reset_index()
                d3["percent"] = (d3["sum"]/d3["sum"].sum())*100
                print(d3.to_markdown(index=False))
                # if there are 2 or more "percent" values
                # below the threshold, group them:
                pthresh_mask = d3["percent"] < percent_threshold
                if d3.loc[pthresh_mask].shape[0]>=2:
                    new_row={
                        "subclass": ["Group_Others"],
                        "colours": ["grey"],
                        "sum": [d3.loc[pthresh_mask, "sum"].sum()],
                        "count": [d3.loc[pthresh_mask, "count"].sum()],
                        "percent": [d3.loc[pthresh_mask, "percent"].sum()],
                    }
                    # remove these rows below threshold
                    d3 = d3.loc[pthresh_mask==False]
                    # add the new grouped row
                    d3 = pd.concat([d3, pd.DataFrame(new_row)])
                d3["megabases_span"]=round(d3["sum"]/10**6, 1)
                # append count of types and sum of bps to the label of the pie
                d3["subclass"] = (
                d3["subclass"]+"_"+d3["count"].astype(str)+"_"+d3["megabases_span"].astype(str)+"Mb")
                axs[j, i].pie(d3["percent"], labels=d3["subclass"],
                              autopct='%.2f%%', colors=d3["colours"])
                axs[j, i].set_title(species+"_"+str(key_subset))
                j+=1
            i+=1
        fig.suptitle("By subclass in "+repeat_class)
        plt.savefig('pies_by'+repeat_class+'_subclass_matplotlib.png', dpi=300)
        plt.close('all')

    # PIES OF ORDER #
    # do not iterate all unique classes, only ones of interest...
    for repeat_subclass in [("DNA", "1|2|NA"), ("Retrotransposon", ".*")]:
        fig, axs = plt.subplots(nrows=3, ncols=len(df["Species"].unique()),
                             figsize=(18,12))
        # subset df by `repeat_subclass` and create `colours`
        # column for each order
        d1 = df.loc[(df["class"].str.contains(repeat_subclass[0])) &
                    (df["subclass"].str.contains(repeat_subclass[1]))]
        d1 = add_coloured_palette(d1, "order")
        i=0
        for species in d1["Species"].unique():
            # remove repetitive fraction sum from rows
            d2 = d1.loc[((d1["class"]=="Repetitive_fraction")==False) &
                      (d1["Species"]==species)]
            j=0
            print("\n# "+species, repeat_subclass)
            for key_subset, regex in regex_to_gensubs.items():
                d3 = d2.loc[d2["sequid_type"].str.contains(regex)]
                d3 = d3.groupby(["order", "colours"])[
                    "algor_bpsum"].agg(["sum", "count"])
                d3 = d3.sort_values(["sum"], ascending=False).reset_index()
                d3["percent"] = (d3["sum"]/d3["sum"].sum())*100
                print(d3.to_markdown(index=False))
                # if there are 2 or more "percent" values
                # below the threshold, group them:
                pthresh_mask = d3["percent"] < percent_threshold
                if d3.loc[pthresh_mask].shape[0]>=2:
                    new_row={
                        "order": ["Group_Others"],
                        "colours": ["grey"],
                        "sum": [d3.loc[pthresh_mask, "sum"].sum()],
                        "count": [d3.loc[pthresh_mask, "count"].sum()],
                        "percent": [d3.loc[pthresh_mask, "percent"].sum()],
                    }
                    # remove these rows below threshold
                    d3 = d3.loc[pthresh_mask==False]
                    # add the new grouped row
                    d3 = pd.concat([d3, pd.DataFrame(new_row)])
                d3["megabases_span"]=round(d3["sum"]/10**6, 1)
                # append count of types and sum of bps to the label of the pie
                d3["order"] = (
                d3["order"]+"_"+d3["count"].astype(str)+"_"+d3["megabases_span"].astype(str)+"Mb")
                axs[j, i].pie(d3["percent"], labels=d3["order"],
                              autopct='%.2f%%', colors=d3["colours"])
                axs[j, i].set_title(species+"_"+str(key_subset))
                j+=1
            i+=1
        fig.suptitle("By orders in "+repeat_subclass[0])
        plt.savefig('pies_by'+repeat_subclass[0]+'_orders_matplotlib.png', dpi=300)
        plt.close('all')

    # PIES OF SUPERFAMILIES #
    # do not iterate all unique classes, only ones of interest...
    for repeat_subclass in [("DNA", ".*", ".*"), ("Retrotransposon", ".*", ".*")]:
        fig, axs = plt.subplots(nrows=3, ncols=len(df["Species"].unique()),
                             figsize=(18,12))
        # subset df by `repeat_subclass` and create `colours`
        # column for each order
        d1 = df.loc[(df["class"].str.contains(repeat_subclass[0])) &
                    (df["subclass"].str.contains(repeat_subclass[1])) &
                    (df["order"].str.contains(repeat_subclass[2]))]
        d1 = add_coloured_palette(d1, "superfam")
        i=0
        for species in d1["Species"].unique():
            # remove repetitive fraction sum from rows
            d2 = d1.loc[((d1["class"]=="Repetitive_fraction")==False) &
                      (d1["Species"]==species)]
            j=0
            print("\n# "+species, repeat_subclass)
            for key_subset, regex in regex_to_gensubs.items():
                d3 = d2.loc[d2["sequid_type"].str.contains(regex)]
                d3 = d3.groupby(["superfam", "colours"])[
                    "algor_bpsum"].agg(["sum", "count"])
                d3 = d3.sort_values(["sum"], ascending=False).reset_index()
                d3["percent"] = (d3["sum"]/d3["sum"].sum())*100
                print(d3.to_markdown(index=False))
                # if there are 2 or more "percent" values
                # below the threshold, group them:
                pthresh_mask = d3["percent"] < percent_threshold
                if d3.loc[pthresh_mask].shape[0]>=2:
                    new_row={
                        "superfam": ["Group_Others"],
                        "colours": ["grey"],
                        "sum": [d3.loc[pthresh_mask, "sum"].sum()],
                        "count": [d3.loc[pthresh_mask, "count"].sum()],
                        "percent": [d3.loc[pthresh_mask, "percent"].sum()],
                    }
                    # remove these rows below threshold
                    d3 = d3.loc[pthresh_mask==False]
                    # add the new grouped row
                    d3 = pd.concat([d3, pd.DataFrame(new_row)])
                d3["megabases_span"]=round(d3["sum"]/10**6, 1)
                # append count of types and sum of bps to the label of the pie
                d3["superfam"] = (
                d3["superfam"]+"_"+d3["count"].astype(str)+"_"+d3["megabases_span"].astype(str)+"Mb")
                axs[j, i].pie(d3["percent"], labels=d3["superfam"],
                              autopct='%.2f%%', colors=d3["colours"])
                axs[j, i].set_title(species+"_"+str(key_subset))
                j+=1
            i+=1
        fig.suptitle("By superfams in "+repeat_subclass[0])
        plt.savefig('pies_by'+repeat_subclass[0]+'_superfams_matplotlib.png', dpi=300)
        plt.close('all')

    return None

def compress_super_fams(df: "`df_complete` TSV",
                        minimum_percent: "int [0-100]" =5):
    """
    Compress superfamilies that are smaller than 5%
    of the repetitive fraction (group as "ClassX - Other"
    whilst keeping tack of grouped superfamilies)

    Moreover, group by "Autosomes", "Allosome" and "Genome"

    minimum_percent: group superfamilies to class_Other which occupy less than
    this value of the sequid length.
    """
    pre_df_return = {
        "Species": [], "genome_subset": [], "class": [],
        "superfam": [], "algor_bpsum": []}
    regex_to_gensubs = {
        "Autosomes": r"chr\d", "Allosome": r"chrX", "Genome": r".*"}
    df["tupled_repclass"] = tuple(zip(
                    list(df["class"]),
                    list(df["subclass"]),
                    list(df["order"]),
                    list(df["superfam"]) ))

    for species in df["Species"].unique():
        d1 = df.loc[df["Species"]==species]
        for genome_subset, regex in regex_to_gensubs.items():
            d2 = d1.loc[d1["sequid_type"].str.contains(regex, regex=True)]
            df_fraction_sum = d2.loc[d2["class"].str.contains("_fraction"),
                                 "algor_bpsum"].sum()
            for repclass in d2["tupled_repclass"].unique():
                d3 = d2.loc[d2["tupled_repclass"]==repclass]
                pre_df_return["Species"].append(species)
                pre_df_return["genome_subset"].append(genome_subset)
                pre_df_return["algor_bpsum"].append(d3["algor_bpsum"].sum())
                # check whether this superfamily is greater than 5% of the
                # repetitive fraction and is not NaN value...
                pre_df_return["class"].append(repclass[0])
                if d3["algor_bpsum"].sum() > (df_fraction_sum*minimum_percent)/100:
                    if pd.isna(repclass[3])==False:
                        pre_df_return["superfam"].append(repclass[3])
                    elif pd.isna(repclass[2])==False:
                        pre_df_return["superfam"].append(repclass[2])
                    elif pd.isna(repclass[1])==False:
                        pre_df_return["superfam"].append(repclass[1])
                    else:
                        pre_df_return["superfam"].append("NA")
                else:
                    pre_df_return["superfam"].append("Other")

    r = pd.DataFrame(pre_df_return)
    r = r.groupby(["Species", "genome_subset", "class",
                "superfam"])["algor_bpsum"].agg(["sum", "count"])
    r = r.sort_values(["sum"], ascending=False).sort_values(
        ["Species", "genome_subset"]).reset_index()
    r = r.loc[(r["class"]=="Repetitive_fraction")==False]

    return r


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

        # create plots
#        plots.boxplots_df_main(folder="Boxplots")
#        plots.histogram_for_species_and_reptype_pairs(folder="Divg_hist")
#        plots.histogram_both_species_and_reptype_pairs(folder="Divg_hist_remove_overlap", filter_overlap_label=True)
#        plots.histogram_both_species_and_reptype_pairs(folder="Divg_hist_both_species")
        plots.histogram_stacked_absolute()

    else:
        print("Send either 'r' or 'p' as first argv "+
              "to specify in which mode the script should run")



#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# resquitllar.py
#
# 14 Nov 2023  <adria@molevol-OptiPlex-9020>

help_msg = """
RM.out files (--rmlist)
=======================

RepeatMasker's output files. They must be paired with a species nickname/label,
which must match the species label found in the index file. For instance:

script.py --rmfiles file1,sp1 file2,sp2

Where file1 and file2 are different files from RM output, and sp1 and sp2 are
species labels (separated by a comma).

Index files (--idxlist)
=======================

They are TSV (Tab-separated Values) files with sequence regexes paired with the
sum of their lengths. The header (first row) must contain an species nickname.
For example:

```idxfile
Species-label
Sequence-regex1\t32
Sequence-regex2\t5600
etc.
```

Another example:

```idxfile
Species       whatever
Chr1          3600
Chr2          5000
...           ...

ChrN          100
Scaffolds     $(sum_len_scaffolds)
```

Where $(sum_len_scaffolds) is the cumulative sum of scaffold lengths. See the
comments below `class Repeats` for more information on regexes and formatting.

Formatting the RM.out files
===========================

The output obtained from RepeatMasker must be modified
Some modifications must be applied to the output obtained from RepeatMasker.

Some rows contain an asterisk (*) as a 16th field. Some other rows do end at the
15th field. To make sure that all rows are of the same length, one can run:

```bash
awk 'BEGIN{OFS="\t"} ; {if ( $16 ) { $16="True" } else { $16="False" } } 1' RM.out
```

Moreover, make sure to correctly label minor scaffolds and chromosomes in column
5. For instance, to substitute "Scaffold_1" by "DtilchrX" in the sequid field:

```bash
awk '$5=="Scaffold_1" {$5="DtilchrX"}' RM.out
```

Lastly, it is possible for the field "position (left)" to contain asterisks.
These obstruct the loading of the file because the column is otherwise of the
integer dtype/category. One solution is to substitute them by zeroes.
"""

import sys, math
# Reading input files' size and creating folders of PNG plots
import os

# Handling input as DataFrames
import pandas as pd, numpy as np
# Plotting figures
import seaborn as sns, matplotlib.pyplot as plt
# Measure the time it takes to run this script
from datetime import datetime
import locale
locale.setlocale(locale.LC_TIME, '') # sets locale to the one used by user

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def prepare_plotting_folder(file, folder=None):
    """
    Makes a new directory if the folder path does not already exist.
    Concatenates the filename (which should include the extension) with the
    folder path.

    Input
    =====

    + file: The filename where the newly created plot should be written to.

    + folder [optional]: The folder name where the newly created plot should be
      written to. If it is not specified the function will return the
      filename, only.

    Output
    ======

    Returns a path where the plot can be written to, as a string. Creates the
    folder if it did not exist.
    """
    if folder:
        # Create the folder if it did not exist previously
        if not os.path.exists(folder):
            os.makedirs(folder)
        # Add an ending forward slash to the directory path if necessary:
        if folder[-1] != "/":
            folder += "/"
        # Return the img path, including the directories:
        return folder + file
    else:
        return file

class Printing:
    """
    Print a given message with a formatting of interest (eg. as a warning,
    status, error...). Includes date and time of printing.

    Use like so:

    Printing("I love fries").status()
    Printing("Don't eat too many fries").warning()
    """

    def __init__(self, message):
        """
        Add the date and time of day to all the messages.
        """
        d = datetime.now().strftime("%d %b, %H:%M")
        self.m = f"{d}: {message}"

    def status(self):
        """
        Print in green and with preceding "STATUS" string
        """
        message = f"STATUS:{self.m}"
        # Change the colour to green (through ANSI colours in the terminal)
        message = "\033[0;32m" + str(message) + "\033[0;0m"
        print(message)
        return message

    def warning(self):
        """
        Print in yellow and with preceding "WARNING" string
        """
        message = f"WARNING:{self.m}"
        # Change the colour to yellow (through ANSI colours in the terminal)
        message = "\033[0;33m" + str(message) + "\033[0;0m"
        print(message)
        return message

    def error(self):
        """
        Print in red and with preceding "ERROR" string
        """
        message = f"ERROR:{self.m}"
        # Change the colour to red (through ANSI colours in the terminal)
        message = "\033[0;31m" + str(message) + "\033[0;0m"
        print(message)
        return message

class Repeats:
    """
    The class `Repeats` reads RepeatMasker's output, standardises the repeat
    type classification and returns summary tables of repeat content.

    RepeatMasker's documentation:
    <https://www.repeatmasker.org/webrepeatmaskerhelp.html>

    Input
    =====

    + files_dict: A dictionary with "species" as keys and "filepaths" as
      values. The "species" of this dictionary and "seqsizes_dict" should
      match. The filepath leads to RepeatMasker's output.

    + seqsizes_dict: A dictionary pairing sequence ID regexes (for
      instance, *chr* to match all chromosomes) with their cumulative length,
      in base pairs. Keys are species and values another dictionary. See
      below for reference:

      dict("species1": dict("chr1": length.int,
                            "chr2": length.int,
                            ...
                            "minor_scaffolds": sum_len_scaffolds.int
                            ),
           "species2": dict(...),
          )

      Sequid strings (ie. "chr1", "chr2", etc.) should be a regex which can
      match multiple sequences.
      Use '|' for OR matching: `string1|string2|string3`, `Scaff|ctg`...
      Use '&' for AND matching: `Sca&ffolds`...

    Output
    ======

    TBD

    Methods
    =======

    TBD
    """
    def __init__(self, files_dict, seqsizes_dict):
        # Store the dictionary with sequence sizes for use down the line.
        self.seqsizes_dict = seqsizes_dict
        # Compute whole genome sizes (sum of all sequences).
        gensizes_dict = dict()
        for species in seqsizes_dict.keys():
            gensizes_dict[species] = 0
            for sequid in seqsizes_dict[species].keys():
                gensizes_dict[species] += seqsizes_dict[species][sequid]
        self.gensizes_dict = gensizes_dict

        # Ensure the filepath to the provided files exist.
        for file in files_dict.values():
            try:
                with open(file) as f:
                    pass
            except:
                Printing("Cannot read the provided path to `"+
                         str(file)+"`").error()
                return None

        # Read input tabulated files. Create dfs with pandas.
        dfs_catalog = dict()
        Printing("Checking file sizes:").status()
        for species, filepath in files_dict.items():
            # Print file size before reading.
            file_size = os.path.getsize(file)
            print(f" * {filepath}:  {file_size} bytes")
            # Load as pd.DataFrame
            dfs_catalog[species] = pd.read_table(filepath,
                header=None,
                skiprows=3,     # skip the first three rows
                sep="\s+",      # fields separated by 'space' regex
                names=[         # use these strings as column names
                    # For an explanation of RepeatMasker's fields, see its
                    # manual.
                    "score_SW", "perc_divg", "perc_del", "perc_ins", "sequid",
                    "begin", "end", "left", "orient", "name",
                    "default_repclass", "begin_match", "end_match",
                    "left_match", "ID", "overlapping"],
                                                 )
            # Compute element length. Compared to consensus' length would be a
            # proxy for divergence. Disabled because we already have "perc_divg".
            ## dfs_catalog[species]["replen"] = abs(df["end"] - df["begin"]) +1

        # Count rows matching each sequid regex in seqsizes_dict. Make sure all
        # rows match with a single regex; not more, not less.
        Printing("Revising sequid regex matches between provided index file "+
                 "and RepeatMasker's output...").status()
        for species in seqsizes_dict.keys():
            df_species = dfs_catalog[species]
            df_species_totnrows = df_species.shape[0]
            df_species_nuniqseq = len(list(df_species["sequid"].unique()))
            df_species_matchrows = 0
            df_species_matchuniqseq = 0
            for sequid_regex in seqsizes_dict[species].keys():
                # Filter dataframe by sequids matching "sequid_regex"
                df_species_and_regex = df_species.loc[
                    df_species["sequid"].str.contains(sequid_regex,
                                    case=False)] # case-insensitive.
                print("--", sequid_regex, "--")
                # Print the number of rows matching and the number of unique
                # sequids that did match.
                print(" + N.rows:", df_species_and_regex.shape[0])
                df_species_matchrows += df_species_and_regex.shape[0]
                print(" + Regex matching count:",
                      len(list(df_species_and_regex["sequid"].unique())),
                      "unique sequids")
                df_species_matchuniqseq += (
                    len(list(df_species_and_regex["sequid"].unique())))
#                print(" + Regex matching list:",
#                      list(df_species_and_regex["sequid"].unique()))
##                # add new col to keep track of these 'sequid_types'
##                df_concat.loc[
##                    (df_concat["Species"] == species) &
##                    (df_concat["sequid"].str.contains(sequid_regex)),
##                    "sequid_type"] = sequid_regex
            print("-----------------------------")
            print(species, "summary: matched", df_species_matchrows,
                  "out of", df_species_totnrows, "rows (",
                  round((df_species_matchrows/df_species_totnrows) *100, ndigits=3),
                  "%)")
            print(" " * len(species + "summary: "), "found",
                  df_species_matchuniqseq, "out of", df_species_nuniqseq,
                  "unique sequids (",
                  round((df_species_matchuniqseq/df_species_nuniqseq) *100,
                        ndigits=3), "%)\n")

        # Reclassify repeat types into 4 new columns.
        [self.reclassify_default_reptypes(df) for df in dfs_catalog.values()]

        # Clean the dataframes (drop columns, change dtypes of columns).
        for df in dfs_catalog.values():
            # Drop some unnecessary columns:
            df.drop(columns=["left", "begin_match", "end_match",
                             "left_match"], inplace=True)
            # Assess the dtypes of all columns and try to downcast numbers
            # (improve memory usage; see dataframe.info(memory_usage="deep")
            float_cols = df.select_dtypes("float").columns
            df[float_cols] = df[float_cols].apply(pd.to_numeric,
                                                  downcast="float")
            int_cols = df.select_dtypes("integer").columns
            df[int_cols] = df[int_cols].apply(pd.to_numeric,
                                              downcast="integer")
            obj_cols = df.select_dtypes("object").columns
            df[obj_cols] = df[obj_cols].astype("category")

        # Check whether the reclassification correctly labelled all default
        # repeat types.
        for sp, df in dfs_catalog.items():
            Printing("Checking the reclassification of "+str(sp)).status()
            self.check_reclassification(df)

        # debug
        [df.info(memory_usage="deep") for df in dfs_catalog.values()]
        [print(df) for df in dfs_catalog.values()]

    def reclassify_default_reptypes(self, df):
        """
        Very long function of many paragraphs. Each paragraph reclassified a
        default repeat type into four standardised clades, which can be easily
        read and analysed.

        Some repeat types are reclassified if the "default_repclass" field
        exactly matches a single string (Y==X), while other repeat types are
        reclassified if the "default_repclass" matches/contains a regex-like
        string (Y.str.contains(X)).

        Input
        =====

        + df: A dataframe from the __init__ function of self (in `dfs_catalog`).

        Output
        ======

        A modified df, with more columns, standardising the repeat types.
        """
        # Unclassified

        df.loc[
            df["default_repclass"]=="__unknown",
            ["class", "subclass", "order", "superfam"]
        ] = ("Unclass.", "NA", "NA", "NA")

        # Enzymatic RNAs.

        df.loc[
            df["default_repclass"]=="tRNA",
            ["class", "subclass", "order", "superfam"]
        ] = ("Other", "tRNA", "NA", "NA")

        df.loc[
            df["default_repclass"]=="rRNA",
            ["class", "subclass", "order", "superfam"]
        ] = ("Other", "rRNA", "NA", "NA")

        df.loc[
            df["default_repclass"]=="srpRNA",
            ["class", "subclass", "order", "superfam"]
        ] = ("Other", "srpRNA", "NA", "NA")

        df.loc[
            df["default_repclass"]=="snRNA",
            ["class", "subclass", "order", "superfam"]
        ] = ("Other", "snRNA", "NA", "NA")

        df.loc[
            df["default_repclass"]=="scRNA",
            ["class", "subclass", "order", "superfam"]
        ] = ("Other", "scRNA", "NA", "NA")

        # Generic repetitive elements.

        df.loc[
            df["default_repclass"]=="Simple_repeat",
            ["class", "subclass", "order", "superfam"]
        ] = ("Tandem_repeat", "Simple_repeat", "NA", "NA")

        df.loc[
            df["default_repclass"]=="Satellite",
            ["class", "subclass", "order", "superfam"]
        ] = ("Tandem_repeat", "Satellite", "NA", "NA")

        df.loc[
            df["default_repclass"]=="Satellite/W-chromosome",
            ["class", "subclass", "order", "superfam"]
        ] = ("Tandem_repeat", "Satellite", "W-chromosome", "NA")

        df.loc[
            df["default_repclass"]=="Satellite/Y-chromosome",
            ["class", "subclass", "order", "superfam"]
        ] = ("Tandem_repeat", "Satellite", "Y-chromosome", "NA")

        df.loc[
            df["default_repclass"].str.contains("Satellite/acro"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Tandem_repeat", "Satellite", "acromeric", "NA")

        df.loc[
            df["default_repclass"].str.contains("Satellite/centr"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Tandem_repeat", "Satellite", "centromeric", "NA")

        df.loc[
            df["default_repclass"].str.contains("Satellite/macro"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Tandem_repeat", "Satellite", "macro", "NA")

        df.loc[
            df["default_repclass"].str.contains("Low_complexity"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Other", "Low_complexity", "NA", "NA")

        df.loc[
            df["default_repclass"].str.contains("Other"),
            ["class", "subclass", "order", "superfam"]
        ] = ("NA", "NA", "NA", "NA")

        # Generic retrotransposons.

        df.loc[
            df["default_repclass"].str.contains("Retroposon") |
            df["default_repclass"].str.contains("__ClassI"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "NA", "NA")

        df.loc[
            df["default_repclass"].str.contains("LINE"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "LINE", "NA")

        df.loc[
            (df["default_repclass"].str.contains("_LTR")) |
            (df["default_repclass"] == "LTR"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "LTR", "NA")

        df.loc[
            df["default_repclass"].str.contains("SINE"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "SINE", "NA")

        df.loc[
            df["default_repclass"].str.contains("Retroposon") &
            (df["default_repclass"].str.contains("SVA") |
            df["default_repclass"].str.contains("L1-dep")),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "NA", "L1-dependent")

        df.loc[
            df["default_repclass"].str.contains("Retroposon/RTE-derived"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "NA", "RTE-derived")

        # Superfamilies in LINE order.

        df.loc[
            (df["default_repclass"].str.contains("LINE") &
            df["default_repclass"].str.contains("Penelope")) |
            (df["default_repclass"].str.contains("ClassI") &
            df["default_repclass"].str.contains("PLE")),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "PLE", "Penelope")

        df.loc[
            df["default_repclass"].str.contains("LINE/CR1") |
            df["default_repclass"].str.contains("LINE/L2") |
            df["default_repclass"].str.contains("LINE/Rex-Babar"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "LINE", "CR1")

        df.loc[
            df["default_repclass"].str.contains("LINE/CRE"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "LINE", "CRE")

        df.loc[
            df["default_repclass"].str.contains("LINE/Deceiver"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "LINE", "Deceiver")

        df.loc[
            df["default_repclass"].str.contains("LINE/Dong-R4"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "LINE", "Dong-R4")

        df.loc[
            df["default_repclass"].str.contains("LINE/I") |
            df["default_repclass"].str.contains("nLTR_LINE_I"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "LINE", "I")

        df.loc[
            df["default_repclass"].str.contains("LINE/I-Jockey"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "LINE", "Jockey")

        df.loc[
            df["default_repclass"].str.contains("LINE/L1") |
            df["default_repclass"].str.contains("nLTR_LINE_L1"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "LINE", "L1")

        df.loc[
            df["default_repclass"].str.contains("LINE/Proto1"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "LINE", "Proto1")

        df.loc[
            df["default_repclass"].str.contains("LINE/Proto2"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "LINE", "Proto2")

        df.loc[
            df["default_repclass"].str.contains("LINE/R1") |
            df["default_repclass"].str.contains("LINE/Tad1"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "LINE", "R1")

        df.loc[
            df["default_repclass"].str.contains("LINE/R2") |
            df["default_repclass"].str.contains("nLTR_LINE_R2"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "LINE", "R2")

        df.loc[
            df["default_repclass"].str.contains("LINE/RTE") |
            df["default_repclass"].str.contains("nLTR_LINE_RTE"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "LINE", "RTE")

        # Superfamilies in LTR order

        df.loc[
            df["default_repclass"].str.contains("LTR") &
            df["default_repclass"].str.contains("DIRS"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "DIRS", "DIRS")

        df.loc[
            df["default_repclass"].str.contains("LTR") &
            df["default_repclass"].str.contains("Ngaro"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "DIRS", "Ngaro")

        df.loc[
            df["default_repclass"].str.contains("LTR/Viper"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "DIRS", "Viper")

        df.loc[
            df["default_repclass"].str.contains("LTR/Pao") |
            df["default_repclass"].str.contains("LTR_BEL"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "LTR", "Bel-Pao")

        df.loc[
            df["default_repclass"].str.contains("LTR") &
            df["default_repclass"].str.contains("Copia"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "LTR", "Copia")

        df.loc[
            df["default_repclass"].str.contains("LTR") &
            df["default_repclass"].str.contains("ERV"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "LTR", "ERV")

        df.loc[
            df["default_repclass"].str.contains("LTR") &
            df["default_repclass"].str.contains("Gypsy"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "LTR", "Gypsy")

        df.loc[
            df["default_repclass"].str.contains("LTR") &
            df["default_repclass"].str.contains("ERV"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "LTR", "ERV")

        df.loc[
            df["default_repclass"].str.contains("LTR") &
            df["default_repclass"].str.contains("Caulimovirus"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "LTR", "Pararetrovirus")

        # Superfamilies in SINE order


        df.loc[
            df["default_repclass"].str.contains("SINE") &
            df["default_repclass"].str.contains("5S"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "SINE", "5S")

        df.loc[
            df["default_repclass"].str.contains("SINE") &
            (df["default_repclass"].str.contains("7SL") |
            df["default_repclass"].str.contains("ALU", case=False) |
            df["default_repclass"].str.contains("B2") |
            df["default_repclass"].str.contains("B4")),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "SINE", "7SL")

        df.loc[
            df["default_repclass"].str.contains("Retroposon/L1-derived"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "SINE", "L1-derived")

        df.loc[
            df["default_repclass"].str.contains("Retroposon/L2-derived"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "SINE", "L2-derived")

        df.loc[
            df["default_repclass"].str.contains("SINE") &
            df["default_repclass"].str.contains("tRNA"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "SINE", "tRNA")

        df.loc[
            df["default_repclass"].str.contains("Retroposon/R4-derived"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "SINE", "R4-derived")

        # Generic DNA class' repeats.

        df.loc[
            (df["default_repclass"] == "DNA") |
            (df["default_repclass"] == "__ClassII"),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "NA", "NA", "NA")

        # DNA, subclass I.

        df.loc[
            df["default_repclass"].str.contains("DNA") &
            df["default_repclass"].str.contains("Crypton"),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "Crypton", "Crypton")

        df.loc[
            df["default_repclass"].str.contains("DNA") &
            df["default_repclass"].str.contains("Zator"),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "NA", "Zator")

        df.loc[
            df["default_repclass"].str.contains("ClassII_MITE"),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "NA", "MITE")

        df.loc[
            df["default_repclass"].str.contains("ClassII_nMITE"),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "NA", "nMITE")

        df.loc[
            df["default_repclass"].str.contains("DNA") &
            df["default_repclass"].str.contains("Academ"),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "TIR", "Academ")

        df.loc[
            (df["default_repclass"].str.contains("DNA") &
            df["default_repclass"].str.contains("CMC")) |
            df["default_repclass"].str.contains("ClassII_DNA_CACTA"),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "TIR", "CACTA")

        df.loc[
            (df["default_repclass"].str.contains("DNA") &
            df["default_repclass"].str.contains("Dada")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "TIR", "Dada")

        df.loc[
            (df["default_repclass"].str.contains("DNA") &
            df["default_repclass"].str.contains("Kolobok")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "TIR", "Kolobok")

        df.loc[
            (df["default_repclass"].str.contains("DNA") &
            df["default_repclass"].str.contains("Ginger")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "TIR", "Ginger")

        df.loc[
            (df["default_repclass"].str.contains("DNA") &
            df["default_repclass"].str.contains("MULE")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "TIR", "MULE")

        df.loc[
            (df["default_repclass"].str.contains("DNA") &
            df["default_repclass"].str.contains("Merlin")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "TIR", "Merlin")

        df.loc[
            (df["default_repclass"].str.contains("DNA/P")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "TIR", "P")

        df.loc[
            (df["default_repclass"].str.contains("DNA/PIF")) |
            (df["default_repclass"].str.contains("ClassII_DNA_Harbinger")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "TIR", "PIF-Harbinger")

        df.loc[
            (df["default_repclass"].str.contains("PiggyBac")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "TIR", "PiggyBac")

        df.loc[
            (df["default_repclass"].str.contains("DNA/Sola")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "TIR", "Sola")

        df.loc[
            (df["default_repclass"].str.contains("DNA/Sola")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "TIR", "Sola")

        df.loc[
            (df["default_repclass"].str.contains("DNA/TcMar") |
            df["default_repclass"].str.contains("DNA_TcMar")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "TIR", "TcMar")

        df.loc[
            (df["default_repclass"].str.contains("DNA") &
            df["default_repclass"].str.contains("hAT")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "TIR", "hAT")

        df.loc[
            (df["default_repclass"].str.contains("DNA_Mutator")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "TIR", "Mutator")

        df.loc[
            df["default_repclass"].str.contains("IS3EU"),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "NA", "IS3EU")

        df.loc[
            df["default_repclass"].str.contains("Novosib"),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "NA", "Novosib")

        df.loc[
            df["default_repclass"].str.contains("Zisupton"),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "NA", "Zisupton")

        # DNA, subclass II

        df.loc[
            (df["default_repclass"].str.contains("Helitron")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "2", "Helitron", "Helitron")

        df.loc[
            df["default_repclass"].str.contains("Maverick"),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "2", "Maverick", "Maverick")

        # MITE or nMITE?

        df.loc[
            df["default_repclass"].str.contains("_MITE"),
            ["mite"]
        ] = (True)

        df.loc[
            df["default_repclass"].str.contains("_nMITE"),
            ["mite"]
        ] = (False)

        return df

    def remove_overlapping(self, intervals):
        """
        Parse the intervals (begin, end) of RE (repetitive elements) of a given
        sequid, computing the cumulative sum whilst removing overlapping base
        pairs. Interval length is computed as `(end-begin)+1`. Intervals with
        the best score are prioritised over repeats with lower score.

        Input
        =====

        + intervals: A list of intervals (which in turn are sublists composed of
          `begin` and `end` coordinates). The pairs of coordinates must be
          integers. Score for each interval must be provided as the third item
          in the sublist. Include a tuple with the classification of the RE.

        intervals = [ [begin_1, end_1, score_1, classification_1], ...,
        [begin_N, end_N, score_N, classification_N] ]

        classification = (class_str, subclass_str, order_str, superfam_str)

        Output
        ======

        The cumulative sum of interval lengths, each computes as
        `(end-begin)+1`. The function removes overlapping segments, counting
        each nucleotide a single time.
        """
        # Init variables.
        # Create a list with all unique repeat type tuples.
        unique_class_tuples = list(set(tuple(x[3]) for x in intervals))
        unique_class_tuples.sort()
        # Convert into pd.MultiIndex()
        idx = pd.MultiIndex.from_tuples(unique_class_tuples, names=[
                'class', 'subclass', 'order', 'superfam'])
        # Create a df which will be returned as answer.
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
        # Sort point list by position (by first value in sublists)
        point_list.sort()
        # Initialise other algorithm variables...
        currentBest = -1
        open_intervals = []

        # Start iterating across the point_list, discounting overlapping base
        # pairs.
        for i in range(0, len(point_list)):
            # If the loop landed on a left point (0) opens an interval
            if point_list[i][1] == 0:

                # If there is no other interval opened:
                if currentBest == -1:
                    # Enters interval 'i'
                    currentBest = point_list[i][2]
                    currentBegin = int(point_list[i][0])
                    currentScore = int(intervals[currentBest][2])

                # Else, there already was an open interval:
                # (and it is of lower score than currentBest)
                elif currentScore > intervals[point_list[i][2]][2]:
                    # Do not close currentBest, because the next interval (repeat)
                    # is of lower score; however, append to open intervals list
                    iden = point_list[i][2]
                    open_intervals.append(
                        [int(intervals[iden][2]),     # score
                         iden])                       # ID
                    open_intervals.sort(reverse=True) # sort by score

                # Else, currentScore should be higher
                else:
                    # Compute length up until now (begin2 -begin1)
                    interval_length = point_list[i][0] -currentBegin
                    # Add this length to the resulting pandas df
                    reptype = intervals[currentBest][3]
                    df_resulting.loc[reptype, "algor_bpsum"] += interval_length
                    # Append index of past Best to open_intervals
                    iden = currentBest
                    open_intervals.append(
                        [int(intervals[iden][2]),     # score
                         iden])                       # ID
                    open_intervals.sort(reverse=True) # sort by score
                    # And change currentBest
                    currentBest = point_list[i][2]
                    currentBegin = point_list[i][0]
                    currentScore = intervals[currentBest][2]

            # If the loop landed on a right point (1) closes an interval
            # Moreover, it has to be the right point of the currently open
            # interval 'i'.
            elif point_list[i][2] == currentBest:
                # DEBUG
    ##            print('point_list[i][2] (', point_list[i], ') == currentBest')
                # Compute length up until now (end -begin +1)
                interval_length = point_list[i][0] -currentBegin +1
                # Add this length to the resulting pandas df
                reptype = intervals[currentBest][3]
                df_resulting.loc[reptype, "algor_bpsum"] += interval_length
                # Close this interval and open the next best one in open_intervals
                if len(open_intervals) == 0:
                    currentBest = -1
                else:
                    # Remove second best from open_intervals and use it as
                    # currentBest
                    c = open_intervals.pop(0)
                    currentBest = c[1]
                    currentScore = c[0]
                    currentBegin = point_list[i][0] +1

            # Otherwise, it is closing an interval listed in `open_intervals`
            else:
                # Remove the interval from `open_intervals`
                iden = point_list[i][2]
                open_intervals.remove(
                    [int(intervals[iden][2]), # score
                    iden])                    # ID

        return df_resulting

    def check_reclassification(self, df):
        """
        Checks how the default repeat types have been reclassified. Returns a
        dataframe pairing every unique default type with their assigned class,
        subclass, order and superfamily.

        Input
        =====

        + df: A dataframe generated in the __init__ function by reading RM
          output.

        Output
        ======

        Prints and returns a pandas dataframe in which default types are paired
        with their reassigned classification.
        """
        # Drop rows with non-unique cell values.
        df_unique_classes = df.drop_duplicates(
            subset=["default_repclass", "class",
                    "subclass", "order", "superfam"])
        # Sort the dataframe by the newly assigned types.
        df_unique_classes = df_unique_classes.sort_values(
            by=["class", "subclass", "order", "superfam"])
        # Reset the dataframe index
        df_unique_classes = df_unique_classes.reset_index()
        # We are only interested in columns with the repeat types;
        # drop the rest of columns
        df_unique_classes = df_unique_classes[["class", "subclass", "order",
                                               "superfam", "mite",
                                               "default_repclass"]]
        # Print the dataframe in markdown format (without writing to a file)
        print(df_unique_classes.to_markdown(None))

        return df_unique_classes


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Si es crida com a script:
if __name__ == '__main__':

    ### PARSE ARGUMENTS ###
    import argparse
    parser = argparse.ArgumentParser(description=help_msg,
        # make sure the 'help_msg' is not automatically
        # wrapped at 80 characters (manually assign newlines).
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--rmlist", nargs="+", required=True,
                        help="A list of RepeatMasker's output files.")
    parser.add_argument("--idxlist", nargs="+", required=True,
                        help="A list of index files, as specified in --help")
    # IMPROVEMENT?
    # Add argument that creates desired output format (this table, that figure)

    # file-name: positional arg.
#    parser.add_argument('filename', type=str, help='Path to ... file-name')
    # choices argument
#    parser.add_argument('-o', '--operacio', 
#            type=str, choices=['suma', 'resta', 'multiplicacio'],
#            default='suma', required=False,
#            help='')

    args = parser.parse_args()
    # call a value: args.operacio or args.filename.

    # Create the `files_dict` dictionary required by the class `Repeats`.
    files_dict = dict()
    for pair in args.rmlist:
        values = pair.split(",")
        if len(values) != 2:
            Printing("One of the specified pairs in --rmlist was not of "+
                     "length equal to two").error()
            Printing("Use --help for an explanation on --rmlist "+
                     "parameter requirements").error()
            sys.exit()
        else:
            file, species = values
            files_dict[species] = file
    # Create the `seqsizes_dict` dictionary required by the class `Repeats`.
    seqsizes_dict = dict()
    for index in args.idxlist:
        df_indexfile = pd.read_table(index, index_col=0, sep="\s+")
        seqsizes_dict[df_indexfile.index.name] = (
            df_indexfile.iloc[:,0].to_dict())

    repeats = Repeats(files_dict, seqsizes_dict=seqsizes_dict)


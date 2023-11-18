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

Species-label
Sequence-regex1\t32
Sequence-regex2\t56
etc.

See the comments below `class Repeats` for more information.
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
                Printing("Cannot read the provided path.").error()
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


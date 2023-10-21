#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# 3era o 4rta versió
#
# 20 d’oct. 2023  <adria@molevol-OptiPlex-9020>

help_msg = """
The purpose of this Python script is to parse (read) a PAF file and create
compelling plots and tables from it.

The script defines a class which creates a more manageable Python3 object from a
PAF mapping file. It stores information in a pandas dataframe.

It also contains a function that subsets a given PAF file. It is useful to
create datasets for testing purposes (before proceeding with the main analyses).
"""

import sys
import random

# Dataframes manipulation, akin to the R project
import pandas as pd
import numpy as np
# Sorting strings "humanely"
from natsort import index_natsorted
# Creating plots from `pandas` dataframes.
import matplotlib.pyplot as plt
import seaborn as sns

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def shuffle_big_into_small_file(path, outpath, nlines=2500):
    """
    This function requires the `random` module to randomly subset N lines from
    the given file.

    Create a sample with testing purposes from randomly selected lines found
    within a bigger data file (PAF mapping file). For instance, see
    <https://stackoverflow.com/q/60623645/21787735>.
    This function does not output an exact number of lines; the specified
    `nlines` will only be approximated in the output file.

    Input
    =====

    + path: A path to the PAF mapping file.

    + outpath: A path to write the output PAF file to.

    + nlines [optional]: The number of lines that will be included in the output
      testing PAF file.

    Output
    ======

    Writes to `outpath` a subset of the given PAF file. It is comprised of
    `nlines` number of lines.
    Returns `outpath`, which can be used as an input to `Mapping` class.
    """
    # Compute total number of lines in the input file.
    with open(path, 'r') as file:
        for count, line in enumerate(file):
            pass
        tlines = count
    # Open once again the file to reset pointer position to the beginning.
    # Moreover, open the output file in 'a' mode (does not delete
    # preexisting files, appends to the end).
    with open(path, 'r') as in_file, open(outpath, 'a') as out_file:
        # Compute a chance of printing a line in order to approximately obtain
        # `nlines` in the output file.
        chance = nlines / tlines
        [out_file.write(line) for count, line in enumerate(in_file)
            if random.random() < chance]

    return outpath

class Mapping:
    """
    Create a Python object from a PAF file.

    Input
    =====

    + path_to_paf: Path to a PAF file in order to read it and create a pandas
      DataFrame from it.

    Output
    ======

    Returns a Mapping() object, hopefully more manageable than the raw data.
    Small documentation of the columns in the dataframe can be found in minimap2
    manual:
    <https://lh3.github.io/minimap2/minimap2.html#10>

    In summary:
    "no_match" + "match" == "ali_len" == "cigM" + "cigI" + "cigD"
    """
    def __init__(self, path_to_paf):
        """
        """
        # Make sure that the provided path points to an existing file
        try:
            with open(path_to_paf):
                pass
        except:
            print("WARNING: The provided path does not point to an existing",
                  "file.")
            return None
        # As stated in the manual page, "PAF may optionally have additional
        # fields ...", i.e. an irregular number of columns (ragged TSV).
        # Try to deduce columns.
        # First, read the columns that are always present (until the 12th)
        colnames = ["Qname", "Qlen", "Qstart", "Qend", "strand",
                   "Tname", "Tlen", "Tstart", "Tend",
                   "matches",  # Number of matching bases
                   "ali_len",  # Sum of matches, mismatches and gaps
                   "mapQ"]     # Mapping quality
        # Create pandas DataFrame reading PAF file in `path_to_paf`
        print("STATUS: Reading the 12 first standard columns")
        df = pd.read_table(path_to_paf, header=None,
            # The regex "\s+" stands for "one or more blank spaces" (includes "\t")
            sep="\s+",
            # Load the appropiate (always present) columns and no more.
            usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
            # Use the following column names (from first to last column).
            names=colnames)

        # Create lists with the optional fields (which may not be included in
        # all rows). See minimap2 reference manual for an explanation of each
        # "tag"/column:  <https://lh3.github.io/minimap2/minimap2.html#10>
        # Keep in mind that the keys or "tags" have been manually sorted by
        # importance (so in the dataframe, important columns appear first).
        tag_to_columns = {"tp": [], "NM": [], "nn": [], "dv": [], "de": [],
                          "cg": [], "cs": [], "SA": [], "ts": [],
                          "cm": [], "s1": [], "s2": [], "MD": [], "AS": [],
                          "ms": [], "rl": [], "zd": [], }
        # Scroll through the lines in the file and populate this list:
        print("STATUS: Reading the optional columns (further than 12th)")
        with open(path_to_paf) as file:
            for count, line in enumerate(file):
                # Get a list with all the fields.
                line = line.strip("\n").split("\t")[12:]
                while line:
                    field = line.pop(0)
                    tag = field.split(":")[0]
                    value = field[5:]
                    tag_to_columns[tag] += [value]
                # Fill with "None" fields that were not found in this line:
                for key, val in tag_to_columns.items():
                    if len(val) < count+1:
                        tag_to_columns[key] += [None]
        # Add the additional/optional columns to the dataframe:
        tag_to_colnames= {"tp": "type_aln", "cm": "num_minimizers",
                          "s1": "chaining_score", "s2": "second_chain_score",
                          "NM": "mismatches_and_gaps",
                          "MD": "regenerate_refseq", "AS": "DP_ali_score",
                          "SA": "list_supp_ali", "ms": "DP_max_score",
                          "nn": "ambiguous_bases", "ts": "transcript_strand",
                          "cg": "cigar_string", "cs": "diff_string",
                          "dv": "divergence", "de": "gap_compr_diverg",
                          "rl": "length_repetitive_seeds",
                          "zd": "zd_unknown",}
        for tag, values_list in tag_to_columns.items():
            df[tag_to_colnames[tag]] = values_list
        # Change these optional columns to numerical type, if it is pertinent:
        df["num_minimizers"] = pd.to_numeric(df["num_minimizers"])
        df["chaining_score"] = pd.to_numeric(df["chaining_score"])
        df["second_chain_score"] = pd.to_numeric(df["second_chain_score"])
        df["mismatches_and_gaps"] = pd.to_numeric(df["mismatches_and_gaps"])
        df["DP_ali_score"] = pd.to_numeric(df["DP_ali_score"])
        df["DP_max_score"] = pd.to_numeric(df["DP_max_score"])
        df["ambiguous_bases"] = pd.to_numeric(df["ambiguous_bases"])
        df["divergence"] = pd.to_numeric(df["divergence"])
        df["gap_compr_diverg"] = pd.to_numeric(df["gap_compr_diverg"])
        df["length_repetitive_seeds"] = pd.to_numeric(df["length_repetitive_seeds"])
        df["zd_unknown"] = pd.to_numeric(df["zd_unknown"])
        # Drop columns with zero non-null entries (empty columns)
        df.dropna(axis="columns", how="all", inplace=True)

        print(df) # debug
        print(df.info())

##        # Read the first line of the PAF file
##        with open(path_to_paf) as file:
##            line = file.readline().split("\t")
##        # The 12 first fields are always present; start after the 12th field
##        fieldnames = [field.split(':')[0] for field in line[12:]]
##        # Add remaining additional colnames to be included in header
##        for tag in fieldnames:
##            colnames += [tag_to_colnames[tag]]
##        # Remove the tag (first five characters in columns from 12th to last).
##        for tag in fieldnames:
##            df[tag_to_colnames[tag]] = df[tag_to_colnames[tag]][5:]

        # Use an index file to adequately name the scaffolds? And their size?
        # Size is already given by the PAF file; but patterns of chr / scaff
        # might be important.

        # Detect groups of scaffolds / chromosomes with pattern matching.


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Si es crida com a script:
if __name__ == '__main__':

    ### PARSE ARGUMENTS ###
    import argparse
    parser = argparse.ArgumentParser(description=help_msg,
        # make sure the 'help_msg' is not automatically
        # wrapped at 80 characters (manually assign newlines).
        formatter_class=argparse.RawTextHelpFormatter)
    # file-name: positional arg.
#    parser.add_argument('filename', type=str, help='Path to ... file-name')
    # integer argument
#    parser.add_argument('-a', '--numero_a', type=int, help='Paràmetre "a"')
    # choices argument
#    parser.add_argument('-o', '--operacio', 
#            type=str, choices=['suma', 'resta', 'multiplicacio'],
#            default='suma', required=False,
#            help='')

    args = parser.parse_args()
    # call a value: args.operacio or args.filename.



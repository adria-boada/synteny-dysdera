#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# resquitllar.py
#
# 14 Nov 2023  <adria@molevol-OptiPlex-9020>

help_msg = """
Describe here the goal, input, output, etc. of the script
as a multi-line block of text.
"""

import sys, math
# Reading input files' size and creating folders of PNG plots
import os

# Handling input as DataFrames
import pandas as pd, numpy as np
# Plotting figures
import seaborn as sns, matplotlib.pyplot as plt

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


prepare_plotting_folder(file, folder=None):
    """
    Makes a new directory if the folder path does not already exist.
    Concatenates the filename (which should include the extension) with the
    folder path.

    Input
    -----

    + file: The filename where the newly created plot should be written to.

    + folder [optional]: The folder name where the newly created plot should be
      written to. If it is not specified the function will return the
      filename, only.

    Output
    ------

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

    def __init__(self, files_dict, seqsizes_dict):
        """
        The class `Repeats` reads RepeatMasker's output, standardises the repeat
        type classification and returns summary tables of repeat content.

        Input
        -----

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
        """
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


    """
    ### LLEGIR INDEX-FILE ###

    Reads an index file in TSV (Tab-Separated Values) format in which sequence
    IDs are paired with their lengths in base pairs. The header (first row)
    contains the species name. For instance:

    Species name
    sequence-1\t32
    sequence-3\t56
    etc.

    Input
    -----

    + filename: The path to the previously explained TSV file.

    Output
    ------

    Returns a dictionary pairing sequids with their sizes
    """
    # script.py --rmlist out1 out2 --idxlist idx1 idx2 --out???

    # second script that does img???



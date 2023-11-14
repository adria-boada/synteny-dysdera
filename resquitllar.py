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


class Printing...


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



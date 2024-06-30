#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# boxplots.py
#
# 25 de juny 2024  <adria@molevol-OptiPlex-9020>

help_msg = """
Describe here the goal, input, output, etc. of the script
as a multi-line block of text.
"""

import sys
import re
import pandas as pd
import matplotlib.pyplot as plt, seaborn as sns

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def line_to_header_items(line: str):
    header = line.strip("\n").split("\t")
    # All fields are tabulated, except for the last field of parameter
    # names. Match it against 'inside parentheses' and split it by ",".
    header = header[:-1] + \
        re.findall(r"\((.*)\)", header[-1])[0].split(",")
    # Strip blank spaces from parameter names.
    header = [hd.strip(" ") for hd in header]
    
    return header

def isinteger(i):
    try:
        int(i)
        return True
    except ValueError:
        return False

def isfloat(i):
    try:
        float(i)
        return True
    except ValueError:
        return False

def reading_prefix_model_optimized_txts(file: str):
    """
    file : str
    The tabulated output produced by 'moments_Run_Optimizations.py'
    """
    # First, skim through the file and find all possible headers.
    # (ragged TSV).
    with open(file, "r") as fl:
        possible_headers = set()
        for line in fl:
            # If 'line' contains "log-likelihood", it's a header.
            if "log-likelihood" in line:
                [possible_headers.add(i) for i in line_to_header_items(line)]

    # Initialise answer dataset.
    answer = dict()
    for hd in possible_headers:
        answer[hd] = list()

    with open(file, "r") as fl:
        for line in fl:
            # If 'line' contains "log-likelihood", it's a header.
            if "log-likelihood" in line:
                header = line_to_header_items(line)
                # Detect keys not included in this 'header' portion of the file.
                not_header = list(
                    filter( lambda x: x not in header, answer.keys()))
            # Else, 'line' is not a header.
            else:
                fields = line.strip("\n").split("\t")
                # The last field '[-1]' contains
                # a list of subfields separated by commas.
                fields = fields[:-1] + fields[-1].split(",")
                for hd, fd in zip(header, fields):
                    # Manage integed/float/str fields with try blocks (bad?).
                    if isinteger(fd):
                        answer[hd].append(int(fd))
                    if isfloat(fd):
                        answer[hd].append(float(fd))
                    else:
                        answer[hd].append(str(fd))

                for nhd in not_header:
                    answer[nhd].append(None)

    return answer

def split_aggregated_columns(prefix_model_optimized: dict):
    """
    prefix_model_optimized : dict
    The dictionary returned by the function
    'reading_prefix_model_optimized_txts' found within this script.
    """
    columna_rep = prefix_model_optimized["Replicate"]
    columna_mod = prefix_model_optimized["Model"]
    # Retalla la informació interessant de les columnes conjugades.
    rounds = [field.split("_")[1] for field in columna_rep]
    replicates = [field.split("_")[3] for field in columna_rep]
    model = ["_".join(field.split("_")[:-1]) for field in columna_mod]
    execs = [field.split("_")[-1] for field in columna_mod]
    # Store them.
    prefix_model_optimized["Replicates"] = replicates
    prefix_model_optimized["Rounds"] = rounds
    prefix_model_optimized["Models"] = model
    prefix_model_optimized["Exec"] = execs

    return prefix_model_optimized

def main_categorical_plots(file: str):
    """
    file : str
    The tabulated output produced by 'moments_Run_Optimizations.py'
    """
    prefix_model_optimized = reading_prefix_model_optimized_txts(file)
    prefix_model_optimized = split_aggregated_columns(prefix_model_optimized)

    # Create a list with column names.
    columns = list(prefix_model_optimized.keys())
    # Create a list with the column names of numerical columns only.
    numerical_cols = list()
    for col in columns:
        if col not in ["Model", "Models", "Rounds",
                       "Replicate", "Replicates", "Exec"]:
            numerical_cols.append(col)

    # Printing the dataset that has been read.
    print(pd.DataFrame(prefix_model_optimized))
    print(pd.DataFrame(prefix_model_optimized).info())
    # Plotting the prior dataset.
    for numerical_var in numerical_cols:
        sns.catplot(
            data=prefix_model_optimized,
            y=numerical_var, x="Rounds", hue="Exec",
            col="Models", # Columnar faceting (also, see 'row').
        )
        plt.show()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Si es crida com a script:
if __name__ == '__main__':
    file = str(sys.argv[1])
    main_categorical_plots(file)


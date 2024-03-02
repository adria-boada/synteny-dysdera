#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# resquitllar.py / repellion.py
#
# repeatmasker_parser_v8
#
# 28 de febr. 2024  <adria@molevol-OptiPlex-9020>

help_msg = """
RM.out files
------------

RepeatMasker output file. RepeatMasker files from multiple species can be feed
to this script/parser at once. Each RepeatMasker file must be associated with an
species label. The script also requires an index file which specifies
chromosome/scaffolds groups with the sum of their sequences basepairs.

Index files
-----------

TSV (Tab-Separated Values) files with two columns: (1) a regex referencing a
group of scaffold/chromosome IDs, and (2) the sum of basepairs for all the
scaffold/chromosomes in the aforementioned group. The regex must be compliant
with `pandas.DataFrame.str.contains()`, which uses `re` regexes.

The first row must contain the species label. The regexes will try to match
values in the `sequid` column case-insensitively. For example:

```idxfile
species
chr.*1\t32
chr2\t5600
Chromosome10\t600
scaff|HRSCAF\t$(sum_len_scaffolds)
```

Where '.*' means one or more of any characters, '|' is the OR operation, etc. Be
careful with regexes; 'chr1' will match the string "chr1" and "chr10".

Formatting the RM.out files
---------------------------

The output obtained from RepeatMasker must be modified before using it. Some
rows contain an asterisk (*) as a 16th field, while others are empty. This
inconsistency makes it impossible for the file to be read by
`pandas.read_table()`.

The function `main_preproc` should create a more usable intermediary file. It
will include the length of the annotated RE as `bp_naive`. Moreover, it will
take into account overlaps and compute `bp_algor`, which is the length of the
RE discounting overlaps. The highest-scoring RE will have their overlapping
basepairs included, while the lower-scoring REs will have their overlapping
basepairs excluded.

Make sure to correctly label minor scaffolds and chromosomes in column 5
(`sequid` column). You can use `awk` or `sed`:

```bash
awk '$5=="Scaffold_1" {$5="DtilchrX"}' RM.out
```

Lastly, it is possible to find missing fields in the column "ID" and asterisks
in the column "position (left)". These obstruct the loading of the file. One
solution would be to substitute them with zeroes. Both columns "ID" and
"position (left)" will later be dropped and not used.
"""

# Reading input files' size and creating folders with figures.
import os

# Handling input as DataFrames.
import pandas as pd
# Plotting figures.
import seaborn as sns, matplotlib.pyplot as plt
# Measure the time it takes to run this script.
from datetime import datetime
import locale
locale.setlocale(locale.LC_TIME, '') # sets locale to the one used by user

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class Printing(object):
    """
    Print a given `message` with a formatting of interest (eg. as a warning,
    status, error...). Includes date and time of printing.

    Examples of use
    ---------------

    Printing("I love fries").status()
    Printing("Don't eat too many fries").warning()
    """
    def __init__(self, message):
        """
        Add the date and time of day to all the messages.
        """
        time = datetime.now().strftime("%d %b, %H:%M")
        self.message = f"{time}: {message}"

    def status(self):
        """
        Print in green and with preceding "STATUS" string
        """
        message = f"STATUS:{self.message}"
        # Change the colour to green (through ANSI colours in the terminal)
        message = "\033[0;32m" + str(message) + "\033[0;0m"
        print(message)
        return message

    def warning(self):
        """
        Print in yellow and with preceding "WARNING" string
        """
        message = f"WARNING:{self.message}"
        # Change the colour to yellow (through ANSI colours in the terminal)
        message = "\033[0;33m" + str(message) + "\033[0;0m"
        print(message)
        return message

    def error(self):
        """
        Print in red and with preceding "ERROR" string
        """
        message = f"ERROR:{self.message}"
        # Change the colour to red (through ANSI colours in the terminal)
        message = "\033[0;31m" + str(message) + "\033[0;0m"
        print(message)
        return message

def preprocess_rmout_asterisks(
    input_rmout: str,
    output_modified_rmout: str, ):
    """
    Preprocesses RM.out files.

    RM.out files contain a 16th field which can be either empty or can contain
    an asterisk. An asterisk indicates, citing RepeatMasker docs, that there
    exists a higher-scoring match whose domain partly (<80%) includes the domain
    of this present match. However, the inconsistency of the 16th field
    (empty/asterisk) makes it difficult to read with pandas.

    A solution is to substitute "asterisks" by "True" and "empty" by "False",
    naming that column "overlapping".

    Input
    -----

    + input_rmout: String pointing to an "RM.out" file.

    + output_modified_rmout: String where a new file will be written to. If
      there exists a file with the same filename, it will be overwritten.
    """
    Printing("Preprocessing the RM output file '"+
             str(input_rmout)+"'").status()
    with open(input_rmout) as inp, open(output_modified_rmout, "w") as out:
        # Skip first three lines (strangely formatted header).
        for i in range(3):
            out.write(inp.readline())
        # Process the rest of the file.
        for num, line in enumerate(inp):
            # Find whether the 16th field is empty or contains asterisk.
            overlapping_bool = (line.split()[-1] == "*")
            # Assemble the new output line.
            out_line = str(line).strip("\n*") + "\t" +\
                str(overlapping_bool) + "\n"
            # Some rows might not have all 16 fields! (eg. missing ID column)
            if len(out_line.split()) != 16:
                Printing("The amount of fields/columns in the line number "+
                         str(int(num) + 4)+
                         " is not 16. Make sure that this line is "+
                         "not missing any fields.").error()
                Printing("The problematic line is..."+"\n"+str(line)).error()
                # Stop writing file if an error is encountered.
                return None
            out.write(out_line)

def read_preprocessed_rmout(
    preproc_rmout: str, ):
    """
    """
    Printing(f"Reading preprocessed RM file {preproc_rmout}.").status()
    # Print file size before reading.
    file_size = os.path.getsize(preproc_rmout)
    Printing(f"Size in bytes: {file_size}.").status()
    # Load the file as pd.DataFrame.
    df = pd.read_table(
        preproc_rmout,
        header=None, # RM.out does not contain header
        skiprows=3,  # Skip the first three rows of header.
        sep="\s+",   # Fields are separated by one or multiple 'spaces'.
        names=[
        # For an explanation of RepeatMasker's fields, search its manual.
        "score_SW", "perc_divg", "perc_del", "perc_ins", "sequid",
        "begin", "end", "left", "orient", "name",
        "default_repclass", "begin_match", "end_match",
        "left_match", "ID", "overlapping"],
    )
    # Decrease memory usage by using efficient 'dtypes'. See:
    # <https://pandas.pydata.org/docs/user_guide/scale.html>)
    # Dropping unnecessary columns also improves memory usage:
    df = df.drop(columns=["left", "begin_match", "ID",
                          "end_match", "left_match"])
    # Object dtypes are more expensive than category dtypes. Try to downcast all
    # eligible object columns to categorical columns.
    obj_cols_to_categ = ["sequid", "orient", "name", "default_repclass"]
    df[obj_cols_to_categ] = df[obj_cols_to_categ].astype("category")

    # Create a new column in which element length (in basepairs) will be
    # computed. RE length could be used as a proxy for divergence. However,
    # substitutions (ie column "perc_divg") is a more straightforward measure;
    # it should be more accurate and direct.
    df["bp_naive"] = abs(df["end"] - df["begin"]) +1

    return df

def algor_overlap_basepairs(
    df: pd.DataFrame(), ):
    """
    Compute the length of repeats in basepairs while accounting for REs
    overlapping. If two REs overlap, take "score_SW" into consideration,
    removing the lower-scoring RE.
    """
    # Initialise a list to temporally store the results of the algorithm.
    answer = list()
    # For each individual sequence (sequid) in the pd.DataFrame():
    for seq in df["sequid"].unique():
        # Filter the main df by each unique `seq`.
        mask_sequid = df["sequid"] == seq
        # Acquire the list of intervals (begin, end) for the current `sequid`,
        # which will be feed to the algorithm. Keep in mind that two pairs of
        # coordinates cannot overlap of they are in different sequences!
        intervals = list(zip(
            list(df.loc[mask_sequid, "begin"]),
            list(df.loc[mask_sequid, "end"]),
            list(df.loc[mask_sequid, "score_SW"]),
            list(df.loc[mask_sequid,].index) ))
        # Sort the interval list by the fourth field (index).
        intervals = sorted(intervals, key=lambda x: x[3])
        # Produce a list consisting in (index, basepairs) tuples.
        answer.extend(remove_overlap_intervals(intervals))

    # Lastly, create a pd.Series from the (index, basepairs) tuples.
    index = [x[0] for x in answer]
    values = [x[1] for x in answer]

    return pd.Series(values, index=index)

def remove_overlap_intervals(
        intervals: list, ):
    """
    Take a list of intervals and remove overlapping basepairs based on the
    Smith-Waterman (SW) score of annotated REs.

    Input
    -----

    + intervals: A list of intervals which, in turn, are sublists themselves.
      Intervals are composed of `begin` and `end` coordinates as the first and
      second items. The pair of coordinates must be integers. Smith-Waterman
      score and an index (df.index) must be provided as the third and fourth
      items in the sublist. For instance:

      [ [beg_1, end_1, score_1, idx_1], ...,
        [beg_N, end_N, score_N, idx_N] ]

    Output
    ------

    Returns a list of tuples. The first item in the tuple is an index, while the
    second item is an amount of non-overlapping basepairs.
    """
    # Init a dict to store and return the results.
    answer = dict()
    # Create a list of points from the list of intervals. For instance,
    # interval-list -> ( [coord_beg, coord_end], [coord_beg, coord_end] )
    # point_list    -> ( [coord, beg, 0], [coord, end, 0], [coord, beg, 1] )
    point_list = []
    for i in range(0, len(intervals)):
        # Left/starting points labeled `0`:
        point_list.append(list( [intervals[i][0], 0, i]))
        # Right/ending points labeled `1`:
        point_list.append(list( [intervals[i][1], 1, i]))
        # Create a dict entry to store the results.
        answer[intervals[i][3]] = 0
    # Sort the list of points by their position (first field of sublists).
    point_list = sorted(point_list, key=lambda x: x[0])
    # Initialise other algorithm variables...
    currentBest = -1
    open_intervals = list()

    # Start iterating across the list of points. Discount overlapping
    # basepairs of REs with the poorer Smith-Waterman (SW) score.
    for point in point_list:
        # If the loop landed on a left point (0), open an interval.
        if point[1] == 0:

            # If there is no other interval opened:
            if currentBest == -1:
                # Then, the current-best-open interval is this one.
                currentBest = int(point[2])
                currentBegin = int(point[0])
                currentScore = int(intervals[currentBest][2])

            # Else, there already was an open interval. Check whether it is
            # of lower or higher score than currentBest.
            elif currentScore > intervals[point[2]][2]:
                # Do not close currentBest, because the next interval is of
                # lower score; however, append to open intervals list.
                open_intervals.append(list([
                    int(intervals[point[2]][2]),  # Score
                    point[2]]))                   # ID
                # Sort open intervals by score (first field).
                open_intervals = sorted(open_intervals,
                                        key=lambda x: x[0],
                                        reverse=True)
            else:
                # The next interval is of higher score. Compute length up
                # until now.
                interval_length = point[0] -currentBegin
                answer[intervals[currentBest][3]] += interval_length
                # Append index of past Best to `open_intervals`.
                pastBest = currentBest
                open_intervals.append(list([
                    int(intervals[pastBest][2]),  # Score
                    pastBest]))                   # ID
                # Sort open intervals by score (first field).
                open_intervals = sorted(open_intervals,
                                        key=lambda x: x[0],
                                        reverse=True)
                # Lastly, update currentBest.
                currentBest = int(point[2])
                currentBegin = int(point[0])
                currentScore = int(intervals[currentBest][2])

        # If the loop landed on a right point (1), close a currently open
        # interval.
        elif point[2] == currentBest:
            #> DEBUG
            #> print('point[2]:', point[2], ' == currentBest')
            # Compute length up until now:
            interval_length = point[0] - currentBegin +1
            answer[intervals[point[2]][3]] += interval_length
            # Close ` currentBest` interval. If there is any `pastBest`
            # interval in `open_intervals`, open the first one of the list
            # (the next best one, because it was sorted by score).
            if len(open_intervals) == 0:
                currentBest = -1
            else:
                # Assign the values of the ` pastBest` interval.
                currentScore, currentBest = open_intervals.pop(0)
                currentBegin = point[0] +1

        else:
            # The point we want to close is not `currentBest`, but one
            # of the points in `open_intervals`.
            open_intervals.remove(list([
                int(intervals[point[2]][2]),  # Score
                point[2]]))                   # ID

    return list(zip(answer.keys(), answer.values()))

def main_preproc(
    input_rmout: str,
    output_modified_rmout: str, ):
    """
    """
    # Preprocess RM.out file.
    preprocess_rmout_asterisks(input_rmout, output_modified_rmout)
    # Read it as pd.DataFrame()
    df = read_preprocessed_rmout(output_modified_rmout)
    # Create a new column with non-overlapping RE basepairs.
    df["bp_algor"] = algor_overlap_basepairs(df)
    Printing("Exporting the preprocessed/modified RM.out file to '"+
             str(os.getcwd()) + "/" + str(output_modified_rmout) + "'")\
        .status()
    df.to_markdown(output_modified_rmout, tablefmt="plain", index=False)

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



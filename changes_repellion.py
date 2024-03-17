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
species label.

Index files
-----------

TSV (Tab-Separated Values) files with three columns: (1) name of genomic subset
(eg. autosomes, sexual_chr, minor_scaffolds), (2) regex pattern which matches
the sequids found in this subset in the `sequid` column of the RM.out file, and
(3) the sum of lengths of all scaffold/chromosomes in this group, in basepairs.

The regex must be compliant with `pandas.DataFrame.str.contains()` function,
which uses regexes from the module `re`. The regexes will try to match values in
the sequid column **case-insensitively**. For example:

```idxfile
genome     .*             50900
sex_chr    chrx           10000
autosomes  chr\d          40000
min_scaff  scaff|ctg        900
```

The regex `.*` of the "genome" subset will match the whole genome - all
scaffolds and chromosomes will be included in this subset. The regex `chr\d`
will match `chr` followed by any possible digit (i.e. numbered autosomes).
The pipe character `|` is the OR operator. Consequently, the regex
`scaff|ctg` will match either `scaff` or `ctg`.

Formatting the RM.out files
---------------------------

The output obtained from RepeatMasker must be modified before using it. Some
rows contain an asterisk (*) as a 16th field, while others are empty. This
inconsistency makes it impossible for the file to be read by
`pandas.read_table()`.

The function `main_preproc` (setting `--preproc`) should create a more usable
intermediary file. It will include the length of the annotated RE as `bp_naive`,
which is computed as the absolute difference between "end" and "begin", plus one
(as RepeatMasker is zero indexed). Moreover, it will take into account overlaps
and compute `bp_algor`, which is the length of the RE discounting overlaps. The
highest-scoring RE will have their overlapping basepairs included, while the
lower-scoring REs will have their overlapping basepairs excluded.

If an index file will be provided, make sure to correctly label minor scaffolds
and chromosomes in the column 5 (`sequid` column). The user can use `awk` or
`sed`:

```bash
awk '$5=="Scaffold_1" {$5="DtilchrX"}' RM.out

sed 's/Scaffold_1\s/DtilchrX/' RM.out
```

Lastly, it is possible to find missing fields in the column "ID" and asterisks
in the column "position (left)". These obstruct the loading of the file. One
solution would be to substitute them with zeroes or "NA". Both columns "ID" and
"position (left)" will later be dropped and not used, so it does not matter how
these missing values are filled.
"""

# Exit the script, stop running...
import sys
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

def files_exist(
    files: str|list, ):
    """
    """
    if type(files) == list:
        for fl in files:
            # Try to open each file.
            try:
                with open(fl) as f:
                    pass
                # If any file is not found, return False.
            except FileNotFoundError:
                Printing("The file "+str(fl)+
                         " was not found.").error()
                return False
        # Else, all files are found, return True.
        return True
    # Otherwise, files must be a single str.
    elif type(files) == str:
        try:
            with open(files) as f:
                pass
        except FileNotFoundError:
            Printing("The file "+str(fl)+
                     "was not found.").error()
            return False
        # Else, the file was found, return True.
        return True

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
        dtype={
            "score_SW": int, "perc_divg": float, "perc_del": float,
            "perc_ins": float, "sequid": "category", "begin": int,
            "end": int, "left": object, "orient": "category",
            "name": "category", "default_repclass": "category",
            "begin_match": object, "end_match": object,
            "left_match": object, "ID": object, "overlapping": bool}
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
    Printing("Preprocessing the RM output file '"+
             str(input_rmout)+"'").status()
    preprocess_rmout_asterisks(input_rmout, output_modified_rmout)
    Printing("Reading preprocessed RM file "+
             f"{output_modified_rmout}.").status()
    # Print file size before reading.
    file_size = os.path.getsize(output_modified_rmout)
    Printing(f"Size in bytes of {output_modified_rmout}: "+
             f"{file_size}.").status()
    df = read_preprocessed_rmout(output_modified_rmout)
    Printing("Computing overlapping basepairs with algorithm "+
             "(time expensive).").status()
    df["bp_algor"] = algor_overlap_basepairs(df)
    Printing("Exporting the preprocessed/modified RM.out file to '"+
             str(os.getcwd()) + "/" + str(output_modified_rmout) + "'")\
        .status()
    df.to_markdown(output_modified_rmout, tablefmt="plain", index=False)

def reading_cmdline_rmout_and_idxfile(
    parameter: list, ):
    """
    """
    # Init a dict to store {"Species": "file"}
    pairs_species_file = dict()
    for pair in parameter:
        try:
            species, file = pair.split("=")
        except ValueError:
            Printing(f"The parameter `{pair}` of --rmout could "+
                     "not be split into two with an equal "+
                     "character (=).").error()
            Printing("Make sure all parameters are formatted "+
                     "as 'SpeciesLabel=RM.out.file'").error()
            return None
        # Store the pair in the dictionary.
        pairs_species_file[species] = file

    return pairs_species_file

class Repeats(object):
    """
    """
    def __init__(
        self,
        pairs_species_rmout: dict,
        # The index file is not always supplied.
        pairs_species_idxfile: dict=dict(), ):
        """
        """
        # First, make sure all provided files exist:
        list_of_files =\
            list(pairs_species_rmout.values()) +\
            list(pairs_species_idxfile.values())
        for file in list_of_files:
            # Try to open all files.
            try:
                with open(file) as f:
                    pass
            # If any file is not found, stop init execution.
            except FileNotFoundError:
                Printing("The file "+str(file)+
                         " was not found.").error()
                return None
        # Check if there are MORE index files than rmout files:
        if len(pairs_species_idxfile.values()) >\
                len(pairs_species_rmout.values()):
            Printing("There are more index files than RM.out "+
                     "files. Make sure there is at least one "+
                     "RM.out file per index file.").error()
            return None
        # Check if there is one rmout file per index file.
        species_idx = pairs_species_idxfile.keys()
        species_rmout = pairs_species_rmout.keys()
        check = [sp in species_idx for sp in species_rmout]
        # If one rmout does not have idx file:
        if not all(check):
            Printing("An RM.out is missing an index file. "+
                     "Index files are required to compute "+
                     "the non-repetitive fraction in the "+
                     "summary tables.").warning()
            # Create a variable to unable the execution of the function which
            # generates summary tables.
            self.correct_idx_files = False
        else:
            self.correct_idx_files = True

        # READING INPUT FILES.
        if self.correct_idx_files:
            Printing("Reading the provided index files...").status()
            df_index = self.read_index_files(pairs_species_idxfile)
            print(df_index)
        Printing("Reading the provided RM.out files...").status()
        df_rmout = self.read_modified_rmout(pairs_species_rmout)

        # UNPACK, STANDARDISE AND RECLASSIFY DEFAULT RE TYPES.
        Printing("Reclassifying RE default types into "+
                 "standardised types...").status()
        df_rmout = self.reclassify_default_reptypes(df_rmout)

        # CHECKING INFO AND RECLASSIFICATION.
        Printing("Printing information/columns about the created "+
                 "dataframe:").status()
        df_rmout.info(memory_usage="deep")
        Printing("Checking the reclassification/standardisation of "+
                 "repeat types:").status()
        self.df_check_reclass = self.check_reclassification(df_rmout)
        print(self.df_check_reclass.to_markdown(tablefmt="plain"))

        # CHECKING OVERALL OVERLAPPING BASEPAIRS.
        Printing("Checking overlapping REs:").status()
        self.df_check_overlap = self.evaluate_whole_overlapping(df_rmout)
        print(self.df_check_overlap.to_markdown(tablefmt="plain"))

        # CHECKING THE MATCHES BETWEEN INDEX SUBSETS AND MAIN DF SEQUIDS.
        if self.correct_idx_files:
            Printing("Checking the provided subsets in the index file:").status()
            self.df_check_matching_subsets = self.check_genome_subsets(
                df_index=df_index, df_rmout=df_rmout)
            print(self.df_check_matching_subsets.to_markdown(tablefmt="plain"))
            self.df_index = df_index

        self.df_rmout = df_rmout

    def read_index_files(
        self,
        pairs_species_idxfile: dict, ):
        """
        Read all the index files and concatenate them in a single pd.DataFrame.
        """
        answer_df_list = list()
        for species, idxfile in pairs_species_idxfile.items():
            df = pd.read_table(
                idxfile,
                sep="\s+",
                header=None,
                names=["Subset", "Pattern", "Bpsum"])
            df["Species"] = species
            answer_df_list.append(df)

        return pd.concat(answer_df_list).reset_index(drop=True)

    def read_modified_rmout(
        self,
        pairs_species_rmout: dict, ):
        """
        Read all the RM.out files and concatenate them in a single pd.DataFrame.
        """
        answer_df_list = list()
        for species, rmout in pairs_species_rmout.items():
            file_size = os.path.getsize(rmout)
            print(f" * {rmout}: {file_size} bytes")
            df = pd.read_table(
                rmout,
                header=0,
                sep="\s+")
            df["Species"] = species
            # Decrease memory usage by using efficient 'dtypes'.
            obj_cols_to_categ = ["sequid", "orient", "name",
                                 "default_repclass", "Species"]
            # Change 'obj' columns to 'category'.
            df[obj_cols_to_categ] = df[obj_cols_to_categ].astype("category")
            answer_df_list.append(df)

        return pd.concat(answer_df_list).reset_index(drop=True)

    def reclassify_default_reptypes(
        self,
        df: pd.DataFrame(), ):
        """
        Very long function with many paragraphs. Each paragraph reclassifies a
        default repeat type into four standardised clades. It helps in
        standardising the analyses.

        Some repeat types are reclassified based on an exact match in the field
        "default_repclass", while others are reclassified based on containing a
        substring/pattern.
        """
        # Unclassified
        df.loc[
            df["default_repclass"]=="__unknown",
            ["class", "subclass", "order", "superfam"]
        ] = ("Unclassified", "NA", "NA", "NA")

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
            (df["default_repclass"].str.contains("^LTR/")) |
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
        ] = ("Retrotransposon", "NA", "YR", "DIRS")
        df.loc[
            df["default_repclass"].str.contains("LTR") &
            df["default_repclass"].str.contains("Ngaro"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "YR", "Ngaro")
        df.loc[
            df["default_repclass"].str.contains("LTR/Viper"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "YR", "Viper-group")
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
        ] = ("DNA", "1", "DDE", "Zator")
        # Both MITE and nMITE are unassignable to "taxons"; these attributes are
        # not synapomorphic (unique to any one clade of repeats). I'd rather
        # classify them as no further than "DNA" (all else unknown).
        df.loc[
            df["default_repclass"].str.contains("ClassII_MITE"),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "NA", "NA")
        df.loc[
            df["default_repclass"].str.contains("ClassII_nMITE"),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "NA", "NA")
        df.loc[
            df["default_repclass"].str.contains("DNA") &
            df["default_repclass"].str.contains("Academ"),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "DDE", "Academ")
        df.loc[
            (df["default_repclass"].str.contains("DNA") &
            df["default_repclass"].str.contains("CMC")) |
            df["default_repclass"].str.contains("ClassII_DNA_CACTA"),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "DDE", "CACTA")
        df.loc[
            (df["default_repclass"].str.contains("DNA") &
            df["default_repclass"].str.contains("Dada")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "DDE", "Dada")
        df.loc[
            (df["default_repclass"].str.contains("DNA") &
            df["default_repclass"].str.contains("Kolobok")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "DDE", "Kolobok")
        df.loc[
            (df["default_repclass"].str.contains("DNA") &
            df["default_repclass"].str.contains("Ginger")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "DDE", "Ginger")
        df.loc[
            (df["default_repclass"].str.contains("DNA") &
            df["default_repclass"].str.contains("MULE")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "DDE", "MULE")
        df.loc[
            (df["default_repclass"].str.contains("DNA") &
            df["default_repclass"].str.contains("Merlin")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "DDE", "Merlin")
        df.loc[
            (df["default_repclass"].str.contains("DNA/P")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "DDE", "P")
        df.loc[
            (df["default_repclass"].str.contains("DNA/PIF")) |
            (df["default_repclass"].str.contains("ClassII_DNA_Harbinger")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "DDE", "PIF-Harbinger")
        df.loc[
            (df["default_repclass"].str.contains("PiggyBac")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "DDE", "PiggyBac")
        df.loc[
            (df["default_repclass"].str.contains("DNA/Sola")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "DDE", "Sola")
        df.loc[
            (df["default_repclass"].str.contains("DNA/Sola")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "DDE", "Sola")
        df.loc[
            (df["default_repclass"].str.contains("DNA/TcMar") |
            df["default_repclass"].str.contains("DNA_TcMar")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "DDE", "TcMar")
        df.loc[
            (df["default_repclass"].str.contains("DNA") &
            df["default_repclass"].str.contains("hAT")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "DDE", "hAT")
        df.loc[
            (df["default_repclass"].str.contains("DNA_Mutator")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "DDE", "Mutator")
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

        # Object dtypes are more expensive than category dtypes. Try to
        # downcast from objects to categories, if possible.
        obj_cols_to_categ = ["class", "subclass", "order", "superfam", "mite"]
        df[obj_cols_to_categ] = df[obj_cols_to_categ].astype("category")

        return df

    def check_reclassification(
        self,
        df_rmout: pd.DataFrame(), ):
        """
        Checks how the default RE types/names have been reclassified.
        Returns a dataframe pairing every unique default type/name with their
        assigned class, subclass, order and superfamily.

        Input
        -----

        + df_rmout: the main dataframe created by the class `Repeats`.

        Output
        ------

        Returns a pandas dataframe in which default types are paired with their
        reassigned classification.
        """
        # Create a subset of the main dataframe by dropping duplicate classes,
        # subclasses, orders, and superfamilies.
        df_unique_classes = df_rmout.drop_duplicates(
            subset=["default_repclass", "class",
                    "subclass", "order", "superfam", "mite"])
        # Sort the dataframe by the newly assigned types. Reset its index.
        df_unique_classes = df_unique_classes.sort_values(
            by=["class", "subclass", "order", "superfam", "mite"],
            ignore_index=True)
        # We are only interested in columns with the repeat types;
        # drop the rest of columns.
        df_unique_classes = df_unique_classes[["class", "subclass", "order",
                                               "superfam", "mite",
                                               "default_repclass"]]

        return df_unique_classes

    def evaluate_whole_overlapping(
        self,
        df_rmout: pd.DataFrame(), ):
        """
        Compute the sum of base pairs obtained by the three different methods:

        1) Sum of all base pairs (column "bp_naive").
        2) Sum of the higher SW-score REs, removing the lower SW-score
           REs (column "bp_naive" filtered by column "overlapping").
        3) Sum of the higher SW-score REs, including the non-overlapping
           fragments of lower SW-score REs (column "bp_algor").

        Input
        -----

        + df_rmout: the main dataframe created by the class `Repeats`.

        Output
        ------

        Prints and returns a pandas dataframe with the sum of basepairs obtained
        from each method.
        """
        answer_df_list = list()
        for species in df_rmout["Species"].unique():
            mask_species = df_rmout["Species"] == species
            method_1_bpsum = df_rmout.loc[
                mask_species, "bp_naive"].sum()
            method_2_bpsum = df_rmout.loc[
                (mask_species) & (~ df_rmout["overlapping"]),
                              "bp_naive"].sum()
            method_3_bpsum = df_rmout.loc[
                mask_species, "bp_algor"].sum()
            # Create dataframe displaying the sum of basepairs per method.
            df = pd.DataFrame({
                "Species": [species]*3,
                "Basepairs": [
                    method_1_bpsum,
                    method_2_bpsum,
                    method_3_bpsum]},
                index=["naive", "best", "algor"])
            # Compute percentages relative to `method_1_bpsum`.
            df["%Comparison"] =\
                (df["Basepairs"] / method_1_bpsum) *100
            # Store each species' overlapping stats.
            answer_df_list.append(df)

        # Concatenate the dataframes of all "Species" into a single one.
        return pd.concat(answer_df_list)

    def check_genome_subsets(
        self,
        df_rmout: pd.DataFrame,
        df_index: pd.DataFrame, ):
        """
        Count how many rows match with the pattern referring to each provided
        genome subset. At the moment it is a poorly optimised function.
        """
        # Initialise return dataframe.
        answer = df_index[["Species", "Subset", "Pattern"]].copy()
        # The following lambda function computes number of rows matching a
        # pattern in the defined species.
        compute_nrows = lambda spec, pat:\
            df_rmout.loc[(df_rmout["sequid"].str.contains(pat, case=False)) &
                         (df_rmout["Species"]==spec)].shape[0]
        # Apply the function to each row (axis=1)
        answer["#Rows"] = answer.apply(lambda x:\
            compute_nrows(x["Species"], x["Pattern"]), axis=1)
        # The following lambda function computes percentage of matching rows.
        compute_perc_rows = lambda spec, pat:\
            (df_rmout.loc[(df_rmout["sequid"].str.contains(pat, case=False)) &
                          (df_rmout["Species"]==spec)].shape[0] /
             df_rmout.loc[df_rmout["Species"]==spec].shape[0]) * 100
        # Apply the function to each row (axis=1)
        answer["%Rows"] = answer.apply(lambda x:\
            compute_perc_rows(x["Species"], x["Pattern"]), axis=1).round(1)
        # The following lambda function computes the number of unique sequids
        # matching a pattern in the defined species.
        compute_nuniq = lambda spec, pat: len(list(
            df_rmout.loc[(df_rmout["sequid"].str.contains(pat, case=False)) &
                         (df_rmout["Species"]==spec), "sequid"].unique()))
        # Apply the function to each row (axis=1)
        answer["#Useq"] = answer.apply(lambda x:\
            compute_nuniq(x["Species"], x["Pattern"]), axis=1)
        # Once again, percentages of unique sequids...
        compute_perc_nuniq = lambda spec, pat: (len(
            df_rmout.loc[(df_rmout["sequid"].str.contains(pat, case=False)) &
                         (df_rmout["Species"]==spec), "sequid"].unique()) /
            len(df_rmout.loc[df_rmout["Species"]==spec,
                             "sequid"].unique())) *100
        answer["%Useq"] = answer.apply(lambda x:\
            compute_perc_nuniq(x["Species"], x["Pattern"]), axis=1)

        answer = answer[["Species", "Subset", "#Rows",
                         "%Rows", "#Useq", "%Useq"]]
        answer["%Rows"] = answer["%Rows"].round(2)
        answer["%Useq"] = answer["%Useq"].round(2)

        return answer

def repeats_content_histogram(
    dict_df_grouped: dict,
    hue_column: str, val_y_column: str, cat_x_column: str,
    yaxis_label: str, title: str,
    savefig: str=None, ):
    """
    + dict_df_grouped: a dict. where keys are titles. values are pd.Dataframes
      for plotting.

    + cat_x_column: column with species.

    + val_y_column: column with count of REs or basepairs.

    + hue_column: column which contains type of repeats. Try to make these
      unique (careful with DNA-NA and RT-NA, if grouped by subclass NA will
      group together).
    """
    # Crea tants `subplots` com items a la llista `dict_df_grouped`.
    fig, axes = plt.subplots(ncols=len(dict_df_grouped), sharey=True)
    if len(dict_df_grouped) == 1:
        axes = [axes]
    # Itera a través de parelles d'axes i dataframes:
    for ax, (df_title, df_loop) in zip(axes, dict_df_grouped.items()):
        g = sns.histplot(data=df_loop, ax=ax, # Utilitza `ax` del "loop".
                     x=cat_x_column, weights=val_y_column,
                     hue=hue_column, multiple="stack",
                     shrink=0.9, legend=True,
                     ).set(title=df_title, ylabel=yaxis_label)
    # Get legend in the last subplot.
    # print(g) # debug
    #handles, labels = ax.get_legend_handles_labels()
    #fig.legend(handles, labels, loc="upper left")
    fig.suptitle(title)

    if savefig:
        Printing(f"Creating a figure at {savefig}").status()
        plt.savefig(savefig, format="svg")
    else:
        plt.show()

    return None

def repeats_content_histo_wrapper(
    df: pd.DataFrame,
    dict_div_masks: dict,
    groupby_colnames: list,
    val_y_column: str,
    yaxis_relative: bool=True,
    title: str="",
    savefig: str=None):
    """
    dict_div_masks = {
        "All": df["perc_divg"] >= 0,
        "div < 1%": df["perc_divg"] < 1,
        "div < 5%": df["perc_divg"] < 5,
        "div >= 5%": df["perc_divg"] >= 5,
    }
    groupby_colnames = ["class"] | ["class", "subclass"]
    """
    # Initialise a dictionary to store pd.DataFrame.groupby().agg() objects.
    dict_df_grouped = {}
    # Iterate across the given df masks found in the `dict_div_masks` parameter:
    for key, mask_divg in dict_div_masks.items():
        # Locate the df rows matching `mask_divg`.
        dict_df_grouped[key] = df.loc[mask_divg] \
            .groupby(groupby_colnames, observed=True, as_index=False) \
            .agg(count=pd.NamedAgg(column=val_y_column,
                                          aggfunc="sum"))
        dict_df_grouped[key] = dict_df_grouped[key].rename(
            columns={"count": val_y_column})

        # Sort the df by `hue_column` (RE types). This will make it easy for
        # RE types to retain the same colour across subplots (otherwise might
        # change colours if the data/order between plots changes).
        dict_df_grouped[key] = dict_df_grouped[key].sort_values(
            by=groupby_colnames)

        # If the parameter `yaxis_relative` is equal to True, divide the column
        # with values by each Species' sum of the values...
        if yaxis_relative:
            d = dict_df_grouped[key]
            # Silence a DeprecationWarning that will happen during the next
            # paragraph...
            d[val_y_column] = d[val_y_column].astype("float")
            for species in d["Species"].unique():
                mask_species = d["Species"]==species
                d.loc[mask_species, val_y_column] = (
                    d.loc[mask_species, val_y_column] /
                    d.loc[mask_species, val_y_column].sum() )
            yaxis_label="bp density"
        else:
            yaxis_label="bp count"
    # Call the plotting function.
    repeats_content_histogram(
        dict_df_grouped=dict_df_grouped,
        hue_column=groupby_colnames[-1], val_y_column=val_y_column,
        cat_x_column="Species",
        yaxis_label=yaxis_label,
        title=title,
        savefig=savefig)

    return None

def repeats_content_kiosk(
    repeats_instance: object,
    folder: str, ):
    """
    """
    # Shortcut to the main `repeats_instance` dataframe.
    df = repeats_instance.df_rmout

    # Contingut de les classes relatiu.
    repeats_content_histo_wrapper(
        df,
        dict_div_masks={
            "All": df["perc_divg"] >=0, },
        val_y_column="bp_algor",
        groupby_colnames=["Species", "class",],
        title="Class relative bp",
        savefig=prepare_plotting_folder(
            "1_class_relbp.svg", folder))
    # Contingut de les classes absolut.
    repeats_content_histo_wrapper(
        df,
        dict_div_masks={
            "All": df["perc_divg"] >=0, },
        val_y_column="bp_algor", yaxis_relative=False,
        groupby_colnames=["Species", "class",],
        title="Class absolute bp",
        savefig=prepare_plotting_folder(
            "1_class_absbp.svg", folder))

    # Contingut de les subclasses relatiu.
    repeats_content_histo_wrapper(
        df,
        dict_div_masks={
            "All": df["perc_divg"] >=0, },
        val_y_column="bp_algor",
        groupby_colnames=["Species", "class", "subclass"],
        title="Subclass relative bp",
        savefig=prepare_plotting_folder(
            "2_subclass_relbp.svg", folder))
    # Contingut de les subclasses absolut.
    repeats_content_histo_wrapper(
        df,
        dict_div_masks={
            "All": df["perc_divg"] >=0, },
        val_y_column="bp_algor", yaxis_relative=False,
        groupby_colnames=["Species", "class", "subclass"],
        title="Subclass absolute bp",
        savefig=prepare_plotting_folder(
            "2_subclass_absbp.svg", folder))

    # Contingut dels ordres relatiu.
    df_orders = df.copy()
    mask_dna = df_orders["class"]=="DNA"
    mask_retro = df_orders["class"]=="Retrotransposon"
    df_orders[["order", "superfam"]] =\
        df_orders[["order", "superfam"]].astype(str)
    df_orders.loc[mask_dna, "order"] = \
        "DNA " + df_orders.loc[mask_dna, "order"].astype(str)
    df_orders.loc[mask_retro, "order"] = \
        "RT " + df_orders.loc[mask_retro, "order"].astype(str)

    repeats_content_histo_wrapper(
        df_orders,
        dict_div_masks={
            "All": df_orders["perc_divg"] >=0, },
        val_y_column="bp_algor",
        groupby_colnames=["Species", "class", "subclass", "order"],
        title="Orders relative bp",
        savefig=prepare_plotting_folder(
            "3_orders_relbp.svg", folder))
    # Contingut dels ordres absolut.
    repeats_content_histo_wrapper(
        df_orders,
        dict_div_masks={
            "All": df_orders["perc_divg"] >=0, },
        val_y_column="bp_algor", yaxis_relative=False,
        groupby_colnames=["Species", "class", "subclass", "order"],
        title="Orders absolute bp",
        savefig=prepare_plotting_folder(
            "3_orders_absbp.svg", folder))

    # Quatre particions segons "edat" de les REs (dependent del valor de
    # divergència respecte RE consensus).
    repeats_content_histo_wrapper(
        df,
        dict_div_masks={
            "All": df["perc_divg"] >=0,
            "div < 1%": df["perc_divg"] < 1,
            "div < 5%": df["perc_divg"] < 5,
            "div >= 5%": df["perc_divg"] >= 5, },
        val_y_column="bp_algor",
        groupby_colnames=["Species", "class",],
        title="Class by divergence subsets",
        savefig=prepare_plotting_folder(
            "1_class_bydiv_relbp.svg", folder))
    repeats_content_histo_wrapper(
        df_orders,
        dict_div_masks={
            "All": df_orders["perc_divg"] >=0,
            "div < 1%": df_orders["perc_divg"] < 1,
            "div < 5%": df_orders["perc_divg"] < 5,
            "div >= 5%": df_orders["perc_divg"] >= 5, },
        val_y_column="bp_algor",
        groupby_colnames=["Species", "class", "subclass", "order"],
        title="Orders by divergence subsets",
        savefig=prepare_plotting_folder(
            "3_orders_bydiv_relbp.svg", folder))

    # Contingut de LTR, LINE, SINE, DDE a escala de superfamilia.
    residual_superfamilies = ["Ngaro"]
    important_orders = ["LTR", "LINE", "DDE", "SINE"]

    # Primer, elimina famílies residuals, com podrien ser "Ngaro", "DADA", etc.
    mask_not_residual = (~ df["superfam"].isin(residual_superfamilies))
    # Segon, localitza els ordres importants (LTR, LINE, SINE, etc.)
    for o in important_orders:
        df_filtered = df.loc[
            (df["order"]==o) &
            (mask_not_residual)].copy()
        df_filtered[["order", "superfam"]] =\
            df_orders[["order", "superfam"]].astype(str)
        repeats_content_histo_wrapper(
            df_filtered,
            dict_div_masks={
                "All": df_filtered["perc_divg"] >=0, },
            val_y_column="bp_algor",
            groupby_colnames=["Species", "class", "subclass",
                              "order", "superfam"],
            title=f"{o} relative bp",
            savefig=prepare_plotting_folder(
                f"4_{o}_superfamilies_relbp.svg", folder))
        repeats_content_histo_wrapper(
            df_filtered,
            dict_div_masks={
                "All": df_filtered["perc_divg"] >=0, },
            val_y_column="bp_algor", yaxis_relative=False,
            groupby_colnames=["Species", "class", "subclass",
                              "order", "superfam"],
            title=f"{o} absolute bp",
            savefig=prepare_plotting_folder(
                f"4_{o}_superfamilies_absbp.svg", folder))
        repeats_content_histo_wrapper(
            df_filtered,
            dict_div_masks={
                "All": df_filtered["perc_divg"] >=0,
                "div < 1%": df_filtered["perc_divg"] < 1,
                "div < 5%": df_filtered["perc_divg"] < 5,
                "div >= 5%": df_filtered["perc_divg"] >= 5, },
            val_y_column="bp_algor",
            groupby_colnames=["Species", "class", "subclass",
                              "order", "superfam"],
            title=f"{o} by divergence subsets",
            savefig=prepare_plotting_folder(
                f"4_{o}_superfamilies_bydiv_relbp.svg", folder))

def div_distrib_recursive_hist(
    df, value_column: str,
    categorical_columns: list,
    hueby_column: str=None,
    yaxis_relative: bool=True, multiple: str="layer",
    plot_style: str="hist",
    folder: str=None, prefix_title: str="", ):
    """
    Plot a column `value_column` of a pandas DataFrame filtering by many
    available categories.
    """
    # For each `category` in the zeroeth column of the list
    # `categorical_columns`:
    for category in df[categorical_columns[0]].unique():
        # Filter `df` by `category` within each for-loop.
        bool_filter = (df[categorical_columns[0]] == category)
        df_filtered = df.loc[bool_filter]
        # Obtain the title that will be used in the current plot.
        current_title = prefix_title + "_" + category
        # Try to find further categorical subdivisions in the list of
        # categorical columns. Count all combinations of unique categories.
        checkval = df_filtered.drop_duplicates(
            subset=categorical_columns)[categorical_columns].shape[0]
        # If `checkval` is above 1 there are further subdivisions of the `df`
        # (more than one combination). Run this function recursively at a lower
        # level (one category down the list) with the filtered `df`.
        if checkval > 1:
            div_distrib_recursive_hist(
                # Run with filtered df.
                df=df_filtered,
                value_column=value_column,
                # Run one category down the list (recursive aspect)
                categorical_columns=categorical_columns[1:],
                hueby_column=hueby_column,
                prefix_title=current_title, folder=folder,
                yaxis_relative=yaxis_relative, multiple=multiple)

        # Now, the plotting section.
        # If a hueby_column was specified, start by iterating across it.
        if hueby_column:
            # Init a DataFrame where value counts will be concatenated.
            df_plot = pd.DataFrame()
            for hue in df_filtered[hueby_column].unique():
                # Filter by `hue` categories and create a Series of
                # value_counts.
                valcounts_series = (
                    df_filtered.loc[df_filtered[hueby_column] == hue,
                                    value_column]
                    ).value_counts(normalize=yaxis_relative).sort_index(
                    ).to_frame()
                # Label these values with their hue and concatenate to main df.
                valcounts_series[hueby_column] = hue
                df_plot = pd.concat([df_plot, valcounts_series])

        # If a hueby_column was not specified...
        else:
            # Create a Series composed of value counts.
            df_plot = df_filtered[value_column].value_counts(
                normalize=yaxis_relative).sort_index().to_frame()

        # Lastly, with the value counts, plot them as an histogram.
        if plot_style == "hist":
            ax = sns.histplot(
                x=df_plot.index, weights=df_plot[value_column],
                hue=df_plot[hueby_column], binwidth=1,
                kde=True, kde_kws={"bw_adjust": .25},
                multiple=multiple, legend=True, # element="step",
                # Colors for each species (it should not be a static
                # variable, but parameterised by the function).
                palette={"Dcat": "#eb8f46", "Dtil": "#30acac",
                         "Dsil": "green", "Dcatv33": "#eb8f46",
                         "Dcatv35": "#eb8f46"},)
        elif plot_style == "kde":
            # Plots only the KDE plot without the underlying histogram.
            ax = sns.kdeplot(x=df_plot.index, weights=df_plot[value_column],
                        hue=df_plot[hueby_column], fill=True, bw_adjust=.25,
                        multiple=multiple, legend=True,
                         # Colors for each species (it should not be a static
                         # variable, but parameterised by the function).
                         palette={"Dcat": "#eb8f46", "Dtil": "#30acac",
                                  "Dsil": "green", "Dcatv33": "#eb8f46",
                                  "Dcatv35": "#eb8f46"},)
            # Avoid plotting values below 0% divergence. Useful for KDE plot
            # `sns.kdeplot()`. Nonetheless, `sns.histplot()` cuts axis
            # automatically. Sets only minimum limit, with `None` maximum limit.
            ax.set_xlim(0, None)

        plt.title(current_title)
        #ax.legend_.set_title("debug")#DEBUG
        # Reduce the linewidth of the patches to be minimalistic.
        for patch in ax.legend_.get_patches():
            patch.set_linewidth(0.5)

        # SVG format does not read `dpi`, only useful for PNG files.
        f = "svg"
        save_fig_as = prepare_plotting_folder(
            file=str(current_title) + "." + f,
            folder=folder)
        plt.savefig(save_fig_as, dpi=300, format=f)
        # Close the plt.axes or they will leak to the next figures.
        plt.close()

    return None

def div_distrib_kiosk(
    repeats_instance: object,
    folder: str, ):
    """
    """
    # Shortcut to the main `repeats_instance` dataframe.
    df = repeats_instance.df_rmout
    div_distrib_recursive_hist(
        df=df, value_column="perc_divg",
        categorical_columns=["class", "order", "superfam"],
        hueby_column="Species", yaxis_relative=True,
        multiple="dodge", plot_style="hist",
        folder=folder, prefix_title="Relative")
    div_distrib_recursive_hist(
        df=df, value_column="perc_divg",
        categorical_columns=["class", "order", "superfam"],
        hueby_column="Species", yaxis_relative=False,
        multiple="dodge", plot_style="hist",
        folder=folder, prefix_title="Absolute")

def div_boxplots_kiosk(
    repeats_instance: object,
    filename: str, ):
    """
    """
    # Shortcut to the main `repeats_instance` dataframe.
    df = repeats_instance.df_rmout.copy()
    # Filter residual TE superfamilies out.
    residual_superfamilies = [
        "Ngaro", "Dada", "Crypton", "5S", "7SL", "R2", "Merlin", "DIRS", "MULE",
        "Ginger", "L1"]
    df = df.loc[ ~ df["superfam"].isin(residual_superfamilies)]
    # Prepare a categorical column. Avoid pooling "NA" superfamilies from
    # different classes (DNA-NA, RT-NA, etc.).
    df["category"] = df["superfam"].astype("object")
    mask = (df["superfam"]=="NA") & (df["class"]=="Tandem_repeat")
    df.loc[mask, "category"] = "Tandem_repeat"
    mask = (df["superfam"]=="NA") & (df["class"]=="DNA")
    df.loc[mask, "category"] = "DNA Unclass."
    mask = (df["superfam"]=="NA") & (df["class"]=="Other")
    df.loc[mask, "category"] = df.loc[mask, "subclass"]
    mask = ((df["superfam"]=="NA") &
            (df["class"]=="Retrotransposon") &
            (df["order"]=="NA"))
    df.loc[mask, "category"] = "RT Unclass."
    mask = ((df["superfam"]=="NA") &
            (df["class"]=="Retrotransposon") &
            (df["order"]!="NA"))
    df.loc[mask, "category"] = "RT " + df.loc[mask, "order"].astype(str)
    mask = (df["superfam"]=="tRNA")
    df.loc[mask, "category"] = "SINE tRNA"
    df.loc[df["class"]=="Unclassified", "category"] = "Unclass."

    # Groupby "category" and "aggregate" by summing the column "bp_algor".
    order_boxes = list(df.groupby("category", observed=True, as_index=False) \
        .agg(bp_sum=pd.NamedAgg(column="bp_algor", aggfunc="sum")) \
        .sort_values("bp_sum", ignore_index=True, ascending=False) \
            ["category"])

    # Plots.
    ax = sns.boxplot(data=df, x="category", y="perc_divg",
                     order=order_boxes, hue="Species",
                     palette={"Dcat": "#eb8f46", "Dtil": "#30acac",
                              "Dsil": "green", "Dcatv33": "#eb8f46",
                              "Dcatv35": "#eb8f46",
                              "Dcat_MegaLTR": "yellow"},)

    # Put more space between categorical X ticks.
    tick_labels = ax.get_xticklabels()
    maxsize = max([t.get_window_extent().width for t in tick_labels])
    m = 0.2
    s = maxsize/plt.gcf().dpi*len(order_boxes)+2*m
    margin = m/plt.gcf().get_size_inches()[0]
    plt.gcf().subplots_adjust(left=margin, right=1.-margin)
    plt.gcf().set_size_inches(s, plt.gcf().get_size_inches()[1])

    # Tight layout, rotate ticks, manually modify figure size. Not necessary
    # with the better previous paragraph, which automatically changes figsize!
    plt.tight_layout()
    #>ax.set_xticklabels(ax.get_xticklabels(), rotation=60)
    #>plt.figure(figsize=(20, 10))

    if filename:
        Printing(f"Creating a figure at {filename}").status()
        plt.savefig(filename, dpi=300)
    else:
        plt.show()
    plt.close()

    return None

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Si es crida com a script:
if __name__ == '__main__':

    ### PARSE ARGUMENTS ###
    import argparse
    parser = argparse.ArgumentParser(
        description=help_msg,
        # make sure the 'help_msg' is not automatically
        # wrapped at 80 characters (manually assign newlines).
        formatter_class=argparse.RawTextHelpFormatter)

    # Opció per preprocessar el fitxer; incompatible amb la resta d'opcions.
    parser.add_argument("--preproc", nargs=2, dest="preproc_files",
                 metavar=("READ-RM.OUT",
                          "WRITE-RM.OUT"),
                 help="A pair of strings. The first string is the "+
                 "input filename of RepeatMasker to be preprocessed."+
                 "The second string is where the output file "+
                 "will be written to.")

    # Opció per analitzar fitxers preprocessats durant el pas anterior.
    parser.add_argument("--rmout", nargs="+",
                 metavar="Species=RM.out",
                 help="A list of preprocessed files containing "+
                 "the annotation of repeats obtained from "+
                 "RepeatMasker.")
    parser.add_argument("--idx", nargs="+",
                 metavar="Species=index.idx",
                 help="A list of index files, as specified in --help. "+
                 "The species labels must be congruent with the "+
                 "ones specified in the --rmout argument. It is not "+
                 "required if the user is not interested in "+
                 "comparing the repetitive fraction with the "+
                 "non-repetitive fraction, or the user is not "+
                 "interested in splitting the genome in data "+
                 "subsets by the 'sequid' column (i.e. "+
                 "chrX, autosomes, etc.)")

    # Anàlisis possibles amb els fitxers preprocessats.
    parser.add_argument("--content_plots", metavar="FIGURE-FOLDER",
                 help="Flag which will call the creation of "+
                 "content plots. Set a folder name where the "+
                 "newly created figure files will be stored")
    parser.add_argument("--distrib_div_hist",
                 metavar="FIGURE-FOLDER",
                 help="Flag which will call the creation of "+
                 "divergence distributions. Set a folder name "+
                 "where the newly created figure files will be "+
                 "stored")
    parser.add_argument("--div_boxplots",
                 metavar="FILENAME",
                 help="Flag which will call the creation of "+
                 "boxplots with the Y axis being 'divergence'. "+
                 "Set a filename where the newly created "+
                 "figure will be stored.")

    # Fitxer d'estil de matplotlib.
    parser.add_argument("--matplotlib_style",
                        help="A file with matplotlib style settings. "+
                        "They are normally stored in the folder "+
                        "$HOME/.config/matplotlib/stylelib/.")

    args = parser.parse_args()

    if args.preproc_files:
        read_rmout, write_rmout = args.preproc_files
        main_preproc(read_rmout, write_rmout)

    elif args.rmout:
        pairs_species_rmout =\
            reading_cmdline_rmout_and_idxfile(args.rmout)
        # Accepta IDX.
        if args.idx:
            pairs_species_idxfile =\
                reading_cmdline_rmout_and_idxfile(args.idx)
            # Before creating a `Repeats` instance, make sure all files exist.
            list_of_files =\
                list(pairs_species_rmout.values()) +\
                list(pairs_species_idxfile.values())
            if not files_exist(list_of_files):
                sys.exit()
            # Create the instance of `Repeats`.
            repeats = Repeats(
                pairs_species_idxfile=pairs_species_idxfile,
                pairs_species_rmout=pairs_species_rmout)
        # Tenir IDX no és necessari per a totes les situacions.
        else:
            # Before creating a `Repeats` instance, make sure all files exist.
            list_of_files =list(pairs_species_rmout.values())
            if not files_exist(list_of_files):
                sys.exit()
            repeats = Repeats(
                pairs_species_rmout=pairs_species_rmout)

        if args.matplotlib_style:
            plt.style.use(args.matplotlib_style)

        if args.content_plots:
            repeats_content_kiosk(repeats, folder=args.content_plots)

        if args.distrib_div_hist:
            div_distrib_kiosk(repeats, folder=args.distrib_div_hist)

        if args.div_boxplots:
            div_boxplots_kiosk(repeats, args.div_boxplots)

    else:
        Printing("Select either the option '--preproc' "+
                 "to preprocess the annotation files, or "+
                 "the option '--rmout' coupled with an "+
                 "available analysis setting, like "+
                 "'--content_plots'.").error()
        Printing("Calling the script with `-h` will show "+
                 "a help message with available "+
                 "settings.").error()


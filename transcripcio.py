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
from natsort import index_natsorted, humansorted
# Creating plots from `pandas` dataframes.
import matplotlib.pyplot as plt
import seaborn as sns
# Measure the time it takes to run this script
from datetime import datetime
import locale
locale.setlocale(locale.LC_TIME, '') # sets locale to the one used by user

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

def cigar_analysis(cigar_string):
    """
    Analyse a CIGAR string and compute its total amount of "M" (matches), "I"
    (insertions to query) and "D" (deletions to query).

    Input
    =====

    + cigar_string: A string of letters and numbers symbolizing the gap
      composition of an alignment/mapping. Accepts a string with "M", "I" and
      "D". Each letter must be preceded by a number (which is the length of the
      fragment). Take a look at the specifications of the "SAM" format for a
      more in-depth explanation:

      <https://samtools.github.io/hts-specs/SAMv1.pdf>

    Output
    ======

    Returns a dictionary with the sum of "M", "I" and "D". Moreover, it returns
    the amount of insertions and deletions instead of their lengths in
    basepairs. It is useful to compute "compressed gap identity/divergence".
    """
    # Split string by separating all letters with a space:
    for i in "MID":
        cigar_string = cigar_string.replace(i, i+" ")
    cigar_list = cigar_string.strip("\n").split(" ")
    # Initialise the return dictionary:
    answer = {"M": 0, "I": 0, "D": 0}
    # Count and store results inside `answer`
    for i in cigar_list:
        if "M" in i:
            answer["M"] += int(i[:-1])
        elif "D" in i:
            answer["D"] += int(i[:-1])
        elif "I" in i:
            answer["I"] += int(i[:-1])
    # Take into account the amount of gaps instead of the length of gaps (i.e.
    # compress gaps into a length of one).
    answer["compressed"] = cigar_string.count("D") + cigar_string.count("I")

    return answer

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
    "matches" + "mismatches_and_gaps" == "ali_len" == "cigM" + "cigI" + "cigD"
        or
    "mismatches" + "matches" == "cigM"
    """
    def __init__(self, path_to_paf, pattern_scaff="scaff", pattern_chr="chr"):
        """
        """
        # Make sure that the provided path points to an existing file
        try:
            with open(path_to_paf):
                pass
            Printing("Opening the file "+str(path_to_paf)).status()
        except:
            Printing("The provided path does not point to an existing "+
                     "file").error()
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
        Printing("Reading the 12 first standard columns").status()
        df = pd.read_table(path_to_paf, header=None,
            # The regex "\s+" stands for "one or more blank spaces" (includes "\t")
            sep="\s+",
            # Load the appropiate (always present) columns and no more.
            usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
            # Use the following column names (from first to last column).
            names=colnames)

        # Create lists with the optional fields (which may not be included in
        # all rows). See minimap2 reference manual for an explanation of each
        # label/column/"tag":  <https://lh3.github.io/minimap2/minimap2.html#10>
        # Keep in mind that the keys or "tags" have been manually sorted by
        # importance (so in the dataframe, important columns appear first).
        tag_to_columns = {"tp": [], "NM": [], "nn": [], "dv": [], "de": [],
                          "cg": [], "cs": [], "SA": [], "ts": [],
                          "cm": [], "s1": [], "s2": [], "MD": [], "AS": [],
                          "ms": [], "rl": [], "zd": [], }
        # Furthermore, analyse the CIGAR string (if any)
        cigar = {"cig_deletions": [], "cig_insertions": [], "cig_matches": [],
                 "cig_compressed": []}
        # Scroll through the lines in the file and populate these lists:
        Printing("Reading the optional columns (further than 12th)").status()
        with open(path_to_paf) as file:
            for count, line in enumerate(file):
                # Get a list with all the fields as items.
                line = line.strip("\n").split("\t")[12:]
                while line:
                    field = line.pop(0)
                    tag = field[0:2]  # "tags" are the first two characters
                    value = field[5:] # values are from the fifth to the end
                    tag_to_columns[tag].append(value)
                # Fill with "None" the fields that were not found in this line:
                for key, val in tag_to_columns.items():
                    if len(val) < count+1:
                        tag_to_columns[key].append(None)
                # If this line contained a CIGAR string, do its analysis:
                if tag_to_columns["cg"][-1] != None:
                    cigar_summary = cigar_analysis(tag_to_columns["cg"][-1])
                    cigar["cig_deletions"].append(cigar_summary["D"])
                    cigar["cig_insertions"].append(cigar_summary["I"])
                    cigar["cig_matches"].append(cigar_summary["M"])
                    cigar["cig_compressed"].append(cigar_summary["compressed"])
                else:
                    cigar["cig_deletions"].append(None)
                    cigar["cig_insertions"].append(None)
                    cigar["cig_matches"].append(None)
                    cigar["cig_compressed"].append(None)

        # Translate these small tags into more descriptive column names:
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
        # Finally, add optional fields to the dataframe.
        for tag, values_list in tag_to_columns.items():
            df[tag_to_colnames[tag]] = values_list
        # Add the sum of lengths of the CIGAR string:
        for colname, values_list in cigar.items():
            df[colname] = values_list

        # Change these optional columns to numerical type, if it is pertinent:
        to_numerical = ["num_minimizers", "chaining_score",
                        "second_chain_score", "mismatches_and_gaps",
                        "DP_ali_score", "DP_max_score", "ambiguous_bases",
                        "divergence", "gap_compr_diverg",
                        "length_repetitive_seeds", "zd_unknown",
                        "cig_deletions", "cig_insertions", "cig_matches"]
        for column in to_numerical:
            df[column] = pd.to_numeric(df[column])
        # Drop columns with zero non-null entries (drop empty columns)
        df.dropna(axis="columns", how="all", inplace=True)

        # Now that the columns have been defined as integer/floats, compute some
        # more parameters and add them to the dataframe:
        df["blast_identity"] = df["matches"] / df["ali_len"]
        df["gap_compr_identity"] = (
            df["matches"] / (df["cig_matches"] + df["cig_compressed"]))
        df["mismatches"] = (
            df["mismatches_and_gaps"] - (df["cig_insertions"]+df["cig_deletions"]))

        # Identify chromosomes and minor scaffolds.
        Printing("Identifying sequence IDs belonging to scaffolds and "+
                     "chromosomes").status()
        nrows_total = df.shape[0]
        for i in ("Qname", "Tname"):
            print(" -- Matches in the", i, "column --")
            # Start by finding rows not matching (tilde is not) the "chr" pattern
            df_copy = df.loc[~ df[i].str.contains(pattern_chr, case=False)]
            nrows_matching_chr = nrows_total - df_copy.shape[0]
            print("Rows matching the chr pattern:", nrows_matching_chr,
                  f"[{round((nrows_matching_chr/nrows_total)*100)} %]")
            df_copy = df_copy.loc[~ df[i].str.contains(pattern_scaff, case=False)]
            nrows_matching_scaff = (nrows_total - df_copy.shape[0] -
                nrows_matching_chr)
            print("Rows matching the scaffold pattern:", nrows_matching_scaff,
                  f"[{round((nrows_matching_scaff/nrows_total)*100)} %]")
            print("Rows not matching any pattern:", df_copy.shape[0])
        # Save these types (sequences are either a minor scaffold or a
        # chromosome)
        ###### Maybe storing sequence type in a column is not needed because
        ###### you can call filter the dataframe at any moment
        ###### So, save these patterns in a variable?
        self.pattern_chr = pattern_chr
        self.pattern_scaff = pattern_scaff

        # Add "Q." and "T." labels to the beginning of Qname and Tname columns
        df["Qname"] = pd.Series(["Q"]*df.shape[0]).str.cat(df["Qname"], sep=".")
        df["Tname"] = pd.Series(["T"]*df.shape[0]).str.cat(df["Tname"], sep=".")

        # Storing modified pandas dataframe in a class variable:
        Printing("The dataframe that has been read has the following "+
                     "columns...").status()
        df.info()  # Print information about the columns that have been read.
        Printing("Printing the first five rows of the dataframe...").status()
        print(df.head(5))
        Printing("Storing modified pandas dataframe in a variable "+
                     "(self.df)").status()
        self.df = df
        # Create an instance of the previous dataframe which only stores
        # mappings including at least one chromosome in columns "Qname" or
        # "Tname"...
        Printing("Storing modified pandas dataframe filtered by chromosome "+
                 "pattern in a variable (self.df_chromosomes and "+
                 "self.df_both_chromosomes)").status()
        self.df_chromosomes = df.loc[
            df["Qname"].str.contains(self.pattern_chr, case=False) |
            df["Tname"].str.contains(self.pattern_chr, case=False)]
        self.df_both_chromosomes = df.loc[
            df["Qname"].str.contains(self.pattern_chr, case=False) &
            df["Tname"].str.contains(self.pattern_chr, case=False)]

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

    def list_chromosomes(self,):
        """
        Prints a list with all the chromosomes for the query and for the target.

        Output
        ======

        Returns a dictionary with two keys: "query" and "target". Their values
        are two lists with their corresponding unique chromosomes.
        """
        results = {}
        for col in (("Qname", "query"), ("Tname", "target")):
            Printing(f"Printing the chromosome IDs for the {col[1]} ("+
                     f"{col[0]} column)").status()
            df = self.df.loc[self.df[col[0]].str.contains(self.pattern_chr,
                                                         case=False)]
            sequences = humansorted(df[col[0]].unique())
            results[col[1]] = sequences

            for sequid in sequences:
                print(f"+ {sequid}")

        return results

    def list_sequences(self,):
        """
        Prints a list with all the sequences for the query and for the target
        (i.e., all chromosomes and minor scaffolds).

        Output
        ======

        Returns a dictionary with two keys: "query" and "target". Their values
        are two lists with their corresponding unique sequences.
        """
        results = {}
        for col in (("Qname", "query"), ("Tname", "target")):
            Printing(f"Printing the sequence IDs for the {col[1]} ("+
                     f"{col[0]} column)").status()
            df = self.df
            sequences = humansorted(df[col[0]].unique())
            results[col[1]] = sequences

            for sequid in sequences:
                print(f"+ {sequid}")

        return results

    def chromosome_coverage(self, sequid):
        """
        Computes the mapping coverage of a single chromosome.

        The mapping coverage consists of the sum of unique queried/targeted
        bases divided by the total length (in bases) of the chromosome. It takes
        overlapping mappings into account, counting them only once.

        Input
        =====

        + sequid: A string pointing to a sequence ID (either query or target).
          Make sure to include the initial "Q." or "T.". The list of beginning
          and ending coordinates of each mapping will be inferred. The total
          length of the sequence is also inferred from the sequid.

        Output
        ======

        Returns a dictionary with the mapping coverage of the selected
        chromosome (in basepairs and relative to the chromosome's total length.
        """
        # Find out whether the sequence is from the query, target or mispelled.
        if sequid[0] == "Q":
            columns = ("Qname", "Qstart", "Qend", "Qlen")
        elif sequid[0] == "T":
            columns = ("Tname", "Tstart", "Tend", "Tlen")
        else:
            Printing("The given sequence ID does not start with neither 'Q' "+
                     "nor 'T'").error()
            return None

        # Infer the interval list from the given sequid.
        df = self.df.loc[self.df[columns[0]] == sequid]
        if df.empty:
            Printing("No sequence exists with the given ID within "+
                     f"the column '{columns[0]}'").error()
            return None
        # Create two lists from beginning and ending coordinates. Items will be
        # paired by index (first in beg. corresponds to first in end).
        begin_list = list(df[columns[1]])
        end_list = list(df[columns[2]])
        intervals = list(zip(begin_list, end_list))
        # Get the sequid's length
        sequid_length = int(df[columns[3]].iloc[0])

        # Create a "point list" from an "intervals list"
        points = []
        for i in range(0, len(intervals)):
            # Left/starting points labeled "0"
            points += [[intervals[i][0], 0, i]]
            # Right/ending points labeled "1"
            points += [[intervals[i][1], 1, i]]
        # Sort points by position (by first value in sublists)
        points.sort()

        # Init variables
        # Keep track of open and closed intervals
        currentOpen = -1
        total_mapped_basepairs = 0

        # For each point in the list
        for i in range(0, len(points)):

            # If the loop landed on a left point (0) which opens an interval:
            if points[i][1] == 0:
                # And there is no other interval opened:
                if currentOpen == -1:
                    # Enter interval "i"
                    currentOpen = points[i][2]
                    currentBegin = int(points[i][0])
                    # print("enters", currentOpen) # debug
                # Otherwise, an open interval already exists:
                else:
                    # From the two mappings that are open, find which has the higher
                    # end coordinate and make it the currentOpen interval.
                    currentEnd = intervals[currentOpen][1]
                    nextEnd = intervals[points[i][2]][1]
                    if nextEnd > currentEnd:
                        currentOpen = points[i][2]

            # If the loop landed on a right point (1)
            # which closes an interval:
            else:
                # And it is the right point of the "currentOpen" interval...
                if points[i][2] == currentOpen:
                    # Compute between begin and end
                    end = points[i][0]
                    mapping_length = end - currentBegin
                    total_mapped_basepairs += int(mapping_length)
                    # Close currentOpen interval
                    currentOpen = -1

        # Compute the mapping coverage (percent of bases covered by at least one
        # mapping)
        percent_mapped = round((total_mapped_basepairs/sequid_length)*100,2)

        Printing(f"The total amount of mapped basepairs for {sequid} was "+
                 f"{total_mapped_basepairs} basepairs [{percent_mapped} "+
                 "% of chr. length]").status()
        return {"mapped_bp": total_mapped_basepairs,
                "mapped_perc": percent_mapped}

    def genome_coverage(self,):
        """
        Computes the mapping coverage of a all the chromosomes in both query and
        target genomes.

        The mapping coverage consists of the sum of unique queried/targeted
        bases divided by the total length (in bases) of the chromosome. It takes
        overlapping mappings into account, counting them only once.

        Input
        =====

        Nothing, reads from `self.df`.

        Output
        ======

        Returns a dictionary. Its keys are the IDs of all chromosomes. Their
        values are nested dictionaries with the total amount of covered
        basepairs and the relative percentage of covered chromosome.
        """
        # Iterate through ALL the chromosomes, in an orderly manner.
        # Exclude minor scaffolds and use `humansorted`.
        sequences = []
        for columns in ["Qname", "Tname"]:
            sequences += list(self.df.loc[self.df[columns].str.contains(
                self.pattern_chr, case=False), columns].unique())
        sequences = humansorted(sequences)

        results = {}
        for seq in sequences:
            results[seq] = self.chromosome_coverage(seq)

        return results

    def type_alignments(self,):
        """
        Prints the number of rows of each alignment type in the PAF file.

        Input
        =====

        Nothing, reads from `self.df`

        Output
        ======

        Returns a dictionary. Its keys are alignment types. Values are tuples
        containing a pair of items: number (0) and decimal percentage (1) of
        rows of each alignment type.
        """
        Printing("Searching for the amount of rows with each type of "+
                 "alignment / mapping...").status()
        df = self.df  # shorten variable
        # Find the total amount of rows in the DataFrame
        n_rows = int(df.shape[0])
        Printing("The total amount of rows in the DataFrame is "+
                 str(n_rows)).status()
        ## I am unsure these are all of the possibilities... Better read the
        ## unique values of the `type_aln` column.
        ## possible_types = {"P": "primary", "S": "secondary", "I": "inversion",
        ##                   "i": "inversion"}
        results = {}
        for value in df["type_aln"].unique():
            number = int(df.loc[df["type_aln"] == value].shape[0])
            percent = round((number/n_rows)*100, 2)
            results[value] = [number, percent]
            Printing("The amount of rows for alignment type '"+
                     f"{value}' is {number} [{percent} %]").status()

        return results

    def boxplots(self, numerical_column, xlabel, out_img_path=None,
                 figure_size=(10, 6), df=pd.DataFrame()):
        """
        IMPROVEMENTS:
        X-axis logarithmic for deletion and insertion sizes. Many outliers make
        it difficult to read the boxplot; a logarithmic transformation would be
        better than cutting off these outliers.

        Create a boxplot from a numerical column. Uses matplotlib and seaborn.

        Input
        =====

        + numerical_column: A numerical column (float or int) of the dataframe.

        + xlabel: Label of the X-axis (units, eg. "basepairs" for alignment
          length).

        + out_img_path [optional]: Path where the figure/image/boxplot will be
          written to. If left unspecified it will show the figure through
          `plt.show()`.

        + figure_size [optional]: A tuple of length two (two numbers), width
          and height (in inches). By default is (10, 6).

        + df [optional]: A filtered dataframe generated from the original; eg.
          with only scaffolds or only chromosomes.

        Output
        ======

        Can write to `out_img_path` a PNG figure.
        """
        # Create one axis for "Tname" column and another axis for "Qname"
        fig, (ax1, ax2) = plt.subplots(
            nrows=1, ncols=2, figsize=figure_size)
        # If df is empty (default option) use the full dataframe stored in
        # self.df (make sure it is a copy of the original dataframe).
        if df.empty:
            df = self.df.copy()
        else:
            df = df.copy()
        # Generate plots with seaborn. Write to the previously created axes.
        for i in (
                ("Qname", ax1, "Q."),
                ("Tname", ax2, "T.")):
            # Change sequence name column to distinguish minor scaffolds and
            # chromosomes. Group all minor scaffold names within "Scaffolds".
            df.loc[df[i[0]].str.contains(self.pattern_scaff, case=False),
                   i[0]] = str(i[2]) + "Scaffolds"
            df = df.sort_values(by=i[0],
                                key=lambda x: np.argsort(index_natsorted(
                                    df[i[0]].str[2:])))
            sns.boxplot(ax=i[1], data=df, x=numerical_column, y=i[0],
                    #showmeans=True, meanprops=dict(color="green", marker="*",
                    #                               markeredgecolor="black",
                    #                               markersize="15")
                    )
            # Set xlabel to both axes
            i[1].set_xlabel(xlabel)

        # Finally, show/save figure in a tight layout.
        plt.tight_layout()
        if out_img_path:
            plt.savefig(out_img_path, dpi=300)
        else:
            plt.show()

    def histogram(self, column, xlabel, out_img_path=None,
                 figure_size=(10, 6), df=pd.DataFrame(), xlogscale=False,
                  ylogscale=False):
        """
        Creates an histogram from almost all columns. Uses matplotlib and seaborn.

        Input
        =====

        + column: A column (float, int or even string) of the
          dataframe. This parameter accepts columns of strings, eg. "Qname",
          "strand", etc. Do not mix a string column with xlogscale parameter
          (logarithmic scale of the X axis).

        + xlabel: Label of the X-axis (units, eg. "basepairs" for alignment
          length).

        + out_img_path [optional]: Path where the figure/image/boxplot will be
          written to. If left unspecified it will show the figure through
          `plt.show()`.

        + figure_size [optional]: A tuple of length two (two numbers), width
          and height (in inches). By default is (10, 6).

        + df [optional]: A filtered dataframe generated from the original; eg.
          with only scaffolds or only chromosomes.

        + xlogscale [optional]: By default, "False". If "True", changes the X
          axis to a logarithmic scale.

        + ylogscale [optional]: By default, "False". If "True", changes the Y
          axis to a logarithmic scale.

        Output
        ======

        Can write to `out_img_path` a PNG figure.
        """
        # If df is empty (default option) use the full dataframe stored in
        # self.df (make sure it is a copy of the original dataframe).
        if df.empty:
            df = self.df.copy()
        else:
            df = df.copy()
        # Agglutinate Qname and Tname: change all scaffolds to the same sequid.
        # This paragraph is useful if the columns "Qname" and "Tname" are given
        # instead of a numerical column.......
        for i in (
                ("Qname", "Q."),
                ("Tname", "T.")):
            # Change sequence name column to distinguish minor scaffolds and
            # chromosomes. Group all minor scaffold names within "Scaffolds".
            df.loc[df[i[0]].str.contains(self.pattern_scaff, case=False),
                   i[0]] = str(i[1]) + "Scaffolds"
            df = df.sort_values(by=i[0],
                                key=lambda x: np.argsort(index_natsorted(
                                    df[i[0]].str[2:])))
        # Create figure
        fig, ax = plt.subplots(figsize=figure_size)
        sns.histplot(data=df, ax=ax, x=column, element="step")
        # Set logscale
        if xlogscale:
            ax.set_xscale("log")
        if ylogscale:
            ax.set_yscale("log")
        # Set the label of the X axis.
        ax.set_xlabel(xlabel)
        # Finally, show/save figure in a tight layout.
        plt.tight_layout()
        if out_img_path:
            plt.savefig(out_img_path, dpi=300)
        else:
            plt.show()


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

    # IMPROVEMENTS
    # Read an array of files. Create <Mapping> instances for each given file in
    # the argument list. Export figures and tables for the suite of PAF
    # mappings, making comparisons between themselves.



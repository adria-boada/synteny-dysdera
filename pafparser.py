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

    + pattern_scaff [optional]: A regex pattern which should match all minor
      scaffolds. To introduce multiple patterns, use the OR operand (which is
      the pipe character '|' for pandas). For instance, "scaff|ctg" will include
      minor scaffolds labelled either "Scaffold_33" and "ctg_V123_456".

    + pattern_chr [optional]: A regex pattern which should match all
      chromosomes.

    + memory_efficient [optional]: Boolean. If True, CIGAR and difference strings will be
      deleted from memory in order to free it up. These strings can be
      expansive, and will end up killing Python processes.

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
    def __init__(self, path_to_paf, pattern_scaff="scaff", pattern_chr="chr",
                 memory_efficient=True):
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
        # Store the file name for later use in tables and figures. Remove the
        # path to the file (remove directories between '/').
        self.filename = path_to_paf.strip("\n").split("/")[-1]
        # Remove the extension (.paf). Lower() returns lowercased string.
        l = self.filename.split(".")
        if l[-1].lower() == "paf":
            self.filename = ''.join(l[:-1])

        # As stated in the manual page, "PAF may optionally have additional
        # fields ...", i.e. an irregular number of columns (ragged TSV).
        # Try to deduce columns.
        # First, read the columns that are always present (until the 12th)
        colnames = ["Qname", "Qlen", "Qstart", "Qend", "strand",
                   "Tname", "Tlen", "Tstart", "Tend",
                   "matches",  # Number of matching bases
                   "ali_len",  # Sum of matches, mismatches and gaps
                   "mapQ"]     # Mapping quality
        # Create a pandas DataFrame by reading the PAF file in `path_to_paf`
        Printing("Reading the 12 first standard columns").status()
        df = pd.read_table(path_to_paf, header=None,
            # The regex "\s+" stands for "one or more blank spaces" (includes "\t")
            sep="\s+",
            # Load the appropiate (always present) columns and no more.
            usecols=list(range(12)), # [0, 1, 2, etc, 10, 11]
            # Use the following column names (from first to last column).
            names=colnames)

        # Create lists with the optional fields (which may not be included in
        # all rows). See minimap2 reference manual for an explanation of each
        # label/column/"tag":  <https://lh3.github.io/minimap2/minimap2.html#10>
        # Keep in mind that the keys or "tags" have been manually sorted by
        # importance (so, in the dataframe, important columns appear first).
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
        # Finally, add the previous optional fields to the dataframe.
        for tag, values_list in tag_to_columns.items():
            df[tag_to_colnames[tag]] = values_list
        # Add the sum of lengths of the CIGAR string:
        for colname, values_list in cigar.items():
            df[colname] = values_list

        # Change these optional columns to numerical type, if it is pertinent.
        # Try to use few columns and efficient dtypes by downcasting if possible
        # (alleviate memory usage; eg, instead of float64, try to use float32
        # dtype when possible).
        to_numerical_ints = ["Qlen", "Qstart", "Qend", "Tlen", "Tstart", "Tend",
                             "matches", "ali_len", "mapQ",
                             "mismatches_and_gaps", "ambiguous_bases",
                             "num_minimizers", "chaining_score", "DP_ali_score",
                             "DP_max_score", "length_repetitive_seeds",
                             "cig_deletions", "cig_insertions", "cig_matches",
                             "cig_compressed"]
        # Apply downcasted integer to the previous columns.
        df[to_numerical_ints] = df[to_numerical_ints].apply(pd.to_numeric,
                                                    downcast="integer")
        to_numerical_floats = ["gap_compr_diverg", "second_chain_score",
                               "zd_unknown"]
        # Apply downcasted float to the previous columns.
        df[to_numerical_floats] = df[to_numerical_floats].apply(pd.to_numeric,
                                                    downcast="float")
        # Drop columns with zero non-null entries (drop empty columns)
        df.dropna(axis="columns", how="all", inplace=True)
        # Another great way to efficiently leverage the computer's memory is to
        # use 'category' dtypes instead of text/object dtypes.
        to_categorical = ["Qname", "Tname", "strand", "type_aln"]
        df[to_categorical] = df[to_categorical].astype("category")

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
            # If some rows do not match any given pattern, WARN the user.
            if df_copy.shape[0] != 0:
                Printing(f"Some rows from the {i} column have not "+
                         "matched neither chromosome "+
                         "nor scaffold pattern. Revise the parameters "+
                         "`pattern_scaff` and `pattern_chr` when calling "+
                         "the `Mapping()` class!").warning()

        # Store this categories (scaffold and chromosome) in a categorical
        # column. First, init. a column with all rows equal to the string
        # "Qscf_Tscf", symbolising that both query and target sequids are minor
        # scaffolds:
        df["sequid_types"] = "Qscf_Tscf"
        # Then, fill with other categories:
        df.loc[df["Qname"].str.contains(pattern_chr, case=False),
               "sequid_types"] = "Qchr_Tscf"
        df.loc[df["Tname"].str.contains(pattern_chr, case=False),
               "sequid_types"] = "Qscf_Tchr"
        # Both Qname and Tname columns contain a chromosome...
        df.loc[
            df["Qname"].str.contains(pattern_chr, case=False) &
            df["Tname"].str.contains(pattern_chr, case=False),
            "sequid_types"] = "Qchr_Tchr"
        # Change the dtype of column "sequid_types" to categorical:
        df["sequid_types"] = df["sequid_types"].astype("category")
        # Save the patterns that have been used for sequid type identification:
        self.pattern_chr = pattern_chr
        self.pattern_scaff = pattern_scaff

        # Add "Q." and "T." labels to the beginning of Qname and Tname columns
        df["Qname"] = pd.Series(["Q"]*df.shape[0]).str.cat(df["Qname"], sep=".")
        df["Tname"] = pd.Series(["T"]*df.shape[0]).str.cat(df["Tname"], sep=".")

        # Storing modified pandas dataframe in a class variable:
        Printing("Printing the first five rows of the dataframe...").status()
        print(df.head(5))
        Printing("Storing modified pandas dataframe in a variable "+
                     "(self.df)").status()
        self.df = df

        # Create a dictionary where sequence sizes are stored...
        Printing("Creating a dictionary to store sequence sizes: "+
                 "`self.sequence_sizes`").status()
        queries, targets = (df["Qname"].unique(), df["Tname"].unique())
        self.sequence_sizes = dict()
        for seq in queries:
            size = int(df.loc[df["Qname"]==seq, "Qlen"].head(1))
            self.sequence_sizes[seq] = size
        for seq in targets:
            size = int(df.loc[df["Tname"]==seq, "Tlen"].head(1))
            self.sequence_sizes[seq] = size

        # From the following lists we can obtain means, medians and st. devs.
        Printing("Creating multiple lists of coordinates, distances between "+
                 "mappings and mapping lengths: `self.coordinates`").status()
        self.coordinates = self._list_coordinates()
        # By default, to compute the size of indels we include the whole main
        # dataframe (ie include indels in minor scaffolds!).
        Printing("Computing a dataframe with the size of indels: "+
                 "`self.df_indels`").status()
        self.df_indels = self._list_indels()
        Printing("Finished computing `self.df_indels`. Printing its "+
                 "five first rows:").status()
        print(self.df_indels.head(5))

        if memory_efficient:
            Printing("Because `memory_efficient` setting is set to True, "+
                     "CIGAR and difference strings will be deleted from "+
                     "the DataFrame to save up memory").status()
            self.df = self.df.drop(columns=["cigar_string", "diff_string"])

        Printing("The dataframe that has been read has the following "+
                     "columns...").status()
        self.df.info(memory_usage="deep")

        Printing("Finished initialising "+self.filename+"!\n").status()

    def _list_coordinates(self, df=pd.DataFrame()):
        """
        Hidden function that produces an intermediate list in order to compute
        mean interdistances (ie. mean distance between individual mappings) and
        mean alignment length.

        Input
        =====

        + df [optional]: A filtered dataframe. By default it will use the rows
          of `self.df` where both query and target sequids are chromosomes.

        Output
        ======

        A dictionary where keys are sequence IDs and values are lists of begin
        and end pairs of coordinates.
        """
        if df.empty:
            # Copy the dataframe where query and target are chromosomes.
            df = self.df.loc[self.df["sequid_types"]=="Qchr_Tchr"]
        else:
            df = df.copy()
        answer = dict()
        # Create a list of sorted coordinates indicating regions that are included
        # in the mapping. Get one list per chromosome.
        # Eg. dict[sequence] = ((begin1, end1), (b2, e2), ..., (bN, eN))
        coordinate_dict = dict()
        queries, targets = (df["Qname"].unique(), df["Tname"].unique())
        for seq in queries:
            d = self.df.loc[self.df["Qname"] == seq]
            l = list(zip(list(d["Qstart"]), list(d["Qend"])))
            l.sort(key=lambda x: x[0])
            coordinate_dict[seq] = l
        for seq in targets:
            d = self.df.loc[self.df["Tname"] == seq]
            l = list(zip(list(d["Tstart"]), list(d["Tend"])))
            l.sort(key=lambda x: x[0])
            coordinate_dict[seq] = l
        answer["coordinates_dict"] = coordinate_dict

        # Compute the distances between begin and end of the previous
        # coordinates list, per chromosome...
        answer["alig_lengths"], answer["interdistances"] = (dict(), dict())
        answer["Q.alig_lengths"], answer["T.alig_lengths"] = (list(), list())
        answer["Q.interdist"], answer["T.interdist"] = (list(), list())

        for seq, values in coordinate_dict.items():
            # Get the first pair before starting the loop.
            last_start, last_end = values.pop(0)
            interdistances = []
            alig_lengths = []
            for coords in values:
                # Last end precedes current begin; ie. mappings are not
                # overlapping.
                if coords[0] - last_end > 0:
                    # Distance between mappings.
                    interdistances.append(coords[0] - last_end)
                    # Alignment length of current mapping (end minus begin).
                    alig_lengths.append(coords[1] - coords[0])
                    last_start = coords[0]
                last_end = coords[1]
            if seq[0] == "Q":
                answer["Q.alig_lengths"].extend(alig_lengths)
                answer["Q.interdist"].extend(interdistances)
            else:
                answer["T.alig_lengths"].extend(alig_lengths)
                answer["T.interdist"].extend(interdistances)
            # Finally, add per sequence lists...
            answer["alig_lengths"][seq] = alig_lengths
            answer["interdistances"][seq] = interdistances

        return answer

    def _list_indels(self, df=pd.DataFrame()):
        """
        """
        # If df is empty (default option) use the full dataframe stored in
        # self.df (make sure it is a copy of the original dataframe).
        if df.empty:
            df = self.df.copy()
        else:
            df = df.copy()
        # Start by creating a dataframe with all gaps in the alignment:
        # Initialise "answer" dictionary (afterwards will generate a df)
        answer = {"Qname": [], "Tname": [],
                  "gap_type": [], "gap_size": []}
        # Iterate across `df` rows and analyse CIGAR strings.
        for index, cigar_string in df["cigar_string"].items():
            # Split string into a list by separating letters with a space:
            for i in "MID":
                cigar_string = cigar_string.replace(i, i+" ")
            cigar_list = cigar_string.strip("\n").split(" ")
            # Remove matches from the list (we're interested in indels).
            # Furthermore, remove the ending "D" or "I" and `intize`.
            deletions = (
                [int(indel[:-1]) for indel in cigar_list if "D" in indel])
            insertions = (
                [int(indel[:-1]) for indel in cigar_list if "I" in indel])
            qname, tname = df.loc[index, ["Qname", "Tname"]]
            answer["Qname"].extend([qname]*len(insertions + deletions))
            answer["Tname"].extend([tname]*len(insertions + deletions))
            answer["gap_type"].extend(["Del"]*len(deletions))
            answer["gap_size"].extend(deletions)
            answer["gap_type"].extend(["Ins"]*len(insertions))
            answer["gap_size"].extend(insertions)

        # Translate the dictionary into a pandas DataFrame.
        answer = pd.DataFrame(answer)
        # Use efficient dtypes, like "categorical" for columns with repetitive
        # string data.
        to_categorical = ["Qname", "Tname", "gap_type"]
        answer[to_categorical] = answer[to_categorical].apply(pd.Categorical)
        answer["gap_size"] = pd.to_numeric(answer["gap_size"],
                                           downcast="integer")
        return pd.DataFrame(answer)

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
            df = self.df.copy()
            sequences = humansorted(df[col[0]].unique())
            results[col[1]] = sequences

            for sequid in sequences:
                print(f"+ {sequid}")

        return results

    def coverage_sequid(self, sequid):
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
        sequid_length = int(self.sequence_sizes[sequid])

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
                "mapped_perc": percent_mapped,
                "sequid_length_bp": sequid_length}

    def coverage_all_chromosomes(self,):
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
            results[seq] = self.coverage_sequid(seq)

        return results

    def coverage_genome(self, df=pd.DataFrame()):
        """
        Computes the mapping coverage of both query's and target's genomes.
        Instead of giving the coverage per chromosome, gives the coverage for
        the whole genome, calculated by the sum of basepairs of all sequids.

        To acquire the genome coverage discounting the minor scaffolds, filter
        them out from the dataframe and input it as the variable `df` (see
        pandas `loc` method to filter a dataframe).

        Input
        =====

        + df [optional]: A filtered dataframe. By default reads from the main
          dataframe, ie. `self.df`.

        Output
        ======

        A dictionary with the following values: (1) the sum of uniquely covered
        base pairs in both query and target sequids, (2) the sum of base pairs
        in the whole genome of query and target sequids, (3) and the genome
        coverage for both query and target species.

        For instance:
        + returned_dict["query"]["genome_coverage"]
        + returned_dict["query"]["covered_bp"]
        + returned_dict["query"]["genome_bp"]
        """
        if df.empty:
            df = self.df.copy()
        else:
            # Just in case, create a copy so the original `df` is unmodified.
            df = df.copy()
        # Initialise answer dictionary
        Printing("Computing genome coverage for "+self.filename).status()
        answer = {
            "query":{"covered_bp": 0, "genome_bp": 0},
            "target":{"covered_bp": 0, "genome_bp": 0}}
        # Iterate through ALL the sequids, including minor scaffolds. Read
        # function help for filtering them out.
        query, target = (list(df["Qname"].unique()), list(df["Tname"].unique()))
        for seq in query:
            cov = self.coverage_sequid(seq)
            answer["query"]["covered_bp"] += cov["mapped_bp"]
            answer["query"]["genome_bp"] += cov["sequid_length_bp"]
        for seq in target:
            cov = self.coverage_sequid(seq)
            answer["target"]["covered_bp"] += cov["mapped_bp"]
            answer["target"]["genome_bp"] += cov["sequid_length_bp"]

        # Compute the division of: (uniquely covered bp / all bp)
        for i in ("query", "target"):
            j = answer[i] # short-hand
            j["genome_coverage"] = (j["covered_bp"] /j["genome_bp"])

        return answer

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
        df = self.df.copy()  # shorten variable
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

    def small_indels_histogram(self, stop:int=16, step:int=5, out_img_path=None,
                               df=pd.DataFrame(), ylogscale:bool=False):
        """
        Create a histogram counting the amount of events per indel category,
        resembling the histogram of Petrov 2002 (Mutational Equilibrium Model of
        Genome Size Evolution).

        Indel categories start from one and increase by `step` between them
        (e.g. with step=5, categories are 1-5, 6-10, 11-15, etc.). For
        categories bigger than `stop` in size, they are all included in a single
        bin for indels beyond this specific size (i.e. bigger than `stop` size).

        Input
        =====

        + stop [optional]: Integer. Compute indel categories up to this length, in
          basepairs. If there are indels larger than this integer, all of them
          will be binned in "equal or higher than this integer".

        + step [optional]: Integer. Width of categories or bins, in basepairs.

        + out_img_path [optional]: Path where the figure will be written to. If
          left unspecified it will show the figure through `plt.show()`.

        + df [optional]: A filtered dataframe generated from the original; eg.
          with only scaffolds or only chromosomes.

        + ylogscale [optional]: By default, "False". If "True", changes the Y
          axis to a logarithmic scale.

        Output
        ======

        An indel histogram figure, resembling the previously cited publication
        of Petrov 2002.
        """
        # The dataframe with small indels should've been computed in __init__
        df = self.df_indels
        # Create one axis for a relative plot and another for an absolute count.
        # Moreover, include weighted and non-weighted plots (x and x^2 plots).
        ## cm = 1/2.54 # centimeters to inches conversion factor
        fig, axes = plt.subplots(
            nrows=2, ncols=2, figsize=(10,6))
        # Absolute count, non-weighted...
        sns.histplot(ax=axes[0,0], data=df, x="gap_size", hue="gap_type",
                     multiple="dodge", binwidth=step, binrange=(1,stop),
                     stat="count", common_norm=False, )
        # Relative without common normalisation of "Del" and "Ins",
        # non-weighted...
        sns.histplot(ax=axes[0,1], data=df, x="gap_size", hue="gap_type",
                     multiple="dodge", binwidth=step, binrange=(1,stop),
                     stat="proportion", common_norm=False, )
        # Absolute count, weighted...
        sns.histplot(ax=axes[1,0], data=df, x="gap_size", hue="gap_type",
                     weights="gap_size", multiple="dodge", binwidth=step,
                     binrange=(1,stop), stat="count", common_norm=False, )
        # Relative without common normalisation of "Del" and "Ins", weighted...
        sns.histplot(ax=axes[1,1], data=df, x="gap_size", weights="gap_size",
                     hue="gap_type", multiple="dodge", binwidth=step,
                     binrange=(1,stop), stat="proportion", common_norm=False, )
        if ylogscale:
            [ax.set_yscale("log") for row in axes for ax in row]
        axes[0,0].set_title("Abs. and non-weighted")
        axes[0,1].set_title("Rel. and non-weighted")
        axes[1,0].set_title("Abs. and weighted")
        axes[1,1].set_title("Rel. and weighted")
        plt.tight_layout()
        # Plot the figure.
        if out_img_path:
            plt.savefig(out_img_path, dpi=300, bbox_inches="tight")
        else:
            plt.show()
        # When finished, close the plot
        plt.close()

def table_coordinates(mapping_list: list, ci_factor=1.96, plots: bool=False):
    """
    Input
    =====

    + mapping_list: A list of `Mapping` objects from multiple, differents files.

    + ci_factor [optional]: The confidence interval factor. By default,
      multiplies standard deviation by `1.96` to obtain the 95 % confidence
      interval.

    + plots [optional]: Create histograms portraying the datasets to make sure
      these are normal distribution or else.

    Output
    ======

    A table with interesting values computed from individual mapping
    coordinates.
    """
    # Init. return columns of table as dictionary
    answer = {"Filename": [], "Map as...": [], "Coverage": [],
              "Interdist.": [], "Align. len.": [], "Del.": [], "Ins.": []}
    for mapping in mapping_list:
        for both in (("Q", "query"), ("T", "target")):
            # Get the list of interdistances computed by the Mapping class and
            # create a pandas Series from them, from which median and stdev can be
            # easily computed.
            interdist = pd.Series(mapping.coordinates[both[0]+".interdist"])
            cell_interdist = (str(int(interdist.median()))+" $\pm$ "+
                str(int(interdist.std() * ci_factor)) )
            # Same for the list of alignment lengths.
            alig_lengths = pd.Series(mapping.coordinates[both[0]+".alig_lengths"])
            cell_alig_lengths = (str(int(alig_lengths.median()))+" $\pm$ "+
                str(int(alig_lengths.std() * ci_factor)) )
            # Retrieve the genome coverage and display it as a percentage, rounded
            # to two decimals.
            coverage = mapping.coverage_genome()[both[1]]["genome_coverage"]
            cell_percent_coverage = round(coverage*100, 2)
            # Retrive median deletion and insertion sizes.
            deletions = mapping.df_indels.loc[
                mapping.df_indels["gap_type"]=="Del"]
            cell_deletions = (str(int(deletions.median()))+" $\pm$ "+
                str(int(deletions.std() * ci_factor)))
            insertions = mapping.df_indels.loc[
                mapping.df_indels["gap_type"]=="Ins"]
            cell_insertions = (str(int(insertions.median()))+" $\pm$ "+
                str(int(insertions.std() * ci_factor)))

            # Create plots to make sure these datasets are normally
            # distributed...
            if plots:
                for data, title in zip((interdist, alig_lengths, deletions,
                                insertions), ('Interdist', 'Alig. len.',
                                'Deletions', 'Insertions')):
                    sns.histplot(data=data, element="step")
                    plt.title(title)
                    plt.show()

            # Add `rows` to answer in order to create a dataframe later on:
            for col, val in zip(answer.keys(), (
                    mapping.filename, both[1].capitalize(), cell_percent_coverage,
                    cell_interdist, cell_alig_lengths, cell_deletions,
                    cell_insertions)):
                answer[col].append(val)

    # Create the table from the `answer` dictionary.
    df = pd.DataFrame(answer)

    print("\nKeep in mind that, by default, the following results "+
          "come from a dataframe where all minor scaffolds had been "+
          "completely removed.")
    print(df.to_markdown(index=False))
    return df

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
    parser.add_argument('--paflist', help="A list of PAF filenames", nargs="+")
    # integer argument
#    parser.add_argument('-a', '--numero_a', type=int, help='Paràmetre "a"')
    # choices argument
#    parser.add_argument('-o', '--operacio', 
#            type=str, choices=['suma', 'resta', 'multiplicacio'],
#            default='suma', required=False,
#            help='')

    args = parser.parse_args()
    # call a value: args.operacio or args.filename.

    # Read an array of files. Create <Mapping> instances for each given file in
    # the argument list. Export figures and tables for the suite of PAF
    # mappings, making comparisons between themselves.
    print(args.paflist)
    maplist = [Mapping(file, pattern_scaff="scaf|ctg") for file in args.paflist]
    # Create small indel histograms.
    for m in maplist:
        m.small_indels_histogram(stop=400, step=5,
            out_img_path="small_indels_"+str(m.filename)+".png")
    table_coordinates(maplist, plots=False)



#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# pafparser.py (v5)
#
# 11 de febr. 2024  <adria@molevol-OptiPlex-9020>

help_msg = """
The purpose of this Python script is to parse (read) a PAF file and create
compelling plots and tables from it. The script defines a class which creates a
more manageable Python3 object from a PAF mapping file. It stores information in
a pandas dataframe. It also contains a function that subsets a given PAF file.
It is useful to create datasets for testing purposes (before proceeding with the
main analyses).
"""

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
    # This function is the only fragment of the script which uses the `random`
    # module.
    import random

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

class Printing(object):
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

class Mapping(object):
    """
    Create a Python object from a PAF file. A documentation of the columns in a
    PAF file can be found in the manual of `minimap2`:

    <https://lh3.github.io/minimap2/minimap2.html#10>

    In summary:
    "matches" + "mismatches_and_gaps" == "ali_len" == "cigM" + "cigI" + "cigD"
        or
    "mismatches" + "matches" == "cigM"

    Input
    -----

    + path_to_paf: Path to a PAF file in order to read it and create a pandas
      DataFrame from it.

    + sequence_patterns: A dictionary with patterns expected to match sequences
      found in both the query and target databases (`Qname` and `Tname`
      columns). Divides the alignment in autosomes, genome, etc. partitions. The
      patterns will be matched case-insensitively. For instance:

      dict = {"genome": ".*",      # `.*` matches any sequence ID.
       "autosomes": "chr.*\\d",    # matches `chr` followed by any digit (\\d).
       "sex_chr": "chr.*x",        # matches `chr` followed by "x" or "X".
       "scaffolds": "Scaf|ctg" }   # matches `scaf` or (|) `ctg` strings.

    + memory_efficient [obsolete]: Boolean. If True, CIGAR and difference
      strings will be deleted from memory in order to free it. These strings
      can be expansive and memory intensive, eventually killing the Python3
      process.

    Output
    ------

    Returns a Mapping() object, hopefully more manageable through pandas/python3
    than the raw data.
    """
    def __init__(
        self,
        path_to_paf: str,
        sequence_patterns: dict={
            "genome": ".*",
            "autosomes": "chr.*\\d",
            "sex_chr": "chr.*x",
            "scaffolds": "Scaffold|HRSCAF|ctg", },
        memory_efficient: bool=True, ):
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
            self.filename = ".".join(l[:-1])
        else:
            self.filename = ".".join(l)

        # Read the input PAF file and create a pandas.DataFrame() from it.
        df = self._read_paf_file(path_to_paf)

        # Identify chromosomes and minor scaffolds.
        self._identify_sequence_patterns(
            df=df, patterns=sequence_patterns)

        mapped_regions, seqtype_lens = self._mapped_regions(
            df=df, patterns=sequence_patterns)

        df_stats = self._stats_mapped_regions(
            mapped_regions=mapped_regions, seqtype_lens=seqtype_lens)

        # Store variables into the class.
        self.df = df
        self.mapped_regions = mapped_regions
        self.df_stats = df_stats

    def cigar_analysis(
        self,
        cigar_string: str, ):
        """
        Analyse a CIGAR string and compute its total amount of "M" (matches), "I"
        (insertions to query) and "D" (deletions to query).

        Input
        -----

        + cigar_string: A string of letters and numbers symbolizing the gap
          composition of an alignment/mapping. Accepts a string with "M", "I" and
          "D". Each letter must be preceded by a number (which is the length of the
          fragment). Take a look at the specifications of the "SAM" format for a
          more in-depth explanation:

          <https://samtools.github.io/hts-specs/SAMv1.pdf>

        Output
        ------

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

    def _read_paf_file(
        self,
        path_to_paf: str, ):
        """
        """
        # As stated in the manual page, "PAF may optionally have additional
        # fields ...", i.e. an irregular number of columns (ragged TSV).
        # Try to deduce columns.
        # First, read the columns that are always present (until the 12th)
        colnames = ["Qname",    # Name of a query scaffold/chr.
                    "Qlen",     # Length of above sequence.
                    "Qstart",   # Start of alignment in Qname.
                    "Qend",     # End of alignment in Qname.
                    "strand",   # Complement DNA seq.
                    # Same variables for target genome/assembly.
                    "Tname", "Tlen", "Tstart", "Tend",
                    "matches",  # Number of matching bases.
                    "ali_len",  # Sum of matches, mismatches and gaps.
                    "mapQ"]     # Mapping quality.

        # Create a pandas DataFrame by reading the PAF file in `path_to_paf`.
        Printing("Reading the 12 first standard columns.").status()
        df = pd.read_table(path_to_paf, header=None,
            # The regex "\s+" stands for "one or more blank spaces" (includes "\t")
            sep="\s+",
            # Load the appropiate (always present) columns and no more.
            usecols=list(range(12)), # [0, 1, 2, etc, 10, 11]
            # Use the following column names (from first to last column).
            names=colnames)

        # Create lists with the optional fields (which may not be included in
        # all rows). See minimap2 reference manual for an explanation of each
        # label/column/"tag":
        # <https://lh3.github.io/minimap2/minimap2.html#10>
        # Keep in mind that the keys or "tags" have been manually sorted by
        # importance (so, in the dataframe, important columns appear first).
        tag_to_columns = {"tp": [], "NM": [], "nn": [], "dv": [], "de": [],
                          "cg": [], "cs": [], "SA": [], "ts": [],
                          "cm": [], "s1": [], "s2": [], "MD": [], "AS": [],
                          "ms": [], "rl": [], "zd": [], }
        # Furthermore, analyse the CIGAR string (if any).
        cigar = {"cig_deletions": [], "cig_insertions": [], "cig_matches": [],
                 "cig_compressed": []}
        # Scroll through the lines in the file and populate these lists:
        Printing("Reading the optional columns (further than 12th).").status()
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
                    last_cigar = tag_to_columns["cg"][-1]
                    cigar_summary = self.cigar_analysis(last_cigar)
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
                          "zd": "zd_unknown", }
        # Finally, add the previous optional fields to the dataframe.
        for tag, values_list in tag_to_columns.items():
            df[tag_to_colnames[tag]] = values_list
        # Add the sum of lengths of the CIGAR string:
        for colname, values_list in cigar.items():
            df[colname] = values_list

        # Make sure the dtypes of these columns is set to integer.
        to_numerical_ints = ["Qlen", "Qstart", "Qend", "Tlen", "Tstart", "Tend",
                             "matches", "ali_len", "mapQ",
                             "mismatches_and_gaps", "ambiguous_bases",
                             "num_minimizers", "chaining_score", "DP_ali_score",
                             "DP_max_score", "length_repetitive_seeds",
                             "cig_deletions", "cig_insertions", "cig_matches",
                             "cig_compressed"]
        df[to_numerical_ints] = df[to_numerical_ints].apply(pd.to_numeric,
                                                    downcast="integer")
        # Make sure the dtypes of these columns is set to float.
        to_numerical_floats = ["gap_compr_diverg", "second_chain_score",
                               "zd_unknown"]
        # Apply downcasted float to the previous columns.
        df[to_numerical_floats] = df[to_numerical_floats].apply(pd.to_numeric,
                                                    downcast="float")

        # Drop empty columns (with zero non-null entries).
        df = df.dropna(axis="columns", how="all")
        # Another great way to efficiently leverage the computer's memory is to
        # use 'category' dtypes instead of text/object dtypes.
        object_columns = df.select_dtypes("object").columns
        df[object_columns] = df[object_columns].astype("category")

        # Now that the columns have been defined as integer/floats, compute some
        # more parameters and add them to the dataframe:
        df["blast_identity"] = df["matches"] / df["ali_len"]
        df["gap_compr_identity"] = (
            df["matches"] / (df["cig_matches"] + df["cig_compressed"]))
        df["mismatches"] = \
            df["mismatches_and_gaps"] - \
            (df["cig_insertions"] + df["cig_deletions"])

        # Add the labels "Q." or "T." to the beginning of Qname and Tname
        # columns.
        df["Qname"] = pd.Series(["Q"]*df.shape[0]).str.cat(df["Qname"], sep=".")
        df["Tname"] = pd.Series(["T"]*df.shape[0]).str.cat(df["Tname"], sep=".")

        return df

    def _identify_sequence_patterns(
        self,
        df: pd.DataFrame(),
        patterns: dict, ):
        """
        """
        Printing("Matching target and query sequence IDs to the given "+
                 "patterns of sequence groups (autosomes, scaffolds, etc.)"
                 ).status()

        nrows_tot = df.shape[0]
        for db_column in ("Qname", "Tname"):
            print(" -- Matches in the", db_column, "column --")
            tot_uniq_seq = len(df[db_column].unique())
            for seq_key, seq_pat in patterns.items():
                # Create a mask with the `seq_pat` pattern in the df column
                # `db_column`.
                mask_seqtype = df[db_column].str.contains(seq_pat, case=False)
                nrows_match_pat = df.loc[mask_seqtype].shape[0]
                percent_pat = round((nrows_match_pat/nrows_tot)*100, ndigits=1)
                unique_sequids = df.loc[mask_seqtype, db_column].unique()
                num_uniq_seq = len(unique_sequids)
                percent_uniq_seq = \
                    round((num_uniq_seq/tot_uniq_seq)*100, ndigits=1)
                # Print nrows matched to screen.
                print("Pattern `"+ str(seq_key)+"`",
                      "matching rows:",
                      str(nrows_match_pat),
                      "["+ str(percent_pat) +" %]")
                # Print num_uniq_seq to the screen.
                print(" "*len("Pattern `"+ str(seq_key)+"`"),
                      "unique sequids:",
                      str(num_uniq_seq),
                      "["+ str(percent_uniq_seq) +" %]")

        return None

    def _list_intervals_in_df(
        self,
        df: pd.DataFrame(),
        column_sequid: str,
        column_start: str, column_end: str, ):
        """
        Input
        -----

        + df: The pd.DataFrame() created in the __init__ function, whole or
          filtered.

        + column_sequid: Either "Qname" or "Tname".

        + column_start: Either "Qstart" or "Tstart" (in conformity with
          `column_sequid` parameter).

        + column_end: Either "Qend" or "Tend" (in conformity with
          `column_sequid` parameter).

        Output
        ------

        A dict. where keys are sequences and values are lists of their mapped
        intervals. Intervals are portrayed by their beginning and ending
        coordinates. For instance:

        values = [ (begin_1, end_1), (b_2, e_2), etc. ]
        """
        # Initialise a dict. to store results and return.
        answer = dict()
        # Initialise a dict. to store a list of sorted intervals indicating
        # regions which are included in the mapping. One list per chromosome.
        # Eg. dict[sequence] = ((begin1, end1), (b2, e2), ..., (bN, eN))
        intervals_per_seq = dict()
        for seq in df[column_sequid].unique():
            df_seq = df.loc[df[column_sequid] == seq]
            intervals = list( zip(
                list(df_seq[column_start]),
                list(df_seq[column_end]), ))
            # Sort the intervals by the starting coordinate.
            intervals = sorted(intervals, key=lambda x: x[0])
            # Store these intervals in the dict.
            intervals_per_seq[seq] = intervals

        return intervals_per_seq

    def _remove_overlapping_mapped_intervals(
        self,
        intervals: list, ):
        """
        """
        # Create a list of points from the list of intervals.
        point_list = []
        for i in range(0, len(intervals)):
            # Left/starting points labeled `0`:
            point_list.append(list( [intervals[i][0], 0, i]))
            # Right/ending points labeled `1`:
            point_list.append(list( [intervals[i][1], 1, i]))
        # Sort the list of points by their position.
        point_list = sorted(point_list, key=lambda x: x[0])
        # Initialise other algorithm variables...
        currentOpen = -1
        answer = list()

        # Start iterating across the list of points. Discount overlapping
        # basepairs of mappings found overlapping.
        for point in point_list:
            # If the loop landed on a left point (0), open an interval.
            if point[1] == 0:

                # If there is no other interval opened:
                if currentOpen == -1:
                    # Then, assign this interval as currentOpen.
                    currentOpen = int(point[2])
                    currentBegin = int(point[0])
                # Else, there already was an open interval.
                else:
                    # From the two intervals that are overlapping, find of them
                    # has the largest end coord. and assign it as `currentOpen`.
                    currentEnd = int(intervals[currentOpen][1])
                    nextEnd = int(intervals[point[2]][1])
                    if nextEnd > currentEnd:
                        currentOpen = int(point[2])

            # Otherwise, the loop landed on a right point (1):
            else:
                # If the right point and the `currentOpen` left point belong to
                # the same interval, close the currently open interval.
                if point[2] == currentOpen:
                    # Store the interval removing the overlapping region.
                    end = point[0]
                    answer.append(list([currentBegin, end]))
                    # Close `currentOpen` interval.
                    currentOpen = -1

        return answer

    def _len_mapped_and_unmapped(
        self,
        intervals: list,
        chr_end: int, chr_beg: int=0, ):
        """
        + intervals: a list of non-overlapping intervals.

        + chr_end: ending coordinate of the chr.

        + chr_beg [optional]: beggining coordinate of the chr. By default set to
          zero.
        """
        # Create a list of points from the list of intervals.
        point_list = []
        for i in range(0, len(intervals)):
            # Extend the list with beg/end.
            point_list.extend(list([
                intervals[i][0],
                intervals[i][1], ]))

        # Init. variables.
        coord_end = chr_beg
        interdist, alig_lens = list(), list()
        # Iterate across the point_list.
        while len(point_list) != 0:
            # Pop the first point (interval beggining).
            coord_beg = point_list.pop(0)
            # Compute distance between last and current alignment.
            interdist.append(int(coord_beg - coord_end))
            # Pop the second point (interval ending).
            coord_end = point_list.pop(0)
            # Compute alignment/mapping length.
            alig_lens.append(int(coord_end - coord_beg))
        # Lastly, compute the ending of last interval to ending of chr.
        interdist.append(int(chr_end - coord_end))

        # Convert lists into pd.Series(); they are easier to manipulate.
        answer = {
            "mapped": pd.Series(alig_lens),
            "unmapped": pd.Series(interdist), }

        return answer

    def _mapped_regions(
        self,
        df: pd.DataFrame(),
        patterns: dict, ):
        """
        """
        # Init a dict. to store `ali_len` and `interdist` (distances between
        # mappings, unmapped regions).
        mapped_regions = dict()
        seqtype_lens = dict()
        # Compute pd.Series with alignment lengths and unmapped lengths.
        for prefix in ("Q", "T"):
            # Init. a dict. for lengths in `prefix`.
            mapped_regions[prefix] = dict()
            seqtype_lens[prefix] = dict()
            # Prepare query or target column names.
            column_sequid, column_start, column_end, column_len = \
                str(prefix) + "name", \
                str(prefix) + "start", \
                str(prefix) + "end", \
                str(prefix) + "len"
            for seqtype, patt in patterns.items():
                # Filter df by the prefix + sequid column using the pattern.
                mask_seqtype = df[column_sequid].str.contains(patt, case=False)
                df_seqtype = df.loc[mask_seqtype].copy()
                # Compute the total length of all sequences in `df_seqtype`
                seq_sum_lens = sum(list(df_seqtype[column_len].unique()))
                # Init. a dict. for lengths in `seqtype`.
                mapped_regions[prefix][seqtype] = {
                    "mapped": pd.Series(dtype="int"),
                    "unmapped": pd.Series(dtype="int"), }
                # Store the sum of sequence lengths for this seqtype.
                seqtype_lens[prefix][seqtype] = seq_sum_lens
                # Create a dict. formatted as {"seq": [intervals]}
                intervals_per_seq = self._list_intervals_in_df(
                    df=df_seqtype, column_sequid=column_sequid,
                    column_start=column_start, column_end=column_end)
                for seq, val in intervals_per_seq.items():
                    # Remove overlapping in the list of intervals `val`.
                    val = self._remove_overlapping_mapped_intervals(
                        intervals=val)
                    # Create two pd.Series of mapped regions' lengths and
                    # unmapped regions' lengths.
                    chr_end = df_seqtype.loc[
                        df_seqtype[column_sequid] == seq, column_len] \
                        .iloc[0]
                    region_series = self._len_mapped_and_unmapped(
                        intervals=val, chr_end=chr_end)
                    answer = mapped_regions[prefix][seqtype]
                    answer["mapped"] = pd.concat(
                        [answer["mapped"], region_series["mapped"]])
                    answer["unmapped"] = pd.concat(
                        [answer["unmapped"], region_series["unmapped"]])

        return mapped_regions, seqtype_lens

    def _stats_mapped_regions(
        self,
        mapped_regions: dict,
        seqtype_lens: dict, ):
        """
        """
        # Compute some stats from the series of alignment lengths or distances
        # between alignments.
        dict_stats = {"db": [], "partition": [], "region": [], "len_median": [],
                      "len_mean": [], "len_max": [], "len_min": [],
                      "bpsum": [], "coverage": [], }
        for db, val1 in mapped_regions.items():
            for partition, val2 in val1.items():
                for region, series in val2.items():
                    seqlens = seqtype_lens[db][partition]
                    dict_stats["db"].append(db)
                    dict_stats["partition"].append(partition)
                    dict_stats["region"].append(region)
                    dict_stats["len_median"].append(series.median())
                    dict_stats["len_mean"].append(series.mean())
                    dict_stats["len_max"].append(series.max())
                    dict_stats["len_min"].append(series.min())
                    dict_stats["bpsum"].append(series.sum())
                    dict_stats["coverage"].append(
                        round((series.sum()/seqlens)*100,ndigits=5))

        df_stats = pd.DataFrame(dict_stats)

        return df_stats

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Si es crida com a script:
if __name__ == '__main__':

    ### PARSE ARGUMENTS ###
    import argparse
    parser = argparse.ArgumentParser(description=help_msg,
        # make sure the 'help_msg' is not automatically
        # wrapped at 80 characters (manually assign newlines).
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--paflist", nargs="+", required=True,
                        help="A list of PAF filenames.")

    args = parser.parse_args()
    # call a value: args.operacio or args.filename.

    # Read an array of files. Create <Mapping> instances for each given file in
    # the argument list. Export figures and tables for the suite of PAF
    # mappings, making comparisons between them.
    maplist = [Mapping(path_to_paf=file) for file in args.paflist]
    for mapping in maplist:
        print(" --", mapping.filename, "--")
        print(mapping.df_stats)

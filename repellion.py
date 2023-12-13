#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# resquitllar.py
#
# 14 Nov 2023  <adria@molevol-OptiPlex-9020>

help_msg = """
In development
==============

Top10 value counts for the "perc_divg" column (and other variable columns, like
replen, if possible). Hand a subset of a dataframe (eg. loc DNA). Otherwise,
hand a list of groupby columns (["class", "Species"]).

Divergence histograms. Try to include max-mins (inflection points using SciPy).

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

```idxfile
Species-label
Sequence-regex1\t32
Sequence-regex2\t5600
etc.
```

Another example:

```idxfile
Species       whatever
Chr1          3600
Chr2          5000
...           ...

ChrN          100
Scaffolds     $(sum_len_scaffolds)
```

Where $(sum_len_scaffolds) is the cumulative sum of scaffold lengths. See the
comments below `class Repeats` for more information on regexes and formatting.

Formatting the RM.out files
===========================

The output obtained from RepeatMasker must be modified
Some modifications must be applied to the output obtained from RepeatMasker.

Some rows contain an asterisk (*) as a 16th field. Some other rows do end at the
15th field. To make sure that all rows are of the same length, one can run:

```bash
awk 'BEGIN{OFS="\t"} ; {if ( $16 ) { $16="True" } else { $16="False" } } 1' RM.out
```

Moreover, make sure to correctly label minor scaffolds and chromosomes in column
5. For instance, to substitute "Scaffold_1" by "DtilchrX" in the sequid field:

```bash
awk '$5=="Scaffold_1" {$5="DtilchrX"}' RM.out
```

Lastly, it is possible for the field "position (left)" to contain asterisks.
These obstruct the loading of the file because the column is otherwise of the
integer dtype/category. One solution is to substitute them by zeroes.
"""

import sys, math
# Reading input files' size and creating folders of PNG plots
import os

# Handling input as DataFrames
import pandas as pd, numpy as np
# Plotting figures
import seaborn as sns, matplotlib.pyplot as plt
# Obtain local minima/maxima values to plot them
from scipy.signal import argrelextrema
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

def prepare_input_idx(list_idx_paths):
    """
    Prepare the IDX input files to use interactively with the class `Repeats`.
    Formats the data according to `Repeats` requirements.

    Input
    =====

    + list_idx_paths: A list with one or multiple paths (length higher than
      zero) to IDX files. These files must be formatted according to the
      requirements specified in the 'help' section of `Repeats`.

    Output
    ======

    Returns a dictionary of sequence sizes wich can be directly passed to the
    init function of the `Repeats` class.
    """
    # Create the `seqsizes_dict` dictionary required by the class `Repeats`.
    seqsizes_dict = dict()
    for index in list_idx_paths:
        df_indexfile = pd.read_table(index, index_col=0, sep="\s+")
        seqsizes_dict[df_indexfile.index.name] = (
            df_indexfile.iloc[:,0].to_dict())

    return seqsizes_dict

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
                Printing("Cannot read the provided path to `"+
                         str(file)+"`").error()
                return None
        # Ensure there is one idxfile for each RM.out file:
        if len(files_dict.values()) != len(seqsizes_dict.keys()):
            Printing("Index file and RepeatMasker's output files "+
                     "are not found in a one-to-one proportion; "+
                     "make sure that there is one index file per "+
                     "RM.out file").error()
            return None

        # Read input tabulated files. Create dfs with pandas.
        dfs_catalog = dict()
        Printing("Checking file sizes:").status()
        Printing("Keep in mind that the first three rows are considered to be "+
                 "the header and, consequently, they are stripped...").status()
        for species, filepath in files_dict.items():
            # Print file size before reading.
            file_size = os.path.getsize(filepath)
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
            d = dfs_catalog[species] # Short-hand variable.

            # Decrease memory usage by dropping unnecessary columns:
            d.drop(columns=["left", "begin_match", "ID",
                            "end_match", "left_match"], inplace=True)

            # Create a new column where element length will be computed. As well
            # as percent substitutions ("perc_divg"), it could be used as a
            # proxy for rep. ele. divergence. However, substitutions are more
            # straightforward and should be more accurate. Later on it will be
            # used to compute the cumulative sum of RE's base pairs.
            d["replen"] = abs(d["end"] - d["begin"]) +1
            # In order to concatenate all dataframes in `dfs_catalog` add a new
            # column which contains the species label.
            dfs_catalog[species]["Species"] = species

            # Decrease memory usage by using efficient 'dtypes' (see
            # <https://pandas.pydata.org/docs/user_guide/scale.html>)

            # Select float columns and try to downcast them.
            ##float_cols = d.select_dtypes("float").columns
            ##d[float_cols] = d[float_cols].apply(
            ##        pd.to_numeric, downcast="float")
            ## Downcasting integers makes it impossible to round them?

            # Select integer columns and try to downcast them.
            int_cols = d.select_dtypes("integer").columns
            d[int_cols] = d[int_cols].apply(
                    pd.to_numeric, downcast="integer")

            # Object dtypes are more expensive than category dtypes. Try to
            # downcast from objects to categories, if possible.
            obj_cols = d.select_dtypes("object").columns
            d[obj_cols] = d[obj_cols].astype("category")

        # Concatenate RM.out files from all species into a single df.
        df = pd.concat(list(dfs_catalog.values()))
        Printing("Information about the loaded DataFrame from the "+
                 "given RepeatMasker output:").status()
        df.info(memory_usage="deep")

        # Count rows matching each sequid regex in seqsizes_dict. Make sure all
        # rows match with a single regex; not more, not less.
        Printing("Revising sequid regex matches between provided index file "+
                 "and RepeatMasker's output...").status()
        for species in seqsizes_dict.keys():
            df_species = df.loc[df["Species"] == species]
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
##                print(" + Regex matching list:",
##                      list(df_species_and_regex["sequid"].unique()))
                # Add a new col to keep track of 'sequid_types'
                df.loc[
                    (df["Species"] == species) &
                    (df["sequid"].str.contains(sequid_regex, case=False)),
                    "sequid_type"] = sequid_regex
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
        # As soon as all sequid types have been added, downcast the
        # `sequid_type` column to the category dtype.
        obj_cols = df.select_dtypes("object").columns
        df[obj_cols] = df[obj_cols].astype("category")
        # Round float columns, like "perc_divg". Otherwise pd.Series.mode() is
        # very unprecise. RepeatMasker returns up to one decimal of precision,
        # so... it shouldn't affect much?
        df = df.round(1)

        # Reclassify repeat types into 4 new columns.
        self.reclassify_default_reptypes(df)

        # Check whether the reclassification correctly labelled all default
        # repeat types.
        Printing("Checking the reclassification/standardisation of "+
                 "repeat types").status()
        self.df_reclassification = self.check_reclassification(df)

        # Summarise the dataframe (RM.out) into a table with absolute bp
        # counts per species, sequid and repeat type combinations.
        # Furthermore, use multiple algorithms/methods to obtain bp count.
        # As some columns are categorical, it is important to specify the
        # parameter `observed=True` (otherwise returns all catg. combinations).

        # "NAIVE" count of base pairs (includes all REs):
        Printing("Computing count of base pairs using the 'NAIVE' "+
                 "method").status()
        df_sum_naive = df.groupby([
            "Species", "sequid_type",
            "class", "subclass", "order", "superfam"], observed=True)[
                "replen"].agg(naive_numele="count", # amount of REs
                              naive_bpsum="sum")    # base pairs of REs
        # "BEST" count of base pairs (remove overlapping repeats from the sum):
        Printing("Computing count of base pairs using the 'BEST' "+
                 "method").status()
        df_sum_best = df.loc[df["overlapping"] == False].groupby([
            "Species", "sequid_type",
            "class", "subclass", "order", "superfam"], observed=True)[
                "replen"].agg(best_numele="count", # amount of REs
                              best_bpsum="sum")    # base pairs of REs

        # Merge both summary dataframes into a single summary dataframe.
        df_summary = df_sum_naive.join([df_sum_best]).reset_index()
        # "ALGORITHMICAL" count of base pairs. Remove only the overlapping
        # regions instead of the entire RE.
        Printing("Computing count of base pairs using the 'ALGOR' "+
                 "method").status()
        # `algorithmically_bpsum` updates df_summary in-place.
        df_summary = self.algorithmically_bpsum(df_summary, df)

        # Print the sum of overlapping bps computed by "BEST" and "ALGOR"
        # methods in comparison to "NAIVE".
        Printing("Printing the differences between 'NAIVE', 'BEST' and "+
                 "'ALGOR' methods caused by overlapping repeats (ie. "+
                 "a measure of the prevalence of overlapping repeats).").status()
        self.evaluate_overlapping(df_summary)

        Printing("Computing the total repetitive and nonrepetitive "+
                 "fractions per `species` and `sequid_type` pairs").status()
        df_summary = self.compute_nonrepetitive(df_summary, seqsizes_dict)

        Printing("Computing summary dataframes with different aggregate "+
                 "levels of 'repeat types' (from class to "+
                 "superfamily).").status()
        dict_df_summary = self.branching_summaries(df_summary)

        Printing("Computing estimators of 'sample' divergence: median, "+
                 "mean and mode divergences (including standard dev.)").status()
        for key, dfsum in dict_df_summary.items():
            dict_df_summary[key] = self.estimators_divergence_summary(dfsum, df)
            # Sort the resulting dataframes by repeat classes
            # using a custom key.
            custom_key_order = {"DNA": 10,
                                "Retrotransposon": 11,
                                "Other": 12,
                                "Tandem_repeat": 13,
                                "Unclassified": 14,
                                "Repetitive_fraction": 15,
                                "Nonrepetitive_fraction": 16, }
            dict_df_summary[key] = dict_df_summary[key].sort_values(
                by=["class"],
                key=lambda x: x.replace(custom_key_order))
            # Sort the dataframe by 'Species' and 'sequid_type'. Use the kind of
            # sort 'mergesort' to preserve the previous sort (stable algorithm).
            # Start by reading whether Species/sequid columns are in the df.
            list_cols = dict_df_summary[key].columns
            list_sorting_cols = [x for x in list_cols if x in [
                "Species", "sequid_type"]]
            # Once we have made sure these columns are in the df, sort it:
            dict_df_summary[key] = dict_df_summary[key].sort_values(
                by=list_sorting_cols, kind="mergesort")

            # We want to introduce np.nan to pandas.DataFrame "int" column;
            # Nonrepetitive_fraction cannot have a count of REs (numele)!!!
            # np.nan values are floats. The dtype has to be changed to "Int64"
            # (capitalised "Int"!!!).
            int_cols = dict_df_summary[key].select_dtypes("integer").columns
            dict_df_summary[key][int_cols] = (
                dict_df_summary[key][int_cols].astype("Int64"))
            dict_df_summary[key].loc[
                dict_df_summary[key]["class"]=="Nonrepetitive_fraction",
                ["naive_numele", "best_numele"]] = pd.NA

        # Print the top of the most common "perc_divg" values (for debug
        # purposes; compare with pd.Series.mode)
        #print(df["perc_divg"].value_counts(bins=None, normalize=True))
        #print(df["perc_divg"].value_counts(bins=20, normalize=True))

        Printing("Writing summary DataFrames to the dictionary (Python var"+
                 "iable) `self.dict_df_summary`.").status()
        self.dict_df_summary = dict_df_summary
        Printing("Writing entire RepeatMasker DataFrame to "+
                 "the variable `self.df`.").status()
        self.df = df

    def reclassify_default_reptypes(self, df):
        """
        Very long function of many paragraphs. Each paragraph reclassifies a
        default repeat type into four standardised clades, which can be easily
        read and analysed.

        Some repeat types are reclassified if the "default_repclass" field
        exactly matches a single string (Y==X), while other repeat types are
        reclassified if the "default_repclass" matches/contains a regex-like
        string (Y.str.contains(X)).

        These modifications are done in-place to the specified `df`.

        Input
        =====

        + df: A dataframe from the __init__ function of this class.

        Output
        ======

        A modified df, with more columns, standardising the repeat types.
        Modifications are done in-place.
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
        ] = ("Retrotransposon", "NA", "DIRS", "DIRS")

        df.loc[
            df["default_repclass"].str.contains("LTR") &
            df["default_repclass"].str.contains("Ngaro"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "DIRS", "Ngaro")

        df.loc[
            df["default_repclass"].str.contains("LTR/Viper"),
            ["class", "subclass", "order", "superfam"]
        ] = ("Retrotransposon", "NA", "DIRS", "Viper")

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
        ] = ("DNA", "1", "NA", "Zator")

        df.loc[
            df["default_repclass"].str.contains("ClassII_MITE"),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "NA", "MITE")

        df.loc[
            df["default_repclass"].str.contains("ClassII_nMITE"),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "NA", "nMITE")

        df.loc[
            df["default_repclass"].str.contains("DNA") &
            df["default_repclass"].str.contains("Academ"),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "TIR", "Academ")

        df.loc[
            (df["default_repclass"].str.contains("DNA") &
            df["default_repclass"].str.contains("CMC")) |
            df["default_repclass"].str.contains("ClassII_DNA_CACTA"),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "TIR", "CACTA")

        df.loc[
            (df["default_repclass"].str.contains("DNA") &
            df["default_repclass"].str.contains("Dada")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "TIR", "Dada")

        df.loc[
            (df["default_repclass"].str.contains("DNA") &
            df["default_repclass"].str.contains("Kolobok")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "TIR", "Kolobok")

        df.loc[
            (df["default_repclass"].str.contains("DNA") &
            df["default_repclass"].str.contains("Ginger")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "TIR", "Ginger")

        df.loc[
            (df["default_repclass"].str.contains("DNA") &
            df["default_repclass"].str.contains("MULE")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "TIR", "MULE")

        df.loc[
            (df["default_repclass"].str.contains("DNA") &
            df["default_repclass"].str.contains("Merlin")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "TIR", "Merlin")

        df.loc[
            (df["default_repclass"].str.contains("DNA/P")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "TIR", "P")

        df.loc[
            (df["default_repclass"].str.contains("DNA/PIF")) |
            (df["default_repclass"].str.contains("ClassII_DNA_Harbinger")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "TIR", "PIF-Harbinger")

        df.loc[
            (df["default_repclass"].str.contains("PiggyBac")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "TIR", "PiggyBac")

        df.loc[
            (df["default_repclass"].str.contains("DNA/Sola")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "TIR", "Sola")

        df.loc[
            (df["default_repclass"].str.contains("DNA/Sola")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "TIR", "Sola")

        df.loc[
            (df["default_repclass"].str.contains("DNA/TcMar") |
            df["default_repclass"].str.contains("DNA_TcMar")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "TIR", "TcMar")

        df.loc[
            (df["default_repclass"].str.contains("DNA") &
            df["default_repclass"].str.contains("hAT")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "TIR", "hAT")

        df.loc[
            (df["default_repclass"].str.contains("DNA_Mutator")),
            ["class", "subclass", "order", "superfam"]
        ] = ("DNA", "1", "TIR", "Mutator")

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
        obj_cols = df.select_dtypes("object").columns
        df[obj_cols] = df[obj_cols].astype("category")

        return df

    def remove_overlapping(self, intervals):
        """
        Parse the intervals (begin, end) of RE (repetitive elements) of a given
        sequid, computing the cumulative sum whilst removing overlapping base
        pairs. Interval length is computed as `(end-begin)+1`. Intervals with
        the best score are prioritised over repeats with lower score.

        Input
        =====

        + intervals: A list of intervals (which in turn are sublists composed of
          `begin` and `end` coordinates). The pairs of coordinates must be
          integers. Score for each interval must be provided as the third item
          in the sublist. Include a tuple with the classification of the RE.

        intervals = [ [begin_1, end_1, score_1, classification_1], ...,
        [begin_N, end_N, score_N, classification_N] ]

        classification = (class_str, subclass_str, order_str, superfam_str)

        Output
        ======

        The cumulative sum of interval lengths, each computes as
        `(end-begin)+1`. The function removes overlapping segments, counting
        each nucleotide a single time.
        """
        # Init variables.
        # Create a list with all unique repeat type tuples.
        unique_class_tuples = list(set(tuple(x[3]) for x in intervals))
        unique_class_tuples.sort()
        # Convert into pd.MultiIndex()
        idx = pd.MultiIndex.from_tuples(unique_class_tuples, names=[
                'class', 'subclass', 'order', 'superfam'])
        # Create a df which will be returned as answer.
        df_resulting = pd.DataFrame(index=idx)
        df_resulting['algor_bpsum'] = 0
        # Create a point-list from an interval-list
        # from [begin, end], [begin, end], etc.
        # to [coord, begin, 0], [coord, end, 0], [coord, begin, 1], etc.
        point_list = []
        for i in range(0, len(intervals)):
            # Left/starting points labeled '0'
            point_list += [ [intervals[i][0], 0, i] ]
            # Right/ending points labeled '1'
            point_list += [ [intervals[i][1], 1, i] ]
        # Sort point list by position (by first value in sublists)
        point_list.sort()
        # Initialise other algorithm variables...
        currentBest = -1
        open_intervals = []

        # Start iterating across the point_list, discounting overlapping base
        # pairs.
        for i in range(0, len(point_list)):
            # If the loop landed on a left point (0) opens an interval
            if point_list[i][1] == 0:

                # If there is no other interval opened:
                if currentBest == -1:
                    # Enters interval 'i'
                    currentBest = point_list[i][2]
                    currentBegin = int(point_list[i][0])
                    currentScore = int(intervals[currentBest][2])

                # Else, there already was an open interval:
                # (and it is of lower score than currentBest)
                elif currentScore > intervals[point_list[i][2]][2]:
                    # Do not close currentBest, because the next interval (repeat)
                    # is of lower score; however, append to open intervals list
                    iden = point_list[i][2]
                    open_intervals.append(
                        [int(intervals[iden][2]),     # score
                         iden])                       # ID
                    open_intervals.sort(reverse=True) # sort by score

                # Else, currentScore should be higher
                else:
                    # Compute length up until now (begin2 -begin1)
                    interval_length = point_list[i][0] -currentBegin
                    # Add this length to the resulting pandas df
                    reptype = intervals[currentBest][3]
                    df_resulting.loc[reptype, "algor_bpsum"] += interval_length
                    # Append index of past Best to open_intervals
                    iden = currentBest
                    open_intervals.append(
                        [int(intervals[iden][2]),     # score
                         iden])                       # ID
                    open_intervals.sort(reverse=True) # sort by score
                    # And change currentBest
                    currentBest = point_list[i][2]
                    currentBegin = point_list[i][0]
                    currentScore = intervals[currentBest][2]

            # If the loop landed on a right point (1) closes an interval
            # Moreover, it has to be the right point of the currently open
            # interval 'i'.
            elif point_list[i][2] == currentBest:
                # DEBUG
    ##            print('point_list[i][2] (', point_list[i], ') == currentBest')
                # Compute length up until now (end -begin +1)
                interval_length = point_list[i][0] -currentBegin +1
                # Add this length to the resulting pandas df
                reptype = intervals[currentBest][3]
                df_resulting.loc[reptype, "algor_bpsum"] += interval_length
                # Close this interval and open the next best one in open_intervals
                if len(open_intervals) == 0:
                    currentBest = -1
                else:
                    # Remove second best from open_intervals and use it as
                    # currentBest
                    c = open_intervals.pop(0)
                    currentBest = c[1]
                    currentScore = c[0]
                    currentBegin = point_list[i][0] +1

            # Otherwise, it is closing an interval listed in `open_intervals`
            else:
                # Remove the interval from `open_intervals`
                iden = point_list[i][2]
                open_intervals.remove(
                    [int(intervals[iden][2]), # score
                    iden])                    # ID

        return df_resulting.reset_index()

    def algorithmically_bpsum(self, df_summary, df):
        """
        In conjunction with `remove_overlapping` function, creates a summary
        column which will be added to the dataframe `df_summary`. This column
        will be composed of the cumulative sums of each group of repeats, whilst
        accounting for REs overlapping.
        """
        # Initialise the new column
        df_summary["algor_bpsum"] = 0
        # Subset the dataframe for each pair of sequid, species:
        for species in df["Species"].unique():
            df_species = df.loc[df["Species"] == species]
            # Iterate across all sequids found in `df_species`.
            for seq in df_species["sequid"].unique():
                df_sequid = df_species.loc[df_species["sequid"] == seq]
                # Acquire the `sequid_type` of `seq` by reading the first cell
                # of the column "sequid_type". It will be used as a category to
                # summarize base pairs to (much like "NAIVE" and "BEST" methods).
                sequid_type = df_sequid.iloc[0]["sequid_type"]
                # Furthermore, group by repeat types:
                tuple_repclass = tuple(zip(
                    list(df_sequid["class"]),
                    list(df_sequid["subclass"]),
                    list(df_sequid["order"]),
                    list(df_sequid["superfam"]) ))
                # Acquire the list of intervals which will be feed to the
                # function `remove_overlapping`.
                intervals = list(zip(
                    list(df_sequid["begin"]),
                    list(df_sequid["end"]),
                    list(df_sequid["score_SW"]),
                    tuple_repclass ))
                # Compute bp count for every type of repeat in sequid in
                # species.
                df_sum_algor = self.remove_overlapping(intervals)
                # Insert a column with sequid and species analyzed in for-loop.
                df_sum_algor.insert(0, "Species",
                                    [species]*df_sum_algor.shape[0])
                df_sum_algor.insert(1, "sequid_type",
                                    [sequid_type]*df_sum_algor.shape[0])
                df_summary = pd.merge(right=df_summary, left=df_sum_algor,
                                      # preserve right dataframe's keys
                                      how="right", on=["Species", "sequid_type",
                                    "class", "subclass", "order", "superfam"])
                df_summary["algor_bpsum_x"] = (
                    df_summary["algor_bpsum_x"].fillna(0))
                df_summary["algor_bpsum_y"] = (
                    df_summary["algor_bpsum_x"] + df_summary["algor_bpsum_y"])
                df_summary.drop(columns=["algor_bpsum_x"], inplace=True)
                df_summary.rename(columns={"algor_bpsum_y": "algor_bpsum"},
                                  inplace=True)

        return df_summary

    def check_reclassification(self, df):
        """
        Checks how the default repeat types have been reclassified. Returns a
        dataframe pairing every unique default type with their assigned class,
        subclass, order and superfamily.

        Input
        =====

        + df: A dataframe generated in the __init__ function by reading RM
          output.

        Output
        ======

        Prints and returns a pandas dataframe in which default types are paired
        with their reassigned classification.
        """
        # Drop rows with non-unique cell values.
        df_unique_classes = df.drop_duplicates(
            subset=["default_repclass", "class",
                    "subclass", "order", "superfam"])
        # Sort the dataframe by the newly assigned types.
        df_unique_classes = df_unique_classes.sort_values(
            by=["class", "subclass", "order", "superfam"])
        # Reset the dataframe index
        df_unique_classes = df_unique_classes.reset_index()
        # We are only interested in columns with the repeat types;
        # drop the rest of columns
        df_unique_classes = df_unique_classes[["class", "subclass", "order",
                                               "superfam", "mite",
                                               "default_repclass"]]
        # Print the dataframe in markdown format (without writing to a file)
        print(df_unique_classes.to_markdown(None))

        return df_unique_classes

    def evaluate_overlapping(self, df_summary):
        """
        Compute the sum of base pairs obtained by each method ("NAIVE", "BEST",
        "ALGOR"). and their percentages relative to "NAIVE" (how much
        overlapping is accounted by each method). Of these three methods, "BEST"
        and "ALGOR" should never surpass "NAIVE" in base pair count.
        """
        for species in df_summary["Species"].unique():
            df = df_summary.loc[df_summary["Species"] == species]
            # Compute the sum of these three columns/methods.
            naive_bpsum = df["naive_bpsum"].sum()
            best_bpsum = df["best_bpsum"].sum()
            algor_bpsum = df["algor_bpsum"].sum()
            # Print results.
            print(" --- Overlapping summary in "+species+" ---")
            print(" + Naive sum of base pairs: "+str(naive_bpsum)+
                  " / 100.0%")
            print(" + Best sum of base pairs: "+str(best_bpsum)+
                  " / "+str(round((best_bpsum/naive_bpsum)*100, 1))+"%")
            print(" + Algor. sum of base pairs: "+str(algor_bpsum)+
                  " / "+str(round((algor_bpsum/naive_bpsum)*100, 1))+"%")
            print() # Aesthetic empty line between species and at the end.

        return None

    def compute_nonrepetitive(self, df_summary, seqsizes_dict):
        """
        Compute non-repetitive and total repetitive fractions per `sequid_type`
        (single chromosomes or multiple minor sequids specified in the
        idx-file).
        """
        for species in df_summary["Species"].unique():
            df_species = df_summary.loc[df_summary["Species"] == species]
            for seq in df_species["sequid_type"].unique():
                df_sequid = df_species.loc[df_species["sequid_type"] == seq]
                # Compute rep. fraction for the `located` dataframe cells.
                rep_fractions = {
                    "naive_bpsum": int(df_sequid["naive_bpsum"].sum()),
                    "naive_numele": int(df_sequid["naive_numele"].sum()),
                    "best_bpsum": int(df_sequid["best_bpsum"].sum()),
                    "best_numele": int(df_sequid["best_numele"].sum()),
                    "algor_bpsum": int(df_sequid["algor_bpsum"].sum()), }
                # Add non-repetitive fraction to the previously computed
                # `rep_fractions`.
                new_rows = {
                    "Species": [species] *2,
                    "sequid_type": [seq] *2,
                    "class": ["Repetitive_fraction",
                              "Nonrepetitive_fraction"],
                    "subclass": ["NA"] *2, "order": ["NA"] *2,
                    "superfam": ["NA"] *2,
                    # it is not possible to combine "NA" in an int column;
                    # instead of "NA" num. of ele. use zero.
                    "naive_numele": [rep_fractions["naive_numele"], pd.NA],
                    "naive_bpsum": [rep_fractions["naive_bpsum"],
                        # Non-repetitive = (SEQ.LENGTH - REP.LENGTH)
                                    seqsizes_dict[species][seq] -
                                    rep_fractions["naive_bpsum"]],
                    # instead of "NA" num. of ele. use zero.
                    "best_numele": [rep_fractions["best_numele"], pd.NA],
                    "best_bpsum": [rep_fractions["best_bpsum"],
                                   seqsizes_dict[species][seq] -
                                   rep_fractions["best_bpsum"]],
                    "algor_bpsum": [rep_fractions["algor_bpsum"],
                                    seqsizes_dict[species][seq] -
                                    rep_fractions["algor_bpsum"]],
                }
                df_summary = pd.concat([df_summary, pd.DataFrame(new_rows)])

        return df_summary.sort_values(
            by=["Species", "sequid_type"]).reset_index(drop=True)

    def aggregate_dfsum_by_agglist(self, df_summary, agg_list):
        """
        Aggregate the summary dataframe by a given list of columns.

        Input
        =====

        + df_summary: Dataframe created in the __init__ function.

        + agg_list: A list of columns found in the given `df_summary`; for
          instance, `["Species", "sequid_type", "class"]`.
        """
        # Compute summary by classes:
        df_sum_by_agglist = df_summary.groupby(agg_list)[
                ["naive_numele", "naive_bpsum", "best_numele",
                "best_bpsum", "algor_bpsum"]].agg([
                    "sum"]).reset_index()
        # Remove MultiIndex column names.
        df_sum_by_agglist.columns = df_sum_by_agglist.columns.droplevel(1)

        return df_sum_by_agglist

    def branching_summaries(self, df_summary):
        """
        Uses `aggregate_dfsum_by_agglist` to aggregate the summary dataframe
        into multiple datasets. The information remains the same, but is
        approached at the level of "class", "subclass", "order", and
        "superfamily".
        """
        answer = {"df_aggby_spec_seq_class": self.aggregate_dfsum_by_agglist(df_summary,
            ["Species", "sequid_type", "class"]),
                  "df_aggby_spec_seq_class_subc": self.aggregate_dfsum_by_agglist(df_summary,
            ["Species", "sequid_type", "class", "subclass"]),
                  "df_aggby_spec_seq_class_subc_ord": self.aggregate_dfsum_by_agglist(df_summary,
            ["Species", "sequid_type", "class", "subclass", "order"]),
                  "df_aggby_spec_seq_class_subc_ord_supf": self.aggregate_dfsum_by_agglist(df_summary,
            ["Species", "sequid_type", "class", "subclass", "order", "superfam"]),
                  "df_aggby_spec_class": self.aggregate_dfsum_by_agglist(df_summary,
            ["Species", "class"]),
                  "df_aggby_spec_class_subc_ord_supf": self.aggregate_dfsum_by_agglist(df_summary,
            ["Species", "class", "subclass", "order", "superfam"]),
                  }
        return answer

    def estimators_divergence_summary(self, df_summary, df):
        """
        Obtain median, mean and mode divergences per category in the summary
        dataframe. It will compute these estimators from "perc_divg" column.

        Input
        =====

        + df: A dataframe from the __init__ function of this class.

        + df_summary: A summary (aggregated) dataframe created in the __init__
          function of this class.
        """
        # Aconsegueix llistat de columnes que calen agrupar.
        list_columns = list(df_summary.columns)
        acceptable_columns = ["Species", "sequid_type", "class", "subclass",
                              "order", "superfam"]
        # Remove non-categorical/non-object/non-acceptable columns.
        # (I'm using a strange reversal range loop; from end to begin).
        for i in range(len(list_columns)-1, -1, -1):
            if list_columns[i] not in acceptable_columns:
                list_columns.pop(i)
        # Groupby the list of columns.
        df_sum_diverg = df.groupby(list_columns, observed=True)[
            "perc_divg"].agg(divg_median="median", divg_mean="mean",
        # The mode() function could return an array of divergence with multiple
        # 'most common' values; take the first and the last one in two columns.
                            divg_modeBeg=lambda x: pd.Series.mode(x)[0],
                            divg_modeEnd=lambda x: pd.Series.mode(x)[
                                len(pd.Series.mode(x))-1])
        # Join the estimators of divergence with the categories in `df_summary`
        df_summary = df_summary.merge(df_sum_diverg.reset_index(),
                                     how="outer",
                                     left_on=list_columns,
                                     right_on=list_columns)
        # Make sure that the dtypes of `df_summary` are integers/floats where
        # necessary.
        dict_assign_dtypes = {"naive_numele": "int", "naive_bpsum": "int",
                              "best_numele": "int", "best_bpsum": "int",
                              "algor_bpsum": "int", "divg_median": "float64",
                              "divg_mean": "float64", "divg_modeBeg": "float64",
                              "divg_modeEnd": "float64"}
        df_summary = df_summary.astype(dict_assign_dtypes)
        # Approximate amount of substitutions (base pairs * divergence):
        df_summary["divg_median_x_algor_bpsum"] = (
            df_summary["divg_median"] *
            df_summary["algor_bpsum"] )/100

        return df_summary

def plot_histolike_recursively(
    df, value_column: str,
    categorical_columns: list, hueby_column: str=None,
    title: str="", yaxis_relative=True, minmaxing: int=10):
    """
    Plot a column `value_column` of a pandas DataFrame. Plot local minima and
    maxima as red and green (respectively) scatter points.

    Input
    =====

    + title: Serves as the plot title. Furthermore, it is added as a prefix to
      the filename of all generated PNG files.

    + minmaxing: An integer >=1 or zero/False. It is passed to the parameter
      `order` of `scipy.signal.argrelextrema`. It determines the sensitivity of
      the algorithm; the higher the `minmaxing` value, the less sensitive (less
      local minima or maxima). Adjust this value depending on the prevalence of
      noise in your datasets. To disable the search of local minima and maxima,
      pass a value of zero/False.
    """
    # For each `category` in the zeroeth column of the list
    # `categorical_columns`:
    for category in df[categorical_columns[0]].unique():
        # Filter `df` by `category` within each for-loop.
        bool_filter = (df[categorical_columns[0]] == category)
        df_filtered = df.loc[bool_filter]
        # Obtain the title that will be used in the current plot.
        current_title = title + category
        # Try to find further categorical subdivisions in the list of
        # categorical columns. Count all combinations of unique categories.
        checkval = df_filtered.drop_duplicates(
            subset=categorical_columns)[categorical_columns].shape[0]
        # If `checkval` is above 1 there are further subdivisions of the `df`
        # (more than one combination). Run this function recursively at a lower
        # level (one category down the list) with the filtered `df`.
        if checkval > 1:
            plot_histolike_recursively(
                # Run with filtered df.
                df=df_filtered,
                value_column=value_column,
                # Run one category down the list.
                categorical_columns=categorical_columns[1:],
                hueby_column=hueby_column,
                title=current_title + "_",
                minmaxing=minmaxing,
                yaxis_relative=yaxis_relative)

        # Now, the plotting section.
        # If a hueby_column was specified, start by iterating across it.
        if hueby_column:
            for hue in df_filtered[hueby_column].unique():
                # Filter by `hue` categories and create a Series of
                # value_counts.
                valcounts_series = (
                    df_filtered.loc[df_filtered[hueby_column] == hue,
                                    value_column]
                    ).value_counts(normalize=yaxis_relative).sort_index(
                    ).to_frame()
                # Try to find local min/max points.
                if minmaxing:
                    valcounts_series["min"] = valcounts_series.iloc[argrelextrema(
                        valcounts_series[value_column].values, np.less_equal,
                        order=minmaxing,)][value_column]
                    valcounts_series["max"] = valcounts_series.iloc[argrelextrema(
                        valcounts_series[value_column].values, np.greater_equal,
                        order=minmaxing,)][value_column]
                    # Plot max as green and min as red scatter points.
                    plt.scatter(valcounts_series.index,
                                valcounts_series["min"], c="r")
                    plt.scatter(valcounts_series.index,
                                valcounts_series["max"], c="g")
                # Lastly, with the value counts, plot them as a line plot.
                plt.plot(valcounts_series.index, valcounts_series[value_column],
                         label=hue)
        else:
            pass # only works hueing atm
            # Es podria crear una columna buida '1' i destriar segons aquesta
            # columna (irresponsable?).

        plt.title(current_title)
        plt.legend()
        plt.tight_layout()
        plt.savefig(current_title+".png", dpi=500,)
        # Close the plt.axes or they will leak to the next figures.
        plt.close()

    return None

def plot_histogram_recursively(
    df, value_column: str,
    categorical_columns: list, hueby_column: str=None,
    title: str="", yaxis_relative=True, ):
    """
    Plot a column `value_column` of a pandas DataFrame. Plot local minima and
    maxima as red and green (respectively) scatter points.

    Input
    =====

    + title: Serves as the plot title. Furthermore, it is added as a prefix to
      the filename of all generated PNG files.

    + minmaxing: An integer >=1 or zero/False. It is passed to the parameter
      `order` of `scipy.signal.argrelextrema`. It determines the sensitivity of
      the algorithm; the higher the `minmaxing` value, the less sensitive (less
      local minima or maxima). Adjust this value depending on the prevalence of
      noise in your datasets. To disable the search of local minima and maxima,
      pass a value of zero/False.
    """
    # For each `category` in the zeroeth column of the list
    # `categorical_columns`:
    for category in df[categorical_columns[0]].unique():
        # Filter `df` by `category` within each for-loop.
        bool_filter = (df[categorical_columns[0]] == category)
        df_filtered = df.loc[bool_filter]
        # Obtain the title that will be used in the current plot.
        current_title = title + category
        # Try to find further categorical subdivisions in the list of
        # categorical columns. Count all combinations of unique categories.
        checkval = df_filtered.drop_duplicates(
            subset=categorical_columns)[categorical_columns].shape[0]
        # If `checkval` is above 1 there are further subdivisions of the `df`
        # (more than one combination). Run this function recursively at a lower
        # level (one category down the list) with the filtered `df`.
        if checkval > 1:
            plot_histogram_recursively(
                # Run with filtered df.
                df=df_filtered,
                value_column=value_column,
                # Run one category down the list (recursive aspect)
                categorical_columns=categorical_columns[1:],
                hueby_column=hueby_column,
                title=current_title + "_",
                yaxis_relative=yaxis_relative)

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
                valcounts_series['hue'] = hue
                df_plot = pd.concat([df_plot, valcounts_series])
            # Lastly, with the value counts, plot them as an histogram.
            sns.histplot(x=df_plot.index, weights=df_plot[value_column],
                         hue=df_plot["hue"], binwidth=1, multiple="layer",
                         # Colors for each species (it should not be a static
                         # variable, but parameterised by the function).
                         palette={"Dcat": "#eb8f46", "Dtil": "#30acac",
                                  "b": "brown"})
        else:
            pass # only works hueing atm
            # Es podria crear una columna buida '1' i destriar segons aquesta
            # columna (irresponsable?).

        plt.title(current_title)
        #plt.legend()
        plt.tight_layout()
        plt.savefig(current_title+".png", dpi=500,)
        # Close the plt.axes or they will leak to the next figures.
        plt.close()

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
    parser.add_argument("--rmlist", nargs="+", required=True,
                        help="A list of RepeatMasker's output files.")
    parser.add_argument("--idxlist", nargs="+", required=True,
                        help="A list of index files, as specified in --help")
    parser.add_argument("--summary",
                        help="A suffix (string) with which a path will be "+
                        "created. Then, multiple TSV summary tables will "+
                        "be written to different paths starting with the "+
                        "given suffix. Do not include the TSV extension!")
    parser.add_argument("--plots",
                        help="A suffix (string) with which a path will be "+
                        "created. Then, multiple PNG plots will be written "+
                        "to different paths starting with the given suffix. "+
                        "Do not include the PNG extension!")

    args = parser.parse_args()
    # call a value: args.operacio or args.filename.

    # Create the `files_dict` dictionary required by the class `Repeats`.
    files_dict = dict()
    for pair in args.rmlist:
        values = pair.split(",")
        if len(values) != 2:
            Printing("One of the specified pairs in --rmlist was not of "+
                     "length equal to two").error()
            Printing("Use --help for an explanation on --rmlist "+
                     "parameter requirements").error()
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
    # Write table summaries
    if args.summary:
        # Store reclassification dataframe
        repeats.df_reclassification.to_csv(
            str(args.summary)+"_reclassification.tsv", # file-name
            sep="\t", na_rep="NA", index=True)
        # Store summary dataframes.
        for key, dfsum in repeats.dict_df_summary.items():
            add_to_path = key.split("_")[2:]
            add_to_path = '_'.join(add_to_path)
            filename = str(args.summary) + "_" + add_to_path + ".tsv"
            dfsum.to_csv(
                filename, sep="\t", na_rep="NA", index=False, decimal=".",
                # Write up to five decimal places; keep trailing zeroes.
                float_format="%.5f")
    # Write PNG figures
    if args.plots:
        plot_histogram_recursively(
            df=repeats.df, value_column="perc_divg",
            categorical_columns=["class", "order", "superfam"],
            hueby_column="Species", title=args.plots + "_", )



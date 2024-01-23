#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# repellion_v2.py
#
# 17 de gen. 2024  <adria@molevol-OptiPlex-9020>

# Development
# -----------
#
# Top10 value counts for the "perc_divg" column (and other variable columns, like
# replen, if possible). Hand a subset of a dataframe (eg. loc DNA). Otherwise,
# hand a list of groupby columns (["class", "Species"]).
# 
# Divergence histograms. Try to include max-mins (inflection points using SciPy).

help_msg = """
RM.out files (--rmlist)
-----------------------

RepeatMasker's output files. They must be paired with a species nickname/label,
which must match the species label found in the index file. For instance:

script.py --rmfiles file1,sp1 [file2,sp2]

Where file1 and file2 are different files obtained as RepeatMasker's output, and
sp1 and sp2 are species labels. Separate file path and species label with a
comma.

Index files (--idxlist)
-----------------------

TSV (Tab-separated Values) files with sequence regexes paired with the sum of
their lengths. The header (first row) must contain an species nickname. For
example:

```idxfile
sp1
Chromosome1\t32
Scaff\t5600
etc.
```

Where `sp1` is an species label (matching the label of an RM.out file) and
`Chromosome1` and `Scaff` are regex patterns that will match sequence IDs  in
the corresponding column of the RM.out file. Another example:

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
---------------------------

The output obtained from RepeatMasker must be modified before using it.

Some rows contain an asterisk (*) as a 16th field. Some other rows do end at the
15th field. To make sure that all rows are of the same length, one can run:

```bash
awk 'BEGIN{OFS="\t"}
{if ( $16 ) { $16="True" } else { $16="False" } }
{print}' RMfile.out > RMfile_modified.out
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
# Reading input files' size and creating folders with figures.
import os

# Handling input as DataFrames.
import pandas as pd, numpy as np
# Plotting figures.
import seaborn as sns, matplotlib.pyplot as plt
# Obtain local minima/maxima values to plot them.
from scipy.signal import argrelextrema
# Measure the time it takes to run this script.
from datetime import datetime
import locale
locale.setlocale(locale.LC_TIME, '') # sets locale to the one used by user

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def read_idxfile_paths(list_idxfiles):
    """
    Read the input IDX files and prepare them to be piped to the class
    `Repeats`, formatting the data according to `Repeats` requirements.

    Input
    -----

    + list_idxfiles: A list with one or multiple paths (strings) pointing to IDX
      files. These files must be formatted according to the requirements
      specified in the `--help` page.

    Output
    ------

    Returns a dictionary of sequence sizes which can be directly passed onto the
    init function of the `Repeats` class.
    """
    seqsizes_dict = dict()
    for idxfile in list_idxfiles:
        df_idxfile = pd.read_table(idxfile, index_col=0, sep="\s+")
        seqsizes_dict[df_idxfile.index.name] = \
            df_idxfile.iloc[:, 0].to_dict()

    return seqsizes_dict

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

class Repeats(object):
    """
    Read RepeatMasker's output, standardise the RE types classification, and
    return summary tables of repeat content.

    RepeatMasker's documentation:
    <https://www.repeatmasker.org/webrepeatmaskerhelp.html>

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
      Use '|' for OR matching, for instance: `string1|string2|string3`,
      `Scaff|ctg`...
      Use '&' for AND matching, for instance: `Sca&ffolds`...

    Output
    ------

    TBD

    Methods
    -------

    TBD
    """
    def __init__(self, files_dict, seqsizes_dict):
        # Store the dictionary with sequence sizes for later use.
        self.seqsizes_dict = seqsizes_dict
        # Compute whole species' genome sizes (sum of bp of all sequences).
        self.gensizes_dict = dict()
        for species in seqsizes_dict.keys():
            self.gensizes_dict[species] = sum(seqsizes_dict[species].values())

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
        Printing("Checking file sizes:").status()
        Printing("Keep in mind that the first three rows are considered to be "+
                 "the header and, consequently, they are stripped...").status()
        df = self._read_files(files_dict)  # Init the definitive df to storage the data.

        Printing("Printing information about the loaded DataFrame from the "+
                 "given RepeatMasker's output files:").status()
        df.info(memory_usage="deep")

        # Count rows matching each sequid regex in `seqsizes_dict`. Make sure
        # that all rows match with a single regex; not more, not less.
        Printing("Revising sequid regex matches between provided index file "+
                 "and RepeatMasker's output:").status()
        df = self._assign_sequid_regexes(df)

        # Round float columns, like "perc_divg". Otherwise pd.Series.mode() is
        # very unprecise. RepeatMasker returns up to one decimal of precision,
        # so... it shouldn't affect much?
        #>df = df.round(1)

        # Reclassify RE types into 4 new columns (unpack and standardise default
        # types).
        Printing("Reclassifying RE types.").status()
        df = self._reclassification(df)

        # Check whether the reclassification correctly labelled all default
        # repeat types.
        Printing("Checking the reclassification/standardisation of "+
                 "repeat types:").status()
        self.df_check_reclassification = self._check_reclassification(df)

        # There are multiple approaches to computing the total amount of base
        # pairs occupied by a RE type/clade. For example, the sum of all base
        # pairs (column "replen"), or the sum of the base pairs of high-scoring,
        # non-overlapping REs, or the sum of non-overlapping base pairs (column
        # "algorithm").
        # Let us compute base pairs with the "algorithm" method.
        Printing("Computing count of base pairs using the 'algorithm' "+
                 "method.").status()
        df = self._algorithm_basepair_sum(df)

        # Print the sum of basepairs computed by these methods.
        Printing("Printing the differences between the three used "+
                 "methods. These differences are caused by artefactually "+
                 "overlapping repeats.").status()
        self.df_evaluate_overlapping = self.evaluate_overlapping(df)

        # Initialise a catalog of summary dataframe. Each one considers a
        # different grouping of the dataset.
        Printing("Grouping the main dataframe by categorical columns, "+
                 "thereby creating summarised dataframes.").status()
        self.catalog_df_summaries = dict()
        for main_categories in [["Species"],
                                ["Species", "sequid_type"]]:
            for secondary_group in [["class"],
                                    ["class", "subclass"],
                                    ["class", "subclass", "order"],
                                    ["class", "subclass", "order", "superfam"]]:
                categories = list(main_categories) + list(secondary_group)
                self.catalog_df_summaries["_".join(categories)] = \
                    self.generate_summary_table(df, categories)

        # Create a new catalog of summaries with a different (more concise)
        # style.
        Printing("Applying a style to the summary dataframes.").status()
        self.catalog_df_summaries_concise = dict()
        for key in self.catalog_df_summaries.keys():
            self.catalog_df_summaries_concise[key] = \
                self.style_julio_dysdera(self.catalog_df_summaries[key])

        Printing("Writing entire RepeatMasker DataFrame to "+
                 "the variable `self.df`.").status()
        self.df = df

    def _read_files(self, files_dict):
        """
        Hidden method that reads the RepeatMasker output/dataset and generates a
        pd.DataFrame() from them.
        """
        df = pd.DataFrame()
        for species, filepath in files_dict.items():
            # Print file size before reading.
            file_size = os.path.getsize(filepath)
            print(f" * {filepath}:  {file_size} bytes")
            # Load the file as pd.DataFrame and temporally store in `d`.
            d = pd.read_table(filepath,
                    header=None,  # Skip the first three rows.
                    skiprows=3,   # Fields separated by 'space' regex.
                    sep="\s+",    # Use these strings as column names.
                    names=[
                    # For an explanation of RepeatMasker's fields, see its
                    # manual.
                    "score_SW", "perc_divg", "perc_del", "perc_ins", "sequid",
                    "begin", "end", "left", "orient", "name",
                    "default_repclass", "begin_match", "end_match",
                    "left_match", "ID", "overlapping"],
                )

            # Decrease memory usage by using efficient 'dtypes' (see
            # <https://pandas.pydata.org/docs/user_guide/scale.html>)

            # Decrease memory usage by dropping unnecessary columns:
            d.drop(columns=[
                "left", "begin_match", "ID", "end_match", "left_match"],
                   inplace=True)

            # Object dtypes are more expensive than category dtypes. Try to
            # downcast from objects to categories, if possible.
            obj_cols = d.select_dtypes("object").columns
            d[obj_cols] = d[obj_cols].astype("category")

            # Create a new column where element length will be computed. RE length
            # could be used as a proxy for divergence. However, substitutions
            # (column "perc_divg") are more straightforward and should be more
            # accurate and direct.
            d["replen"] = abs(d["end"] - d["begin"]) +1

            # Add a new column which labels the dataset by species.
            d["Species"] = species

            # Concatenate this data file to `df`
            df = pd.concat([df, d])

        return df

    def _assign_sequid_regexes(self, df):
        """
        Count rows matching each sequid regex in `seqsizes_dict`. Make sure
        that all rows match with a single regex; not more, not less.
        """
        df = df.copy()

        for species in df["Species"].unique():
            mask_species = df["Species"] == species
            # Print total number of rows and total number of unique sequids.
            species_totnrows = df.loc[mask_species].shape[0]
            species_totuniq = \
                len(list(df.loc[mask_species, "sequid"].unique()))
            # Print the number of rows matching and the number of unique
            # sequids that did match.
            species_matchrows = 0
            species_matchuniq = 0
            for sequid_regex in self.seqsizes_dict[species].keys():
                # Filter dataframe by sequids matching `sequid_regex`
                mask_seqreg = df.loc[mask_species, "sequid"].str.contains(
                    sequid_regex, case=False)
                print("--", sequid_regex, "--")
                print(" + Nrows:", df.loc[mask_seqreg].shape[0])
                species_matchrows += df.loc[mask_seqreg].shape[0]
                print(" + Regex matching count:",
                      len(df.loc[mask_seqreg, "sequid"].unique()),
                      "unique sequids")
                species_matchuniq += len(
                    df.loc[mask_seqreg, "sequid"].unique())
                #>print(" + Regex matching list:",
                #>      list(df.loc[mask_seqreg, "sequid"].unique()))

                # Add a new column to keep track of sequid types/regexes.
                df.loc[mask_seqreg, "sequid_type"] = sequid_regex

            # Print a summary of the species, including all sequids.
            print("-----------------------------")
            print(species, "summary: matched", species_matchrows,
                  "out of", species_totnrows, "rows (",
                  round((species_matchrows/species_totnrows) *100, ndigits=3),
                  "%)")
            print(" " * len(species + "summary: "), "found",
                  species_matchuniq, "out of", species_totuniq,
                  "unique sequids (",
                  round((species_matchuniq/species_totuniq) *100, ndigits=3),
                  "%)\n")
        # As soon as all sequid types have been added, downcast the
        # `sequid_type` column to the category dtype.
        obj_cols = df.select_dtypes("object").columns
        df[obj_cols] = df[obj_cols].astype("category")

        return df

    def _reclassification(self, df):
        """
        Very long function with many paragraphs. Each paragraph reclassifies a
        default repeat type into four standardised clades, whcih can be easily
        read and analysed.

        Some repeat types are reclassified based on an exact match in the field
        "default_repclass", while others are reclassified based on containing a
        substring/pattern.
        """
        df = df.copy()

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
        obj_cols = df.select_dtypes("object").columns
        df[obj_cols] = df[obj_cols].astype("category")

        return df

    def _algorithm_basepair_sum(self, df):
        """
        Compute the cumulative sum of basepairs while accounting for REs
        overlapping.
        """
        df = df.copy()

        # Initialise a dictionary to temporally store the results of the
        # algorithm.
        answer = dict()
        # For species and sequid:
        for species in df.loc[:, "Species"].unique():
            mask_species = df["Species"] == species
            for sequid in df.loc[mask_species, "sequid"].unique():
                # Acquire the list of intervals (begin, end) for the current
                # `sequid`, which will be feed to the algorithm.
                # Keep in mind that two pairs of coordinates cannot overlap if
                # they are in different sequences!
                mask_sequid = df.loc[mask_species, "sequid"] == sequid
                intervals = list(zip(
                    list(df.loc[mask_sequid, "begin"]),
                    list(df.loc[mask_sequid, "end"]),
                    list(df.loc[mask_sequid, "score_SW"]),
                    list(df.loc[mask_sequid,].index), ))
                # Sort interval list by fourth field (index).
                intervals = sorted(intervals, key=lambda x: x[3])
                # Produce a dict. from `intervals` with the function
                # `_remove_intervals_overlappings`. Merge it to `answer` with
                # the operand "|" or "merge dictionaries", ie `dict | dict`.
                answer = answer | \
                    self._remove_intervals_overlappings(intervals)

        # Lastly, create a pd.Series from the previous dictionary. Merge it by
        # index with the main dataframe, as if it were a new column.
        answer = pd.Series(answer, name="algor_bp", dtype="int")

        return df.merge(answer, left_index=True, right_index=True)

    def _remove_intervals_overlappings(self, intervals):
        """
        Take a list of intervals and remove overlapping base pairs based on
        Smith-Waterman (SW) score.

        Input
        -----

        + intervals: A list of intervals (which, in turn, are sublists
          themselves). Intervals are composed of `begin` and `end` coordinates
          as the first and second items. The pair of coordinates must be
          integers. Smith-Waterman score and `df.index` must be provided as the
          third and fourth items in the sublist. For instance:

          [ [beg_1, end_1, score_1, idx_1], ...,
            [beg_N, end_N, score_N, idx_N] ]

        Output
        ------

        Returns a dictionary where keys are a row's index and values are a row's
        amount of basepairs using the algorithm method.
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

        return answer

    def _check_reclassification(self, df):
        """
        Checks how the default RE types/names have been reclassified.
        Returns a dataframe pairing every unique default type/name with their
        assigned class, subclass, order and superfamily.

        Input
        -----

        + df: the main dataframe created by the class `Repeats`.

        Output
        ------

        Prints and returns a pandas dataframe in which default types are paired
        with their reassigned classification.
        """
        # While using the main dataframe (self.df), create a subset of it by
        # dropping non-unique cell values.
        df_unique_classes = df.drop_duplicates(
            subset=["default_repclass", "class",
                    "subclass", "order", "superfam"])
        # Sort the dataframe by the newly assigned types.
        df_unique_classes = df_unique_classes.sort_values(
            by=["class", "subclass", "order", "superfam"])
        # Reset the dataframe index.
        df_unique_classes = df_unique_classes.reset_index()
        # We are only interested in columns with the repeat types;
        # drop the rest of columns.
        df_unique_classes = df_unique_classes[["class", "subclass", "order",
                                               "superfam", "mite",
                                               "default_repclass"]]
        # Print the dataframe in markdown format (without writing to a file)
        print(df_unique_classes.to_markdown(None))

        return df_unique_classes

    def evaluate_overlapping(self, df):
        """
        Compute the sum of base pairs obtained by the three different methods:

        1) Sum of all base pairs (column "replen").
        2) Sum of the higher SW-score REs, removing the lower SW-score
           REs (column "replen" filtered with column "overlapping").
        3) Sum of the higher SW-score REs, including the non-overlapping
           fragments of lower SW-score REs (column "algor_bp").

        Input
        -----

        + df: the main dataframe created by the class `Repeats`.

        Output
        ------

        Prints and returns a pandas dataframe with the sum of basepairs obtained
        from each method.
        """
        answer = list()
        # Compute sum per species.
        for species in df["Species"].unique():
            mask_species = df["Species"] == species
            method_1_bpsum = df.loc[mask_species, "replen"].sum()
            method_2_bpsum = df.loc[
                (mask_species) &
                (~ df["overlapping"]), "replen"].sum()
            method_3_bpsum = df.loc[mask_species, "algor_bp"].sum()
            genome_size = self.gensizes_dict[species]
            # Create a dataframe with the sum of base pairs.
            df_answer = pd.DataFrame(
                data={
                    "Species": [species]*3,
                    "basepairs": [
                    method_1_bpsum,
                    method_2_bpsum,
                    method_3_bpsum]},
                index=["naive", "best", "algor"])
            # Compute percentages relative to `method_1_bpsum` (naive) and
            # relative to the species' genome.
            df_answer["perc_comparison"] = \
                (df_answer["basepairs"]/method_1_bpsum)*100
            df_answer["perc_genome"] = (df_answer["basepairs"]/genome_size)*100
            # Store the species' summary to return at the end.
            answer.append(df_answer)

        # Concatenate the dataframes of all "Species" into a single one.
        answer = pd.concat(answer)
        # Print basepair sum per species as a df in markdown format (without
        # writing to a file).
        print(" --- Overlapping summary ---")
        print(answer.to_markdown(None))

        return answer

    def generate_summary_table(self, df, categorical_columns):
        """
        + categorical_columns: A list of categorical column names (strings)
          found in `df`. They will be used to group the dataframe. Info:
          <https://pandas.pydata.org/docs/user_guide/groupby.html>

        """
        # Compute the total count of REs and their genome occupancy (in bp) with
        # the "naive" approach (does not discard low scoring REs).
        summary_naive = df.groupby(categorical_columns,
                             as_index=False, observed=True) \
            .agg(
            # Aggregate the categorical groups (sum/count each group's rows)
                naive_numele=pd.NamedAgg(column="replen", aggfunc="count"),
                naive_totbp=pd.NamedAgg(column="replen", aggfunc="sum"),
                naive_algbp=pd.NamedAgg(column="algor_bp", aggfunc="sum"))
        # Compute a partial count of REs and their genome occupancy (in bp) with
        # the "best" approach (discards low scoring REs).
        df_filtered = df.loc[~ df["overlapping"]]
        summary_best = df_filtered.groupby(categorical_columns,
                             as_index=False, observed=True) \
            .agg(
            # Aggregate the categorical groups (sum/count each group's rows)
                best_numele=pd.NamedAgg(column="replen", aggfunc="count"),
                best_totbp=pd.NamedAgg(column="replen", aggfunc="sum"),
                best_algbp=pd.NamedAgg(column="algor_bp", aggfunc="sum"))

        # Compute median, mean and mode estimators of divergences per category.
        # It will compute these estimators from the column "perc_divg".
        summ_div_naive = df.groupby(categorical_columns,
                             as_index=False, observed=True) \
            .agg(
            # Aggregate the categorical groups.
                naive_div_medi=pd.NamedAgg(column="perc_divg",
                    aggfunc="median"),
                naive_div_mean=pd.NamedAgg(column="perc_divg",
                    aggfunc="mean"),
                # The pd.Series.mode() function could return an array of
                # divergence with multiple 'most common' values; take the first
                # (Beg) and the last one (End) in two columns.
                naive_div_modeBeg=pd.NamedAgg(column="perc_divg",
                    aggfunc=lambda x: pd.Series.mode(x)[
                        0]),
                naive_div_modeEnd=pd.NamedAgg(column="perc_divg",
                    aggfunc=lambda x: pd.Series.mode(x)[
                        len(pd.Series.mode(x))-1]),
                # Apply the floor of the series of values (np.floor) and,
                # afterwards, get its mode (mode from the intervals/bins of one
                # percent of divergence).
                naive_div_modeFloor=pd.NamedAgg(column="perc_divg",
                    aggfunc=lambda x: pd.Series.mode(np.floor(x))[0]))

        # Merge the multiple summaries into a single table/df.
        df_summary = summary_naive.merge(
            summary_best).merge(
                summ_div_naive)

        # Aggregate the total RE fraction (genome occupancy) of each species.
        # Consider all RE types as a single RE category. Discover whether
        # 'species' and 'sequid_type' were included as categories.
        main_agg_groups = \
            [x for x in categorical_columns if \
             x in ["Species", "sequid_type"]]
        secondary_group = \
            [x for x in categorical_columns if \
             x not in ["Species", "sequid_type"]]
        if len(main_agg_groups)>0:
            # Same block as above, but for `main_agg_groups` only instead of the
            # whole of `categorical_columns`.
            totfrac_summary_naive = df.groupby(main_agg_groups,
                                 as_index=False, observed=True) \
                .agg(
                    naive_numele=pd.NamedAgg(column="replen", aggfunc="count"),
                    naive_totbp=pd.NamedAgg(column="replen", aggfunc="sum"),
                    naive_algbp=pd.NamedAgg(column="algor_bp", aggfunc="sum"))
            totfrac_summary_best = df_filtered.groupby(main_agg_groups,
                                 as_index=False, observed=True) \
                .agg(
                    best_numele=pd.NamedAgg(column="replen", aggfunc="count"),
                    best_totbp=pd.NamedAgg(column="replen", aggfunc="sum"),
                    best_algbp=pd.NamedAgg(column="algor_bp", aggfunc="sum"))

            totfrac_summ_div_naive = df.groupby(main_agg_groups,
                                 as_index=False, observed=True) \
                .agg(
                    naive_div_medi=pd.NamedAgg(column="perc_divg",
                        aggfunc="median"),
                    naive_div_mean=pd.NamedAgg(column="perc_divg",
                        aggfunc="mean"),
                    naive_div_modeBeg=pd.NamedAgg(column="perc_divg",
                        aggfunc=lambda x: pd.Series.mode(x)[
                            0]),
                    naive_div_modeEnd=pd.NamedAgg(column="perc_divg",
                        aggfunc=lambda x: pd.Series.mode(x)[
                            len(pd.Series.mode(x))-1]),
                    naive_div_modeFloor=pd.NamedAgg(column="perc_divg",
                        aggfunc=lambda x: pd.Series.mode(np.floor(x))[0]))

            totfrac_summary = totfrac_summary_naive.merge(
                totfrac_summary_best).merge(
                    totfrac_summ_div_naive)
            for column in secondary_group:
                if column == "class":
                    totfrac_summary["class"] = "Repetitive_fraction"
                else:
                    totfrac_summary[column] = pd.NA

            # Compute non-repetitive basepairs (complimentary of repetitive
            # fraction).
            nonrep_summary = totfrac_summary.copy()
            if "class" in nonrep_summary.columns:
                nonrep_summary["class"] = "Nonrepetitive_fraction"
            nonrep_summary[["naive_numele", "best_numele",
                            "naive_div_medi", "naive_div_mean",
                            "naive_div_modeBeg", "naive_div_modeEnd",
                            "naive_div_modeFloor"]] = pd.NA
            if main_agg_groups == ["Species"]:
                # For each column with a base pair count, and for each row in
                # the summary dataframe, compute the difference between species'
                # genome and repetitive fraction... with list comprehension!
                for basepair_column in ["naive_totbp", "naive_algbp",
                                        "best_totbp", "best_algbp"]:
                    nonrep_summary[basepair_column] = [
                        self.gensizes_dict[sp] - repfrac \
                        for sp, repfrac in \
                        zip(totfrac_summary["Species"],
                            totfrac_summary[basepair_column])]
            elif main_agg_groups == ["Species", "sequid_type"]:
                # Same as the previous `if` statement, but includes sequid where
                # necessary.
                for basepair_column in ["naive_totbp", "naive_algbp",
                                        "best_totbp", "best_algbp"]:
                    nonrep_summary[basepair_column] = [
                        self.seqsizes_dict[sp][sequid] - repfrac \
                        for sp, repfrac, sequid in \
                        zip(totfrac_summary["Species"],
                            totfrac_summary[basepair_column],
                            totfrac_summary["sequid_type"])]

            # Concatenate the total repetitive fraction to `df_summary`.
            df_summary = pd.concat([df_summary, totfrac_summary, nonrep_summary])

        # Revise mode columns: is there multiple "most frequent" value?
        # If so, assign them `pd.NA`.
        df_summary["naive_div_modeRaw"] = 0  # Initialise column.
        df_summary.loc[
            df_summary["naive_div_modeEnd"]==df_summary["naive_div_modeBeg"],
                       "naive_div_modeRaw"] = (
            df_summary.loc[
            df_summary["naive_div_modeEnd"]==df_summary["naive_div_modeBeg"],
                       "naive_div_modeBeg"] )
        df_summary.loc[
            df_summary["naive_div_modeEnd"]!=df_summary["naive_div_modeBeg"],
                       "naive_div_modeRaw"] = pd.NA
        # Drop the now redundant "Beg" and "End" mode columns.
        df_summary = df_summary.drop(columns=["naive_div_modeBeg",
                                              "naive_div_modeEnd"])

        return df_summary.reset_index(drop=True) \
            .sort_values(by=categorical_columns)

    def style_julio_dysdera(self, df_summary):
        """
        Applies the table aesthetics chosen by a dearest PI of mine, JR, to a
        dearest spider genus of mine, *Dysdera*.

        Input
        -----

        + df_summary: Summary dataframe obtained by grouping by categorical
          columns using the function `generate_summary_table`.

        Output
        ------

        Returns a stylized table.
        """
        # Create a copy of the dataframe. Avoid overwriting variables/data.
        df_summary = df_summary.copy()
        # Use a custom sorting key for classes and species.
        custom_key_order = {"DNA": 10, "Retrotransposon": 11, "Other": 12,
                            "Tandem_repeat": 13, "Unclassified": 14,
                            "Repetitive_fraction": 15,
                            "Nonrepetitive_fraction": 16,
                            "Dcatv33": 10, "Dcatv35":11, "Dcat": 12,
                            "Dtil": 13, "Dsil": 14}
        # First of all, sort by the column "class" with the key
        # "custom_key_order".
        df_summary = df_summary.sort_values(
            by="class",
            key=lambda x: x.replace(custom_key_order))
        if "sequid_type" in df_summary.columns:
            # Use the kind of sort 'mergesort' to preserve the previous sort
            # (stable algorithm). Otherwise subsequent sorts do not preserve,
            # for instance, the "class" sort.
            df_summary = df_summary.sort_values(
                by="sequid_type", kind="mergesort")
        df_summary = df_summary.sort_values(
            by="Species", kind="mergesort",
            key=lambda x: x.replace(custom_key_order))

        # Drop excess information. Try to be concise. Find the trascendental
        # datum!
        df_summary = df_summary.drop(columns=["best_algbp"])

        # The column "mode_Floor" points to the floor of the interval. Change it
        # to the midpoint of the interval by summing 0.5, which makes it more
        # intuitive to understand.
        df_summary["naive_div_modeFloor"] = \
            0.5 + df_summary["naive_div_modeFloor"]

        # Change column names.
        df_summary = df_summary.rename(columns={
            "Species": "Species", "sequid_type": "Chromosome",
            "class": "RE_class", "subclass": "RE_subclass",
            "order": "RE_order", "superfam": "RE_superfam",
            "naive_numele": "#TotRE", "naive_totbp": "#TotBP",
            "naive_algbp": "#AlgBP",
            "best_numele": "#ResRE", "best_totbp": "#ResBP",
            "naive_div_medi": "div_median", "naive_div_mean": "div_mean",
            "naive_div_modeRaw": "div_mode_raw",
            "naive_div_modeFloor": "div_mode_bins", })

        return df_summary.reset_index(drop=True)


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
    parser.add_argument("--matplotlib_style",
                        help="The basename of a file in ``")

    args = parser.parse_args()
    # call a value: args.operacio or args.filename.

    # Create the `seqsizes_dict` dictionary required by the class `Repeats`.
    seqsizes_dict = read_idxfile_paths(args.idxlist)
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

    repeats = Repeats(files_dict, seqsizes_dict=seqsizes_dict)


#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# recursive_reptype.py
#
# 07 de des. 2023  <adria@molevol-OptiPlex-9020>

help_msg = """
Describe here the goal, input, output, etc. of the script
as a multi-line block of text.
"""

import sys
from scipy.signal import argrelextrema
import numpy as np
import matplotlib.pyplot as plt

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def recursive_histolike_plots(
    df, value_column: str,
    categorical_columns: list, hueby_column: str=None,
    title: str="", yaxis_relative=True, minmaxing: int=10):
    """
    Input
    =====

    + minmaxing: An integer >=1. It is passed to the parameter `order` of
      `scipy.signal.argrelextrema`. It determines the sensitivity of the
      algorithm; the higher the `minmaxing` value, the less sensitive (less local
      minima or maxima). Adjust this value depending on the prevalence of noise in
      your datasets.
    """
    # For each categorical column in the dataframe (specified by
    # `categorical_columns`:
    for category in df[categorical_columns[0]].unique():
        bool_filter = (df[categorical_columns[0]] == category)
        # Filter the given dataframe by the obtained categories.
        df_filtered = df.loc[bool_filter]
        # Obtain the title that will be used in the plot.
        current_title = title + category
        # Try to find further categorical subdivisions in the list of
        # categorical columns. Count all combinations of unique categories.
        checkval = df_filtered.drop_duplicates(
            subset=categorical_columns)[categorical_columns].shape[0]
        # If `checkval` is above 1 there are further subdivisions of the df
        # (more than one combinations). Run this function recursively at a lower
        # level (one category down the list).
        if checkval > 1:
            recursive_histolike_plots(
                df_filtered, value_column,
                categorical_columns[1:], hueby_column,
                current_title + "_")

        # Now, the plotting section.
        # If a hueby_column was specified, start by iterating across it.
        if hueby_column:
            for hue in df_filtered[hueby_column].unique():
                # Filter by `hue` categories and create Series of value_counts.
                valcounts_series = (
                    df_filtered.loc[df_filtered[hueby_column] == hue,
                                    value_column]
                    ).value_counts(normalize=yaxis_relative).sort_index(
                    ).to_frame()
                # Try to find local min/max points.
                valcounts_series["min"] = valcounts_series.iloc[argrelextrema(
                    valcounts_series[value_column].values, np.less_equal,
                    order=minmaxing,)][value_column]
                valcounts_series["max"] = valcounts_series.iloc[argrelextrema(
                    valcounts_series[value_column].values, np.greater_equal,
                    order=minmaxing,)][value_column]
                # Now, with the value counts, plot them as a line plot.
                plt.plot(valcounts_series.index, valcounts_series[value_column])
                # Plot max as green and min as red scatter points.
                plt.scatter(valcounts_series.index,
                            valcounts_series["min"], c="r")
                plt.scatter(valcounts_series.index,
                            valcounts_series["max"], c="g")
        else:
            pass # only works hueing atm
            # Es podria crear una columna buida '1' i destriar segons aquesta
            # columna (irresponsable).

        plt.title(current_title)
        plt.tight_layout()
        plt.savefig(current_title+".png", dpi=500,)
        # Close the plt.axes or they will leak to the next figures.
        plt.close()


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



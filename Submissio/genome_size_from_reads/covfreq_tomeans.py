#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# coverage.py
#
# 05 de nov. 2022  <adria@molevol-OptiPlex-9020>

"""Reads one or more *.csv-like files, with a single column of frequencies of
coverage ordered by coverage.

Example:
...
covg-file,freq-file
1,1050
2,504
3,306
4,607
5,1272
6,1443
7,1149
etc
...

The first number (1050) is the frequency of reads which had a mapping coverage
of 1, etc.
"""

import sys
# DataFrame reading, parsing, manip...
import pandas as pd
# Plotting
import matplotlib.pyplot as plt
# Square root for stdev
import math

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def freq_to_rows(dataframe, freq_col_name, str_col_name):
    """
    This function is useful to create a dataframe compatible with a Poisson
    distribution. Nevertheless, coverage data is massive and the newly created
    DataFrames consume insanely high amounts of RAM.

    Use only for small frequecy DataFrames.

    INPUT
    =====

    + dataframe: a pandas.DataFrame object.
    + freq_col_name: the name of a column which measures the frequency of
    'str_col_name' column occurrences.
    + str_col_name: the name of a column which displays any str/int. The
    function will create a new pandas.DataFrame with as many as 'freq_col_name'
    rows of each row of 'str_col_name'.
    """
    return pd.DataFrame(
        [row[str_row_name]  # select row to be returned
        for row in df.to_dict(orient='records')
        for _ in range(row[freq_row_name])  # as many times as this range
    ])

def frequency_pandasdf_mean(df, values_colname, freqs_colnames):
    """
    The first line sums all rows of (Cov*Freq).
    The second line sums all rows of (Freq).
    The division is the mean of the distribution.
    """
    return ((df[values_colname] * df[freqs_colnames]).sum() /
        df[freqs_colnames].sum())

def frequency_pandasdf_median(df, values_colname, freqs_colnames):
    """
    """
    return_point = df[freqs_colnames].sum() / 2
    df['cum_sum'] = df[freqs_colnames].cumsum()
    df_median = df[df['cum_sum'] >= return_point]
    return df_median[values_colname].iloc[0]

def frequency_pandasdf_max(df, values_colname, freqs_colnames):
    """
    """
    # Create a subset of the dataframe with all rows where the frequency reaches
    # the maximum(s) values. In case there are multiple frequency maximums.
    max_subset = df[df[freqs_colnames] == df[freqs_colnames].max()]
    # Return the list of coverages with maximum frequency.
    return max_subset[values_colname].tolist()

def frequency_pandasdf_min(df, values_colname, freqs_colnames):
    """
    """
    # Create a subset of the dataframe with all rows where the frequency reaches
    # the minimum(s) values. In case there are multiple frequency minimums.
    min_subset = df[df[freqs_colnames] == df[freqs_colnames].min()]
    # Return the list of coverages with minimum frequency.
    return min_subset[values_colname].tolist()

def frequency_pandasdf_stdev(df, values_colname, freqs_colnames):
    """
    The first line sums all values in the 'freq' row.
    The second line invokes the mean of the DataFrame.
    The third line computes the squared differences between values and mean.
    The fourth line returns the square root of the third line divided by the
    first line (square root of variance).
    """
    denominador = df[freqs_colnames].sum()
    mean = frequency_pandasdf_mean(df, values_colname, freqs_colnames)
    numerador = (((df[values_colname]-mean)**2) * df[freqs_colnames]).sum()
    return math.sqrt(numerador/denominador)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Si es crida com a script:
if __name__ == '__main__':
    # Init the main DF
    df = pd.read_csv(sys.argv[1])[:-1]
    # Make sure all columns are integers (coverage could be binned into
    # >1000 for high values, which is a nuisance for our case).
    # [:-1] drops the last row, usually the one that is binned.
    for (columnName, columnData) in df.items():
        df = df.astype({columnName: 'int'})

    # Make sure that each csv contains a column tagged with 'freq'.
    #~print(df.dtypes)

    # Translate columns.
    for col in df.columns:
        if 'freq' in col.lower():
            df = df.rename(columns={col: 'freq'})
        else:
            df = df.rename(columns={col: 'covg'})
    # Make sure the translation was correct:
    print(df.dtypes)

    # Do analyses.
    print(f"+ Covg. values range: {df['covg'].iloc[0]} -- {df['covg'].iloc[-1]}")
    print(f"+ Mean: {frequency_pandasdf_mean(df, 'covg', 'freq')}")
    print(f"+ Median: {frequency_pandasdf_median(df, 'covg', 'freq')}")
    print(f"+ Max: {frequency_pandasdf_max(df, 'covg', 'freq')}")
    print(f"+ Min: {frequency_pandasdf_min(df, 'covg', 'freq')}")
    print(f"+ Stdev: {frequency_pandasdf_stdev(df, 'covg', 'freq')}")


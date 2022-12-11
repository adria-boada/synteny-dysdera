#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# cas especial: shebang modificat per a hercules...
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
    return max_subset[values_colname].tolist()[0]

def frequency_pandasdf_min(df, values_colname, freqs_colnames):
    """
    """
    # Create a subset of the dataframe with all rows where the frequency reaches
    # the minimum(s) values. In case there are multiple frequency minimums.
    min_subset = df[df[freqs_colnames] == df[freqs_colnames].min()]
    # Return the list of coverages with minimum frequency.
    return min_subset[values_colname].tolist()[0]

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
    numerador = (((df[values_colname] - mean)**2) * df[freqs_colnames]).sum()
    return math.sqrt(numerador/denominador)

def frequency_pandasdf_localmodes(df, values_colname, freqs_colnames):
    """
    Get the local maxima/minima from an array.
    Max.: a value preceded and succeded by other values lower than itself.
    Min.: a value preceded and succeded by other values higher than itself.

    Example:
    [1, 3, 5, 4, 9] -> Third and last values are local maxima.
    [1, 3, 5, 4, 9] -> First and second to last (fourth)
    values are local minima.
    """
    idx_local_maxima=[]
    idx_local_minima=[]
    for i, row in df.iterrows():
        if i==0:
            if df.loc[i,freqs_colnames] < df.loc[i+1,freqs_colnames]:
                idx_local_minima.append(i)
            elif df.loc[i,freqs_colnames] > df.loc[i+1,freqs_colnames]:
                idx_local_maxima.append(i)
        elif i == len(df)-1:
            if df.loc[i,freqs_colnames] < df.loc[i-1,freqs_colnames]:
                idx_local_minima.append(i)
            elif df.loc[i,freqs_colnames] > df.loc[i-1,freqs_colnames]:
                idx_local_maxima.append(i)
        else:
            if df.loc[i,freqs_colnames] < df.loc[i-1,freqs_colnames] and df.loc[i,freqs_colnames] < df.loc[i+1,freqs_colnames]:
                idx_local_minima.append(i)
            elif df.loc[i,freqs_colnames] > df.loc[i-1,freqs_colnames] and df.loc[i,freqs_colnames] > df.loc[i+1,freqs_colnames]:
                idx_local_maxima.append(i)
    return {'idx_max': idx_local_maxima,
            'idx_min': idx_local_minima,
            'val_max': [df.loc[val, values_colname] for val in idx_local_maxima],
            'val_min': [df.loc[val, values_colname] for val in idx_local_minima],
    }

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

    # Make sure that each csv contains two columns: one tagged with 'freq'
    # and the other will be interpreted as coverage at the given 'freq'.
    #~print(df.dtypes)

    # Translate columns. If it cannot translate to integer, throw error.
    for col in df.columns:
        if 'freq' in col.lower():
            df = df.rename(columns={col: 'freq'})
        else:
            df = df.rename(columns={col: 'covg'})
    # Make sure the translation was correct, and columns are
    # entirely of integer type.

    # Do analyses.
    print(f"+ Covg. values range: {df['covg'].iloc[0]} -- {df['covg'].iloc[-1]}")
    print(f"+ Mean: {round(frequency_pandasdf_mean(df, 'covg', 'freq'), 3)}")
    print(f"+ Median: {frequency_pandasdf_median(df, 'covg', 'freq')}")
    print(f"+ First abs. max.: {frequency_pandasdf_max(df, 'covg', 'freq')}")
    print(f"+ First ten Local max.: {frequency_pandasdf_localmodes(df, 'covg', 'freq')['val_max'][0:9]}")
    print(f"+ First abs. min.: {frequency_pandasdf_min(df, 'covg', 'freq')}")
    print(f"+ First ten Local min.: {frequency_pandasdf_localmodes(df, 'covg', 'freq')['val_min'][0:9]}")
    print(f"+ Stdev: {round(frequency_pandasdf_stdev(df, 'covg', 'freq'), 3)}")


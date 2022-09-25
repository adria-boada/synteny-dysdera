#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# statmaths.py
#
# 29 may 2022  <adria@molevol-OptiPlex-9020>

""" Gather simple mathematical and statistical methods,
in order to easily import when needed.
"""

import sys
import math as math

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

K = n = Ex = Ex2 = 0.0


def reinitialize():
    """All variables back to zero. Change sample/population study.
    """
    global K, n, Ex, Ex2
    K = n = Ex = Ex2 = 0.0

def add_measure(x):
    """Add a single measure from the sample to include in
    the analysis. Use inside a 'for i in data' loop in
    order to add all the data to analyze.
    """
    global K, n, Ex, Ex2
    if n == 0:
        K = x
    n += 1
    Ex += x - K
    Ex2 += (x - K) * (x - K)

def remove_variable(x):
    global K, n, Ex, Ex2
    n -= 1
    Ex -= x - K
    Ex2 -= (x - K) * (x - K)

def get_mean():
    """Returns the mean of the sample of added measures
    with `add_variable('measure')` function.
    """
    global K, n, Ex
    return K + Ex / n

def get_variance():
    """Returns the variance of the sample of added measures
    with `add_variable('measure')` function.
    """
    global n, Ex, Ex2
    return (Ex2 - (Ex * Ex) / n) / (n - 1)

def get_stdeviation():
    """Returns the standard deviation of the sample of added
    measures with `add_variable('measure')` function.
    """
    global n, Ex, Ex2
    return math.sqrt( (Ex2 - (Ex * Ex) / n) / (n - 1) )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Si es crida com a script:
if __name__ == '__main__':
    # Instruccions respecte els arguments necessaris per cridar l'script:
    # (substituir arg1, arg2... pel nom dels arguments).
    try:
        measures = [int(x) for x in sys.argv[1:]]
    except (ValueError, IndexError):
        sys.exit("\nCrit script: %s <measure-1> <measure-2> (...) <measure-N> \n" % sys.argv[0])

    # Afageix totes les mesures dels arguments a l'anàlisi;
    for i in measures:
        add_measure(i)

    print('Mean of measures:', get_mean())
    print('Standard Deviation of measures:', get_stdeviation())
    


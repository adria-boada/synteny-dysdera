#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# distance_between_mappings.py
#
# 21 jun 2022  <adria@molevol-OptiPlex-9020>

"""Compute the mean distance between .paf mappings.
A measure of how much interspersed are they.
"""

import sys, statmaths as sms, pafparser as pf, tabulate

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def listing_coordinates(file):
    """List the coordinates of all main chr inside a paf mapping file.
    """
    out={}
    out2={}
    lengths={}
    ali_len={}

    with open(file) as fn:
        for line in fn:
            obj = pf.Mapping(line)

            if ('chr' in obj.qname.casefold()) or ('nc' in obj.qname.casefold()):
                try:
                    out[obj.qname] += [[obj.qstart, obj.qend]]
                except:
                    out[obj.qname] = [[obj.qstart, obj.qend]]
                    lengths[obj.qname] = obj.qlen

            if ('chr' in obj.tname.casefold()) or ('nc' in obj.tname.casefold()):
                try:
                    out[obj.tname] += [[obj.tstart, obj.tend]]
                except:
                    out[obj.tname] = [[obj.tstart, obj.tend]]
                    lengths[obj.tname] = obj.tlen

    for crm, lst in out.items():
        # sort by the first value, ie beginning of each mapping:
        lst.sort(key=lambda x: x[0])
        # re-write sorted list over unsorted one.
        out[crm] = lst

    for crm, lst in out.items():
        # get the first mapping before loop.
        lastS, lastE = lst.pop(0)
        d = []
        ali = []
        for crd in lst:
            if crd[0]-lastE>0:
                # distància entre mapatges.
                d += [crd[0]-lastE]
                # llargada mapatge.
                ali += [lastE-lastS]
                lastS = crd[0]
            lastE = crd[1]
        out2[crm]=d
        ali_len[crm]=ali

    return out2, lengths, ali_len


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Si es crida com a script:
if __name__ == '__main__':

    ### PARSE ARGUMENTS ###
    import argparse, pathlib

    parser = argparse.ArgumentParser(
        description='Calculates mean distance between mappings from a .paf file.')
    # file-name: positional arg.
    parser.add_argument('paf', type=pathlib.Path, help='Path to a PAF mapping file.')
    # integer argument
#    parser.add_argument('-a', '--numero_a', type=int, help='Paràmetre "a"')
    # choices argument
#    parser.add_argument('-o', '--operacio', 
#            type=str, choices=['suma', 'resta', 'multiplicacio'],
#            default='suma', required=False,
#            help='')

    args = parser.parse_args()
    # call a value: args.operacio or args.filename.

    differences_list, lengths, ali_len_list = listing_coordinates(args.paf)

    # Ordena segons les claus del diccionari.
    differences_list = dict(sorted(differences_list.items()))
    ali_len_list = dict(sorted(ali_len_list.items()))

    # header of the table...
    capD = ["Chromosomes", "Min. dist.", "Max. dist.", "Mean dist.",
              "...by Gsize", "St.dev. dist.", "...by Gsize"]
    tableD = []

    capL = ["Chromosomes", "Min. ali. len.", "Max. ali. len.", "Mean ali. len.",
              "...by rel. chr size", "St. dev. ali. len.", "...by rel. chr size"]
    tableL = []

    for crm, lst in differences_list.items():
        for val in lst:
            sms.add_measure(val)

        tableD += [[crm,
                   min(lst),
                   max(lst),
                   sms.get_mean(),
                   str(round((sms.get_mean()*100)/lengths[crm], 5))+" %",
                   sms.get_stdeviation(),
                   str(round((sms.get_stdeviation()*100)/lengths[crm], 5))+" %",
                   ]]
        # reinicia les mesures.
        sms.reinitialize()

    for crm, lst in ali_len_list.items():
        for val in lst:
            sms.add_measure(val)

        tableL += [[crm,
                   min(lst),
                   max(lst),
                   sms.get_mean(),
                   str(round((sms.get_mean()*100)/lengths[crm], 5))+" %",
                   sms.get_stdeviation(),
                   str(round((sms.get_stdeviation()*100)/lengths[crm], 5))+" %",
                   ]]
        # reinicia mesures
        sms.reinitialize()

    print("## Mapping interdistances per chromosome\n")
    print(tabulate.tabulate(tableD, tablefmt='pipe', headers=capD))
    print( # table footer in markdown format...
        ": "+  # these footers open with colons as first character
        "**Min. dist.**: length of the shortest alignment (subsetting by chr). "+
        "**Max. dist.**: length of the longest alignment (subsetting by "+
        "chr). "+
        "**Mean dist.**: mean of alignment lengths (subsetting by chr). "+
        "**St. dev. dist.**: standard deviation of alignment lengths "+
        "(subsetting by chr). "+
        "**By rel. chr size**: variable in previous column (mean or std. dev.) "+
        "divided by corresponding chr size.")
    print()
    print("## Alignment lengths per chromosome\n")
    print(tabulate.tabulate(tableL, tablefmt='pipe', headers=capL))
    print( # table footer in markdown format...
        ": "+  # these footers open with colons as first character
        "**Min. ali. len.**: length of the shortest alignment (subsetting by chr). "+
        "**Max. ali. len.**: length of the longest alignment (subsetting by "+
        "chr). "+
        "**Mean ali. len.**: mean of alignment lengths (subsetting by chr). "+
        "**St. dev. ali. len.**: standard deviation of alignment lengths "+
        "(subsetting by chr). "+
        "**By rel. chr size**: variable in previous column (mean or std. dev.) "+
        "divided by corresponding chr size.")
    print()



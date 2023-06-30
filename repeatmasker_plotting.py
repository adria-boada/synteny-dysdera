#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# repeatmasker_final_countdown.py
#
# 21 de juny 2023  <adria@molevol-OptiPlex-9020>

"""
"""

import sys
import pandas as pd

# some csv are too big to read; they must be split beforehand
# (failed attempt)
import os, csv # reading and splitting csv files

# plotting
import seaborn as sns, matplotlib.pyplot as plt

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def histoing(path_to_some_df):
    """
    """
    df = pd.read_table(path_to_some_df, sep='\t')
    colors = ["grey", "C0", "C3", "C2"]
    df['Species & Genome subset']=df['Species']+" _ "+df['seqtype']

    # remove edgecolor from bars
    plt.rcParams['patch.edgecolor'] = 'none'
    # set the size of the plot (wider than tall)
    plt.figure(figsize=(8, 14))
    # plot
##    plt.margins(y=0.6) # push the bars to the center;
##    ax = sns.histplot(data=df, y="Species", hue="Repeat type", weights="sum",
##                 multiple="stack", shrink=0.3)
##    plt.margins(x=0.1) # push the bars to the center;
    ##sns.set_context("talk")
    ax = sns.histplot(data=df, x='Species & Genome subset', hue="class", weights="bp",
                multiple="stack", shrink=1, palette=colors)
    plt.grid(axis='y')
    # move legend outside of the plotting box
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1,1))
    ax.tick_params(axis='y', labelsize=20, rotation=90)
    ax.tick_params(axis='x', labelsize=15)
    plt.title("Distribution of repetitive and nonrepetitive fraction", y=1.02,
              fontsize=20)
    # adjust margins of figure, so legend, axis, etc.
    # has enough space to be drawn
    plt.subplots_adjust(left=0.3, bottom=0.18, right=0.59, top=0.95)
    plt.ylabel("Gb", fontsize=18)
    plt.xticks(rotation=90)
    plt.xlabel("")
    plt.savefig('histo_alt_stacked_totalbp_v2_matplotlib.png', dpi=300)
    plt.close('all')

    return None

def pieing(path_to_some_df):
    """
    """
    df = pd.read_table(path_to_some_df, sep='\t')
    colors = ["grey", "C0", "C3", "C2"]

    df['perc']=0
    for sp in df["Species"].unique():
        for seq in df["seqtype"].unique():
            mask = (df["Species"]==sp) & (df["seqtype"]==seq)
            tot_bp = int(df.loc[mask,'bp'].sum())
            df.loc[mask, 'perc'] = (df.loc[mask, 'bp']/tot_bp)*100
    df = df.round(1)

    for sp in df["Species"].unique():
        for seq in df["seqtype"].unique():
            mask = (df["Species"]==sp) & (df["seqtype"]==seq)
            fig, ax = plt.subplots()
            ax.pie(df.loc[mask,'perc'], labels=df.loc[mask, 'class'],
                   colors=colors, autopct='%1.1f%%')
            ax.set_title(sp+" - "+seq)
            plt.savefig(sp+"_"+seq+'_v2_pies.png', dpi=300)

    return None

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Si es crida com a script:
if __name__ == '__main__':

    ### PARSE ARGUMENTS ###
    import argparse

    parser = argparse.ArgumentParser(description='')
    # file-name: positional arg.
    parser.add_argument('filename', type=str, help='Path to ... file-name')
    # integer argument
#    parser.add_argument('-a', '--numero_a', type=int, help='Paràmetre "a"')
    # choices argument
#    parser.add_argument('-o', '--operacio', 
#            type=str, choices=['suma', 'resta', 'multiplicacio'],
#            default='suma', required=False,
#            help='')

    args = parser.parse_args()
    # call a value: args.operacio or args.filename.

    histoing(args.filename)
    pieing(args.filename)

    # arithmetics
    df = pd.read_table(args.filename, sep='\t')
    df['perc']=0

    for sp in df["Species"].unique():
        for seq in df["seqtype"].unique():
            mask = (df["Species"]==sp) & (df["seqtype"]==seq)
            tot_bp = int(df.loc[mask,'bp'].sum())
            df.loc[mask, 'perc'] = (df.loc[mask, 'bp']/tot_bp)*100
    df = df.round(1)
    print(df.to_markdown())


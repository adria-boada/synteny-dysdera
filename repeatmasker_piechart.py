#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# repeatmasker_piechart.py
#
# 08 Jun 2022  <adria@molevol-OptiPlex-9020>

"""Import the table resulting from repeatmasker_parser3.py,
analyze it and create a pie-chart.
"""

import sys
import repeatmasker_tandem as rmt
import repeatmasker_parser3 as rm3

# Save figs inside wsl.
import matplotlib
matplotlib.use('Agg') # no UI backend in WSL. 
import matplotlib.pyplot as plt


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def resume_rmtable(rmtable):
    """Outputs a summary of the table from 
    repeatmasker_parser3.py.

    [0]: list of labels
    [1]: list of amount for each label
    [2]: list of length for each label

    List of items in the lists should be ordered in
    corresponding positions.

    There is no 'list of masked length'; you only
    need to add the difference between sum(rep)
    and gsize as a value in the list, labeled as
    functional or non-masked region.
    """
    out = []

    with open(rmtable) as tbl:
        # skip until you find 'Class' in line:
        while 0!=1:
            # readline of file
            line = tbl.readline()
            # if there is 'Class' inside, break the loop.
            if 'Class' in line:
                break
        past = 0
        for line in tbl:
            cls = line.split('|')[1].strip(' ')
            # skip hyphens separators (after header, before total):
            if '-' in cls:
                continue

            # First round of computing:
            if past==0:
                past = cls
                amount = int(line.split('|')[5].strip(' '))
                bp_length = int(line.split('|')[6].strip(' '))
            # From second-row onwards, check that you're summing
            # values from the same class.
            elif cls==past:
                amount += int(line.split('|')[5].strip(' '))
                bp_length += int(line.split('|')[6].strip(' '))
            # If class changes, reinitialize vars and change class-name.
            else:
                out += [ [past, amount, bp_length] ]
                past = cls
                amount, bp_length = int(line.split('|')[5].strip(' ')), int(line.split('|')[6].strip(' '))
                if cls=='Total':
                    out += [ [cls, amount, bp_length] ]

    return out


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Si es crida com a script:
if __name__ == '__main__':
    ### VALORS D'ARGUMENTS A LA TERMINAL ###
    import argparse

    parser = argparse.ArgumentParser(description='Input a table of repeats, index-file and chromosome name. Outputs a piechart.')
    # file-name: positional arg.
    parser.add_argument('repeat_tbl', type=str, 
            help='Table-of-repeats file-name.')
    # file-name: str flagged arg.
    parser.add_argument('-i', '--index_file', type=str, required=True,
            help='FASTA.idx file-name.')
    # file-name: str flagged arg.
    parser.add_argument('-c', '--chr_name', type=str, required=True, 
            help='A chromosome name for Python to match to.')

    args = parser.parse_args()
    # call a value: args.operacio or args.filename.

    # DEBUG func():
    #~print('FUNC()-DEBUG', resume_rmtable(args.repeat_tbl) )

    # Create the piechart for length and amount...
    values = resume_rmtable(args.repeat_tbl)

    # Delete small values and pool them in a new class called 'All the rest':
    all_the_rest = 0
    new_values = []
    for key, am, length in values:
        if (am/values[-1][1])*100 <= 5:
            all_the_rest += am
        elif key=='Total':
            pass
        else:
            new_values += [ [key, am, length] ]
    new_values += [ ['All the rest', all_the_rest, 0.0] ]

    # first plot: get amount (x[1]) and plot with labels (x[0]).
    plt.pie([x[1] for x in new_values], 
            # get all first values
            labels=[x[0] for x in new_values],
            # formatting of numbers in the piechart.
            autopct=lambda p: '{:.2f}% ({:,.0f})'.format(p, p* sum([x[1] for x in new_values])/100)
            )
    plt.title(f'Amount of repeats for Crm. {args.chr_name}')
    plt.savefig(f"amount_pchart.png", bbox_inches='tight')

    # second plot: get length (x[2]) and plot with labels (x[0]).
    plt.pie([x[2] for x in new_values],
            # zip with labels:
            labels=[x[0] for x in new_values],
            # formatting of numbers in the piechart.
            autopct=lambda p: '{:.2f}% ({:,.0f})'.format(p, p* sum([x[1] for x in new_values])/100)
            )
    plt.title(f'Length of repeats for Crm. {args.chr_name}')
    plt.savefig(f"len_pchart.png", bbox_inches='tight')
    

    # Aconsegueix la mida cromosòmica per calcular %masked vals.
    gsize = rmt.refid_to_length(refid_name=args.chr_name, indexed_fasta=args.index_file)

    for val in resume_rmtable(args.repeat_tbl): pass
        # the function resume_rmtable outputs a list of
        # sublists: each sublist contains:
            # [0]: list of labels
            # [1]: list of amount for each label
            # [2]: list of length for each label

    


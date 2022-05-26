
#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# paf_parser.py
#
# 2022-02-21  <adria@molevol-OptiPlex-9020>

""" Count how many matches have length of one,
two, three, etc., until max match length.

Do the same procedure for ins and dels.

Create a function which could be imported.

Add a main call which creates a matplotlib
histogram with the acquired data (from
the latter function).

"""

import sys

# Plots:
import matplotlib.pyplot as plt

# Assigning colours on data-points by value...
import numpy as np


# Instruccions respecte els arguments necessaris per cridar l'script:
if len(sys.argv) < 1:
    sys.exit('\nCrit script: script.py <paf-file.paf>\n')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def paf_cighist (cigar: "A cigar string"):

    """ Input a CIGAR string and output a
    histogram like dictionary, with keys being
    the length of an item and its value the
    amount of items of such length.
    """

    # Create a dict in which to store definitive results.
    results = {
            'M': {},
            'I': {},
            'D': {}
            }

    # Split the cig.
    for i in "MID":
        cigar = cigar.replace(i, i+' ' )
    
    # Remove last space and split by the added separations.
    cigar = cigar.split(' ')[:-1]

    # Compute counts time...
    for i in cigar:
        try:
            # i[-1] should be type.
            # int(i[:-1]) should be that type length.
            results[i[-1]][int(i[:-1])] += 1
        except:
            results[i[-1]][int(i[:-1])] = 1
            
    for key, val in results.items():
        
        results[key] = sorted( val.items() )

    return results


def paf_parse_hist (paf_file: "A *.paf minimap2 alignment"):

    """ Input a paf_file and outputs the whole
    histogram list of categorical and quantity.

    """

    with open(paf_file) as paf:
            
        # El material s'acumula a results.
        results = {
                "M": {},
                "I": {},
                "D": {}
                }

        # Per cada línia del fitxer.
        for line in paf:

            # Cut CIGAR string out.
            c = line.split("\t")[-1][5:-1]

            # Sum histogram-like data.
            for hit_type, histo in paf_cighist( c ).items():
                for x in histo:
                    try:
                        results[hit_type][x[0]] += int(x[1])
                    except:
                        results[hit_type][x[0]] = int(x[1])
        
        for hit_type, histo in results.items():
            # Sort the resulting histogram-like data.
            results[hit_type] = (sorted(histo.items()))

    return results


def plot_cighist (title: "plot super title"= "CIGAR histogram",
        **cig_histlike_dict: "Data in an 'arg= histogram-like' format" ):

    """ Input CIGAR-histogram'ed like data and output
    a histogram plot from matplotlib.
    The input is a list of tuples of length two,
    formatted like (length, amount).
    It is possible to supply more than one histogram in
    a list, and they will be plotted in the same axes
    with different tags.
    """

    # Prepare figure calling matplotlib's .plt
    # Create as many windows as len(input); as many inputs are given.
    fig, axs = plt.subplots(1, len(cig_histlike_dict.keys() ), figsize=(9, 9/len(cig_histlike_dict.keys() ) ), sharey=True)

    counter=0

    # Add each arg to the histogram
    # as a set of data-points with different tags.
    for cigar_type, d in cig_histlike_dict.items():
        
        names = [] # list(data.keys())
        values = [] # list(data.values())

        # separa categoria i quantitat.
        # mantent-les al mateix index. 
        for i in d:
            names += [ i[0] ]
            values += [ i[1] ]

        # Other types of plots, not as pleasing to the eye...
        # axs[0].bar(names, values)
        # axs[1].plot(names, values)

        # Send data to plot number 'counter'.
        axs[counter].scatter(names, values)
        # Send sub-title to plot number 'counter'.
        axs[counter].set_title(cigar_type)
        axs[counter].set_yscale('log')
        axs[counter].set_xscale('log')
        axs[counter].set_ylabel('Frequency')
        axs[counter].set_xlabel('Length of match/gap, in bp')
        counter+=1
        

    # Add title and show the plot.
    fig.suptitle(title)
    plt.show()


def chr_dict_paf (paf_file: "path/to/file to a minimap2 paf alignment output"):
    """
    INPUT
    ------

    A path to a PAF file.

    OUTPUT
    ------
    
    A list with all the scaffolds/contigs/chromosomes 
    included in the alignment, paired in
    sublists with their chromosomal length. 

    [ [DsilChr1, length1], [DsilChr2, length2], *[...] ]   

    """

    out = []

    with open(paf_file) as paf:
        for line in paf:

            # Query and Target name+len for this line
            # (line = one alignment/hit)
            qname = "Q." + line.split("\t")[0]
            qlen = int( line.split("\t")[1] )
    
            tname = "T." + line.split("\t")[5]
            tlen = int( line.split("\t")[6] )

            # If they are not in the out-list,
            # add these chromosomes. 
            if not [qname, qlen] in out:
                out += [ [qname, qlen] ]
            if not [tname, tlen] in out:
                out += [ [tname, tlen] ]

    return out


def locate_transposon (paf_file: "pa/th/file to a minimap2 alignment output" ):
    """
    INPUT
    ------

    A path to a PAF file.


    OUTPUT
    ------

    A list which enumerates hits here and
    hits there, tagged with their position and type.

    [ [position, type, length, chromosome], *[...] ]

    """
    
    # Output list.
    out = []

    with open(paf_file) as paf:
        
        # Per cada línia del fitxer.
        for line in paf:

            # Cut CIGAR string out.
            c = line.split("\t")[-1][5:-1]
            # Cut starting coordinate for the line=hit
            start_coord = {'Q': int(line.split("\t")[2]), 'T': int(line.split("\t")[7]) }
            # Cut out name of aligned CRM.
            qname = "Q.{}".format( line.split("\t")[0] )
            tname = "T.{}".format( line.split("\t")[5] )
            
            # Split the cig. by adding spaces.
            for i in "MID":
                c = c.replace(i, i+' ' )
    
            # Remove last space and split by the added separations.
            c = c.split(' ')[:-1]

            # Guarda cada CIGAR type independentment.
            for i in c:
                # format 'i' és '123M'
                # i[-1] -> 'M'
                # i[:-1] -> '123'
                out += [ [start_coord['Q'], i[-1], int( i[:-1] ), qname], # Query -> [0]
                        [start_coord['T'], i[-1], int( i[:-1] ), tname] ] # Target -> [1]
    
    return out



# If it's called as the main script;
# run and plot a cig histogram.
if __name__ == '__main__': 

    # Switch 'where is my indel' on/off.
    if True:
        # Create test dataset.
        prova = locate_transposon( sys.argv[1] )
        # Format:
        # [ [position, type, length, chromosome], *[...] ]

        # Chromosomes in paf-file:
        # print( chr_dict_paf( sys.argv[1] ))

        # Plot the dataset.
        # for each [crm, llarg.crm]:
        for i in chr_dict_paf( sys.argv[1] ):
            # Separa segons tipus dins
            # del cromosoma selecte...

            # Prepara tres diccionaris amb dades pels eixos X i Y:
            M={'x': [], 'y': []} ; D={'x': [], 'y': []} ; I={'x': [], 'y': []}

            # For each data-point:
            # [coord, type, ll, crm]
            for j in prova:
                # Si crm_i == crm_j:
                if i[0] in j:

                    # separa per tipus:
                    if j[1]=='M':
                        M['x'] += j[0]
                        M['y'] += j[2]

                    elif j[1]=='D':
                        D['x'] += j[0]
                        D['y'] += j[2]

                    else:
                        I['x'] += j[0]
                        I['y'] += j[2]

            # Fer les tres gràfiques,
            # una per a cada tipus:
            # print("PER A ", i, "\n\n", M, "\n", I, "\n", D, sep="")
            for obj in [M, D, I]:
                # Assigna colors depenent del valor Y del punt:
                colour = np.where( obj['y'] < 110, # Condició:
                        'k', # Si y per sota de 110, negre
                        np.where( obj['y'] > 240,
                            'k', # Si y per sobre de 240, negre
                            'r') # else (110 < r < 240) vermell
                        )
                plt.scatter( obj['x'], obj['y'], c=colour )
                plt.show()


    # Switch histogram on/off.
    if False:

        # Create test dataset
        prova = paf_parse_hist( sys.argv[1] )
        # Format of the test.
        #print(prova['M'])
        
        # Create plots from the dataset. 
        plot_cighist( "CIGAR histogram", Match= prova['M'], Deletion= prova['D'], Insertion= prova['I'] )
        plot_cighist( "CIGAR histogram", Deletion= prova['D'], Insertion= prova['I'] )
  




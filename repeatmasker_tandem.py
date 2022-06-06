#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# repeatmasker_tandem.py
#
# 19 may 2022  <adria@molevol-OptiPlex-9020>

""" Divideix un fitxer de repeticions de RepeatMasker en finestres, i assigna a
cada finestra el nombre de repeticions de cada família. 

matplotlib plot: how do we do it?
    * option1: uses R.
    <https://www.biostars.org/p/336999/>
    * option2: needs ideogram.
    <https://www.biostars.org/p/147364/#147637>
    * option3: simple histogram.
    <https://stackoverflow.com/questions/33203645/how-to-plot-a-histogram-using-matplotlib-in-python-with-a-list-of-data>
    * hatched histo: 
    <https://matplotlib.org/stable/gallery/lines_bars_and_markers/filled_step.html?highlight=stacked%20histogram>
    * ideas for sliding window representation:
    <https://bernatgel.github.io/karyoploter_tutorial/Examples/GeneDensityIdeograms/images/Figure9-1.png>

"""

import sys
import os.path
import numpy as np
import repeatmasker_parser3 as rm3
from collections import defaultdict
# colours for matplotlib.
from matplotlib.pyplot import cm

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def unique_refid_list(rmfile):
    """Funció per extreure els refIDs únics d'un fitxer de
    repeticions de RepeatMasker (rmfile)
    """

    out=[]
    
    with gzip.open(rmfile) if rmfile.endswith('.gz') else open(rmfile) as fh:
        for line in fh:
            # Create a Repeat obj.
            obj = rm3.Repeat(line)
            # Si el refId extret és nou:
            if not obj.refid in out:
                out += [obj.refid]

    out.sort()
    return out


def unique_main_refid(rmfile):
    """Funció per extreure els refIDS únics que pertanyin
    a un cromosoma principal; filtrar a fora els cromosomes secundaris.    
    """
    unique_refids = unique_refid_list(rmfile)
    for crm in unique_refid_list(rmfile):
        # Si el cromosoma dins la llista unique_refids
        # no conté 'chr' (és 2ndari):
        if not ('chr'.casefold() in crm.casefold() ):
            # elimina aquests (2ndaris) de la llista:
            unique_refids.remove(crm)

    return unique_refids


def refid_to_length(refid_name, indexed_fasta):
    """
    Inputs: chromosome/scaffold name and an indexed fasta;
    a file with tabulated pairs of (chr, length).

    Returns the length of the inputted chromosome.
    """
    # Obre el fitxer indexed_fasta:
    with open(indexed_fasta) as idx:
        for line in idx:
            # Si la línia conté el 'chr' d'interès:
            if refid_name in line:
                return line.split()[1]


#def matplotlib_abline(slope, intercept):
#    """ Draws an abline with some determined slope and intercept
#    into a matplotlib plot. 
#    """
#    axes = plt.gca()
#    x_vals = np.array(axes.get_xlim())
#    y_vals = intercept + slope * x_vals
#    plt.plot(x_vals, y_vals, '--')


def matplotlib_cutoffs(x, y):
    """Returns which point should be used as a cutoff
    in a matplotlib plot.

    Inputs:
     * X and Y lists of points.

    Outputs:
     * A list of 3 horizontal lines:
     [sup-cutoff, mean, inf-cutoff, out].

     * Out: a list composed of pairs of
     (x_values, y_values) which can be plotted.
    """

    out=[]

    # Calcula mitjana per l'axis Y.
    mean = sum(y) / len(y)
    # Calcula percentils sup. i inf.
    # divideix entre dos el percent (meitat superior, meitat inferior).
    inf_cutoff = (len(y)/100) * (percent/2)
    sup_cutoff = (len(y)/100) * (100-(percent/2))
    # Ordena els punts, per recollir els que estiguin per sobre
    # i per sota del llindar calculat segons percent.
    sortedY = [] + y
    sortedY.sort()
    sortedX = [] + x
    sortedX.sort()
    # Busca el punt amb el qual es farà una línia horitzontal
    # que intercepti amb l'eix Y i separi segons percentils:
    sup_line = sortedY[int(sup_cutoff)]
    inf_line = sortedY[int(inf_cutoff)]

    # Crea els punts per a les tres línies horitzontals:
    # La mitjana al centre, el cutoff superior i l'inferior.
    for intercept in (sup_line, mean, inf_line):
        x_vals = np.array((0, sortedX[-1]))
        y_vals = intercept + 0 * x_vals
        out += [(x_vals, y_vals)]

    return (sup_line, mean, inf_line, out)


#def scatterplot_repeats(x, y, ylab, lbl, axes):
#    """ Create a scatterplot figure with matplotlib.
#
#    Inputs: 
#     * X and Y lists of points,
#     * Name of the ylab (length or amount)
#     * Name of the main title for the plot. 
#     * The created axes with fig, ax = plt.subplots()
#    """
#
#    # Calcula mitjana
#    mean = sum(y) / len(y)
#    # Calcula percentils sup. i inf.
#    # divideix entre dos el percent (meitat superior, meitat inferior).
#    sup_cutoff = (len(y)/100) * (percent/2)
#    inf_cutoff = (len(y)/100) * (100-(percent/2))
#    # Ordena els punts, per recollir els que estiguin per sobre
#    # i per sota del llindar calculat segons percent.
#    sortedY = [] + y
#    sortedY.sort()
#    sortedX = [] + x
#    sortedX.sort()
#    # Busca el punt amb el qual es farà una línia horitzontal
#    # que intercepti amb l'eix Y i separi segons percentils:
#    sup_line = sortedY[int(sup_cutoff)]
#    inf_line = sortedY[int(inf_cutoff)]
#    
#    # Dibuixa les tres línies horitzontals:
#    # La mitjana al centre, el cutoff superior i l'inferior.
#    for intercept in (sup_line, mean, inf_line):
#        x_vals = np.array((0, sortedX[-1]))
#        y_vals = intercept + 0 * x_vals
#        axes.plot(x_vals, y_vals, '--')
#
#    # Si hi han pocs punts (<100), fes línia.
#    if len(x) < 100:
#        data, = axes.plot(x, y, '.-', label=lbl)  
#        # '.-' is the format string: points joint by lines.
#    # Si hi han molts punts, scatterplot sense línies:
#    else:
#        data, = axes.plot(x, y, '.', label=lbl)  # '.' is the format string.
#
#    axes.legend(handles=[data])
#    axes.set_xlabel('Posició dins el cromosoma (b)')
#    axes.set_ylabel(ylab)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Si es crida com a script:
if __name__ == '__main__':
    ### VALORS DELS ARGUMENTS VIA TERMINAL ###
    # Instruccions respecte els arguments necessaris per cridar l'script:
    import argparse

    parser = argparse.ArgumentParser(description='Analyze RepeatMasker output by windows')
    # rmfile name.
    parser.add_argument('rmfile', type=str,
            help='RepeatMasker file-name')
    # chr-len file name.
    parser.add_argument('fasta_idx', type=str,
            help="Tabulated index file with 2 columns: main refids/chromosomes and their length")
    # Window-size.
    parser.add_argument('-w', '--window', type=int, required=True,
            help='Window analysis size')
    # Jump between windows.
    parser.add_argument('-j', '--jump', type=int, 
            help='Distance between neighbours windows',
            default=0, required=False)
    # percent cutoffs.
    parser.add_argument('-p', '--percent', 
            help="A horizontal line will be cutting off the values above and below '-p' %% percentils",
            type=int, required=False, default=5)
    # save or show the figure with matplotlib?
    parser.add_argument('-s', '--microsoft', type=str,
            help="Save the matplotlib figure as a file or show it directly",
            required=False, default='savefig',
            choices=['savefig', 'showfig', 'showfull'])

    args = parser.parse_args()

    # Correcció dels arguments marxa enrere.
    rmfile = args.rmfile
    window_size = args.window
    jump = args.jump
    microsoft = args.microsoft
    percent = args.percent

    # Save figs inside wsl.
    import matplotlib
    matplotlib.use('Agg') # no UI backend in WSL. 
    import matplotlib.pyplot as plt

    ### CREATE A 'RESULTS' DICT AND POPULATE WITH DATA ###
    # create a default dict which stores values
    # 1er, dos cerques,   per repeatclass,    per chromosome,      per posició dins chr: length & amount.
    results = {'general':defaultdict(lambda : defaultdict(lambda : defaultdict(lambda :{'l': 0, 'a': 0}) )),
            'especific': defaultdict(lambda : defaultdict(lambda : defaultdict(lambda :{'l': 0, 'a': 0}) ))
}

    # permet fer el següent sense necessitar generar prèviament claus:
    # results[repeat-type][repeat-class][chr-refid][chr-position][data-type] += 1
    # (això és degut a que inicialitza amb valor = 0)

    # si la llista es troba buida...
    # (no s'ha trobat cap cromosoma dins la llista unique_refids):
    if unique_main_refid(rmfile) == []:
        sys.exit("No s'ha trobat cap cromosoma principal (amb etiqueta 'chr' a la columna refid.) dins el fitxer RepeatMasker.")
    else:
        # Procedeix a fer anàlisi: ensenya els main chr que ha agafat, per si de cas.
        print("Els següents refids contenen 'chr' i se'ls hi farà l'anàlisi:")
        for x in unique_main_refid(rmfile): print(x)

    # Per a cada cromosoma principal que hagis trobat:
    # si els vols tots, substituir unique_main_refid() per la funció unique_refid_list().
    for crm in unique_main_refid(rmfile):
        # parse the file... 
        with gzip.open(rmfile) if rmfile.endswith('.gz') else open(rmfile) as fn:
            for line in fn:
                repobj = rm3.Repeat(line)
                
                # if you find a line with a refid match: 
                if crm in repobj.refid:
                    # repeat position in query => repobj.begin
                    # repeat class => repobj.repclass
                    ### Find the number of window we are found in...
                    modulo = repobj.begin % (jump+window_size)
                    # make sure that the modulo is inside the range [0, Window_Size]:
                    if modulo <= window_size:
                        n = int((repobj.begin-modulo)/(window_size+jump))
                        finestra = n*(jump+window_size)
                        # Pel complet de superfamilies:
                        results['especific'][repobj.repclass[-1]][crm][finestra]['a'] += 1
                        results['especific'][repobj.repclass[-1]][crm][finestra]['l'] += repobj.replen
                        # Per un overview general de classes:
                        results['general'][repobj.repclass[0]][crm][finestra]['a'] += 1
                        results['general'][repobj.repclass[0]][crm][finestra]['l'] += repobj.replen 


    ### AFAGEIX AL DICCIONARI DE RESULTATS LA LLARGADA MÀXIMA DE CADA CRM. ###
    for rep_type, rep_type_val in results.items():
        for rep_class, rep_class_val in rep_type_val.items():
            for crm, crm_val in rep_class_val.items():
                # find max-chr-len with the below function.
                crm_len = refid_to_length(crm, args.fasta_idx)
                # create a new position as the max_crm_len.
                results[rep_type][rep_class][crm][crm_len]['a'] = 0
                results[rep_type][rep_class][crm][crm_len]['l'] = 0


    ### SCATTER-PLOT OF REPEAT FAM LENGTH FOR THE WHOLE GENOME ###
    # agafa els resultats de classes generals...
    for cls, cls_val in results['general'].items():
        # initialize vars to capture results we're interested in...
        last_x_pos = 0
        # Create a new figure for each repeat class.
        fig, axes = plt.subplots()
        # List with cutoffs for each chr.
        cutoffs = []
        # Create an iterable object with all the required colours.
        # You need a colour for each main chromosome/scaffold (~ refid).
        colour = iter(cm.rainbow(np.linspace(0, 1, len(unique_main_refid(rmfile)))))

        # Include all chromosomes in a single figure:
        # Plot the length by position of each chromosome.
        for crm, crm_val in cls_val.items():
            x_pos = []
            y_len = []
            y_amt = []
            for pos, val in crm_val.items():
                # Recopilar valors:
                x_pos += [int(pos)]
                y_len += [int(val['l'])]
                y_amt += [int(val['a'])]
            
            cutoffs += [matplotlib_cutoffs(x_pos, y_len)]
            j = []
            for i in x_pos:
                j += [i + last_x_pos]
            last_x_pos = int(j[-1])

            # Plot either length or amount...
            #~scatterplot_repeats(j, y_len, 'Length of repeats found in window', crm, axL)

            # Get a colour to plot this specific chromosome:
            c = next(colour)
            
            # El bloc de sota feia línies per a gràfiques amb pocs punts...
            ### ### ### ###
#            # Si hi han pocs punts (<100), fes línia.
#            if len(j) < 100:
#                axes.plot(j, y_len, '.-', label=crm, c=c)  
#                # '.-' is the format string: points joint by lines.
#            # Si hi han molts punts, scatterplot sense línies:
#            else:
            ### ### ### ###
            # Millor que sempre siguin scatter-plots.
            axes.plot(j, y_len, '.', label=crm, c=c)  # '.' is the format string.

        # Plot the cutoffs (superior perc., mean, inferior perc.)
        # color for each cutoff:
        colour = iter(cm.rainbow(np.linspace(0, 1, len(unique_main_refid(rmfile)))))
        for cut in cutoffs:
            # L'ordre de les llistes de color i cutoffs haurien de coincidir;
            # el primer valor de cada llista és pel cromosoma A,
            # el segon valor de cada llista és pel cromosoma B...
            c = next(colour)

            # Dibuixa les tres línies horitzontals:
            # La mitjana al centre, el cutoff superior i l'inferior.
            # intercept: a on intercepta la línia horitzontal i eix Y.
            for intercept in (cut[0], cut[1], cut[2]):
                # valors x: la gràfica sencera (de 0 a last val.)
                x_vals = np.array((0, last_x_pos))
                # valors y: el mateix per a tots els punts, l'intercept.
                y_vals = intercept + 0 * x_vals
                axes.plot(x_vals, y_vals, '--', c=c)

        # Crea una llegenda a partir de les etiquetes subministrades 
        # (variable label dins la funció axes.plot()).
        axes.legend()
        axes.set_xlabel('Posició dins el cromosoma (b)')
        axes.set_ylabel('Length of repeats found in window (b)')
        
        title = f"Repeat: {cls.replace('/', '')}\nW. Size: {window_size}; Space between W.: {jump}; %% Cutoffs: {percent}"
        axes.set_title(title)

        # Save a fig for each Repeat class in the results['gen'] dict. 
        plt.savefig(f"length_{cls.replace('/', '')}.png", bbox_inches='tight')

#    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    last_x_pos = 0
#
#    figL, axL = plt.subplots()
#    figA, axA = plt.subplots()
#
#    for crm, crm_val in results.items():
#        for cls, cls_val in crm_val.items():
#            x_pos = []
#            y_len = []
#            y_amt = []
#            for pos, val in cls_val.items():
#                # Recopilar valors:
#                x_pos += [int(pos)]
#                y_len += [int(val['a'])]
#                y_amt += [int(val['l'])]
#
#            j = []
#            for i in x_pos:
#                j += [i + last_x_pos]
#            last_x_pos = int(j[-1])
#
#            title = f"Chromosome: {crm}, Repeat: {cls.replace('/', '')}\nWindow Size: {window_size}, Space between Wndws: {jump}"
#            scatterplot_repeats(j, y_len, 'Length of repeats found in window', title, axL)
#            scatterplot_repeats(j, y_amt, 'Amount of repeats found in window', title, axA)
#    
#    plt.show()
#
#    
#    ### SCATTER-PLOT OF MAIN CHRS AND MAIN CLASSES:
#    for crm, crm_val in results.items():
#        for cls, cls_val in crm_val.items():
#            # initialize vars:
#            X_posicions=[]
#            Y_amount=[]
#            Y_length=[]
#            for pos, val in cls_val.items():
#                # Recopilar valors:
#                X_posicions += [pos]
#                Y_amount += [int(val['a'])]
#                Y_length += [int(val['l'])]
#
#            for Y in ((Y_amount, 'Amount of repeats'), (Y_length, 'Length of repeats found in window')):
#                # Y[0] -> List of data.
#                # Y[1] -> Name of data.
#                fig, ax = plt.subplots()
#                scatterplot_repeats(X_posicions, Y[0], Y[1], ax)
#                # la funció replace() elimina les barres ("/") que confonen amb directoris 
#                if microsoft=='savefig':
#                    # to name the saved png figures:
#                    # use the first word in the 'name' variable.
#                    plt.savefig( f"{Y[1].split(' ')[0]}_{crm}_{cls.replace('/', '')}.png")
#                elif microsoft=='showfull':
#                    plt.get_current_fig_manager().full_screen_toggle()  # toggle fullscreen mode
#                    plt.show()
#                elif microsoft=='showfig':
#                    plt.show()




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
import numpy as np
import repeatmasker_parser3 as rm3
from collections import defaultdict

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


def matplotlib_abline(slope, intercept):
    """ Draws an abline with some determined slope and intercept
    into a matplotlib plot. 
    """
    axes = plt.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, '--')


def scatterplot_repeats(x, y, name, axes):
    """ Create a scatterplot figure with matplotlib.
    Inputs: 
     * X and Y lists of points,
     
     * Whether you'd like to create a png file or
    show the plot directly without saving. 
    
     * Name of the analysis (length or amount?)
    """

    # Calcula mitjana
    mean = sum(y) / len(y)
    # Calcula percentils sup. i inf.
    # divideix entre dos el percent (meitat superior, meitat inferior).
    sup_cutoff = (len(y)/100) * (percent/2)
    inf_cutoff = (len(y)/100) * (100-(percent/2))
    sortedY = [] + y
    sortedY.sort()
    sortedX = [] + x
    sortedX.sort()
    # intercepts amb l'eix Y de les línies de percentil:
    sup_line = sortedY[int(sup_cutoff)]
    inf_line = sortedY[int(inf_cutoff)]
    
    # Dibuixa les tres línies horitzontals:
    for intercept in (sup_line, mean, inf_line):
        x_vals = np.array((0, sortedX[-1]))
        y_vals = intercept + 0 * x_vals
        axes.plot(x_vals, y_vals, '--')

    # Si hi han pocs punts, fes línia
    if len(x) < 100:
        axes.plot(x, y, '.-')  # '.-' is the format string.
    # Si hi han molts punts, scatterplot:
    else:
        axes.plot(x, y, '.')  # '.' is the format string.

    axes.set_title(f"Chromosome: {crm}, Repeat: {cls.replace('/', '')}\nWindow Size: {window_size}, Space between Wndws: {jump}")
    axes.set_xlabel('Posició dins el cromosoma (b)')
    axes.set_ylabel(name)


# Funció (refid, fileRM.out, val_finestra ) --> dict{key: refid, val: where_reps}

# where_reps = {'0--val_finestra': amount, 'v_f--2*v_f': amount, ...}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Si es crida com a script:
if __name__ == '__main__':
    ### VALORS DELS ARGUMENTS VIA TERMINAL:
    # Instruccions respecte els arguments necessaris per cridar l'script:
    try:
        rmfile, window_size, microsoft = (sys.argv[1], int(sys.argv[2]), sys.argv[3] )
    except (ValueError, IndexError):
        print("\nCrit script: ",
                "\n%s <RepeatMasker.out> <window-size> <Show-or-Save?> [<space-between-windows>] [<percent_cutoffs>] \n" % sys.argv[0],
                "\n  * Show-or-Save: do you want the figure to be saved or shown directly with matplotlib?\n",
                "       - savefig: save the figure to the current directory.\n",
                "       - showfig: display the figure directly. \n",
                "       - showfull: display the figure in full-screen mode\n")
        sys.exit(1)
    # Optional Jump between windows:
    try:
        jump = int(sys.argv[4])
    except:
        # default: no jump.
        jump = 0
    # Setting Percent 'ablines'
    try:
        percent = int(sys.argv[5])
    except:
        # default: 2.5% superior and inferior cutoffs.
        percent = 5

    # switch si estas dins a WSL ja que no 
    # te accés a la UI de plt.show():
    if microsoft=='savefig':
        import matplotlib
        matplotlib.use('Agg') # no UI backend in WSL. 
    import matplotlib.pyplot as plt

    # create a default dict which stores values
    #         per refid,           per window,          per repeatclass.  length & amount.
    results = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda :{'l': 0, 'a': 0}) ))

    # permet fer el següent sense necessitar generar prèviament claus:
    # results['refid']['class'][position] += 1
    # (això és degut a que inicialitza amb valor = 0)


    # si la llista es troba buida
    # (no s'ha trobat cap cromosoma dins la llista unique_refids):
    if unique_main_refid(rmfile) == []:
        sys.exit("No s'ha trobat cap cromosoma principal (amb etiqueta 'chr' a la columna refid.) dins el fitxer RepeatMasker.")

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
                        #~results[crm][repobj.repclass[-1]][finestra]['a'] += 1
                        #~results[crm][repobj.repclass[-1]][finestra]['l'] += repobj.replen
                        # Per un overview general de classes:
                        results[crm][repobj.repclass[0]][finestra]['a'] += 1
                        results[crm][repobj.repclass[0]][finestra]['l'] += repobj.replen 

    # create a plot.
    # for each chr:
        # sum all repeat values
        # sum the last chr length to the X axis.
        # plot them
    # show plot. 

    ### SCATTER-PLOT OF WHOLE GENOME, SUMMED FAMILIES...
    
    ### SCATTER-PLOT OF LENGTH:
    # initialize vars to capture results we're interested in...
    last_x_pos = 0

    for crm, crm_val in results.items():
        for cls, cls_val in crm_val.items():
            x_pos = []
            y_len = []
            y_amt = []
            for pos, val in cls_val.items():
                # Recopilar valors:
                x_pos += [int(pos)]
                y_len += [int(val['a'])]
                y_amt += [int(val['l'])]

            j = []
            for i in x_pos:
                j += [i + last_x_pos]
            last_x_pos = int(j[-1])

            scatterplot_repeats(j, y_len,'Length of repeats found in window' , axL)
            scatterplot_repeats(j, y_amt, 'Amount of repeats found in window', axA)
    
    plt.show()

    
    ### SCATTER-PLOT OF MAIN CHRS AND MAIN CLASSES:
    for crm, crm_val in results.items():
        for cls, cls_val in crm_val.items():
            # initialize vars:
            X_posicions=[]
            Y_amount=[]
            Y_length=[]
            for pos, val in cls_val.items():
                # Recopilar valors:
                X_posicions += [pos]
                Y_amount += [int(val['a'])]
                Y_length += [int(val['l'])]

            for Y in ((Y_amount, 'Amount of repeats'), (Y_length, 'Length of repeats found in window')):
                # Y[0] -> List of data.
                # Y[1] -> Name of data.
                fig, ax = plt.subplots()
                scatterplot_repeats(X_posicions, Y[0], Y[1], ax)
                # la funció replace() elimina les barres ("/") que confonen amb directoris 
                if microsoft=='savefig':
                    # to name the saved png figures:
                    # use the first word in the 'name' variable.
                    plt.savefig( f"{Y[1].split(' ')[0]}_{crm}_{cls.replace('/', '')}.png")
                elif microsoft=='showfull':
                    plt.get_current_fig_manager().full_screen_toggle()  # toggle fullscreen mode
                    plt.show()
                elif microsoft=='showfig':
                    plt.show()




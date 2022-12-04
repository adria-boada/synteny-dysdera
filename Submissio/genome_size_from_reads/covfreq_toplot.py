#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# covfreq_toplot.py
#
# 03 de des. 2022  <adria@molevol-OptiPlex-9020>

"""
Creates an interactive matplotlib plot with multiple csv dataframes of paired
coverage and frequency values.

The naming process is unstable and prone to bugs and errors; naming should be
handled through the file header instead of the filename (processing) itself.

Not ready to run with Hercules, it's better to run the script manually, once the
.csv files have already been created.
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Si es crida com a script:
if __name__ == '__main__':

    ### PARSE ARGUMENTS ###
    files = []
    for f in sys.argv[1:]:
        if '.csv' not in f:
            sys.exit('A file without .csv extension was given as argv.')
        else:
            files += [f]
    print('-- A plot will be created from the following files: --')
    [print(f) for f in files]

    ### CREATE A PLOT ###
    # initialize the dataframe:
    df=pd.DataFrame(index=pd.read_csv(files[0]).index)

    print('\n-- The plot will be tagged with the following names: --')
    for f in files:
        # Create tags for each file to be placed in their corresponding
        # plotted line in the histogram.
        # Cut out the '_covg_histogram.csv' extension of the file.
        name = f.split('_covg')[0]
        # Cut the pathname (until last '/')
        name = name.split('/')[-1]
        # Replace all underscores with spaces.
        name = name.replace('_', ' ')
        print(name)
        # A la matriu inicialitzada, fes pd.merge() de cada fitxer csv.
        # Llegeix els arxius csv amb read_csv()
        # Assigna names=[...] a les columnes
        # Elimina la primera fila i queda't sols amb la segona columna gràcies a .iloc[1:,1]
        # Que el merge sigui entre l'index del df anterior i
        #     el df que has llegit amb read_csv() gràcies a left/right_index=True
        df = df.merge(pd.read_csv(f, names=['cov', name]).iloc[1:,1],
                      left_index=True, right_index=True)
    # transforma els valors del df. a int
    df = df.astype('int')
    print('\n-- The df has the following shape: --')
    print(df)
    df.plot()
    plt.show()


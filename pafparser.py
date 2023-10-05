#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# .paf_parser.py
#
# 27 abr 2022  <adria@molevol-OptiPlex-9020>

"""Parse a PAF file. Uses a class to store values inside a Python object.
Transforms a PAF file into a more manageable Python3 object.

Afterwards, create a pandas dataframe from the refined and translated Python
objects.

--30-09-22--
Intenció: emprar 'pandas' per gestionar les anàlisis. Crear un dataframe similar
al PAF cru, que difereixi en el seu refinament exquisit.
  1. class Mapping(): eat a raw line from a PAF file and spit out a refined obj.
  2. parse file: create a list of refined obj, one for each line.
  3. get results: print out data thanks to pandas from the created dataframe
  (mean, max, plots...).

  Pandas can eat a list of dictionaries with keys as column names, and values as
cells. Each dictionary in the list is a row.

  Def function_parser(): open the file. for each line, create a Mapping() obj.
Append these obj to a list. From it, create a dataframe. Now: which columns am
I interested in? Two chromosome tags, query and target. Ali. length in query and
target (end - start). BLAST-id and compressed-id Number of bases hit (col11).
Matches (col10). Mismatches (NM tag). Add CIGAR matches, so you can compare
between them and see if there is any errors. Deletion should mean gap in the
target sequence. Insertion should mean gap in the query sequence. Map quality.

  I'd be rather interested into creating a script which iterates over ranges and
finds the non-overlapping ranges. Then we could compute the real coverage for chr A.
"""

import sys

# Dataframes manipulation, akin to the R project
import pandas as pd
import numpy as np
from natsort import index_natsorted
import seaborn as sns

# Plotting data into barplots, histograms, etc.
import matplotlib.pyplot as plt

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def round_to_Mbps(number:"Amount of bases", decims:"Nº of trailing decimals"=2):
    """Transforma una mesura qualsevol al prefixe d'unitats `Mega`, que denota
    un factor de 10⁶. Mostra fins a `decims` decimals a l'hora d'arrodonir.
    Facilita la visualització de les dades, tot i que la funció és ximplona i
    tant sols és un embolcall de la funció `round`.
    """
    return round( number/(10**6), decims)

class Mapping(object):
    """ Parse fields given a line from a PAF formatted file. Creates an
    instance of the "Mapping" class from a single line of a PAF file.
    """
    def cig_analysis (self, cig: "CIGAR string"):
        """
        INPUT
        -----

        A cigar string: A string of letters symbolizing the composition of an
        alignment. Accepts a string with 'M', 'I', 'D' (and no other type).

        OUTPUT
        ------

        Returns a dictionary with total length for matches, insertions and deletions
        ('M', 'I', 'D'). It also returns the amount of insertions and deletions as a
        key named 'compressed' (good for calculating compressed sequence identity;
        refer to
        <https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity>
        for an in-depth explanation of the 'compressed-identity' concept).
        """

        # Split string by adding ' ' to the end of any letter:
        for i in 'MID':
            cig = cig.replace(i, i+' ')

        # Transform spaces into list separations:
        c = cig.strip().split(' ')

        # Prepara diccionari pel retorn:
        answer = {
                'M': 0,
                'I': 0,
                'D': 0,
                'compressed': 0
                }

        # Recompte i guarda dins 'answer':
        for i in c:
            if 'M' in i:
                answer['M'] += int(i[:-1])
            elif 'D' in i:
                answer['D'] += int(i[:-1])
            elif 'I' in i:
                answer['I'] += int(i[:-1])

        # Tingues en compte el nombre de gaps i no tant la seva llargada
        # (comprimeix-los). Per info revisa sota la def. de la funció.
        answer['compressed'] = cig.count('D') + cig.count('I')

        return answer

    def __init__(self, line):
        """Create Python variables from the given file columns...

        INPUT
        -----

        A line from a PAF formatted file.

        OUTPUT
        ------

        A Python object of the Mapping() class, hopefully more manageable and
        refined than the raw data.

        Small documentation of alignment measures:
         * Column 10 (self.matches) is the number of exactly matching bases in
         the alignment.

         * Column 11 (self.ali_len) are the total number of bases forming part of the
         alignment (including gaps, matches and no-matches).

         * mapQ (self.mapQ) scores how likely is it that a mapping is correctly aligned.
         Low scores mean that there are many possible mappings for the
         read/query inside the targeted database. High scores mean fewer ways in
         which the read/query could be mapped to the targeted database.

         * Non-Matching (self.no_match) is the number of incorrect matches in
         the alignment (the sum of mismatches and gaps).

         * Ambiguous bases (self.ambiguous) are bases left as "N" inside the
         alignment.

        In summary:
        no_match + match = ali_len = cigM + cigI + cigD
        """
        (self.qname,        # query name
            self.qlen,      # query length
            self.qstart,    # query start coord
            self.qend,      # query end coord
            self.strand,    # strand ( + or - )
            self.tname,     # same for target...
            self.tlen,
            self.tstart,
            self.tend,
            self.matches,   # number of matching bases
            self.ali_len,   # nº of matches+ misses+ gaps
            self.mapQ       # mapping quality
        ) = line.split()[:12]

        # int-ize numeric values...
        self.qlen = int(self.qlen)
        self.qstart = int(self.qstart)
        self.qend = int(self.qend)
        self.tlen = int(self.tlen)
        self.tstart = int(self.tstart)
        self.tend = int(self.tend)
        self.matches = int(self.matches)
        self.ali_len = int(self.ali_len)
        self.mapQ = int(self.mapQ)

        # Detect whether the given scaffold is a chromosome or a minor
        # unassembled contig. If it is the later, change its name to a generic
        # tag (`Scaff`) which groups all of them under a single name: 
        if 'scaffold' in self.qname.lower() or 'ctg' in self.qname.lower():
            self.qname = 'Scf'
        if 'scaffold' in self.tname.lower() or 'ctg' in self.tname.lower():
            self.tname = 'Scf'
        # `Scaff` results will be a mean/sum of all minor contigs. The lower
        # function makes the string lowercase (case-insensitive matching).

        # add a small tag to visually distinguish query and target when printed to
        # the terminal:
        self.qname = "Q." + self.qname
        self.tname = "T." + self.tname

        # Find and add the other lines we're interested in (some measures might
        # not always be present in the PAF-file... misterious?):
        for i in line.split()[12:]:
            if 'NM:i:' in i:
                # sum of mismatches and gaps in alig.
                self.no_match = int(i[5:])
            elif 'tp:A:' in i:
                # tipus d'alineament; 1ari o 2ari.
                self.typeA = i
            elif 'nn:i:' in i:
                # Bases ambigues (NNN)
                self.ambiguous = int(i[5:])
            elif 'cg:Z:' in i:
                # CIGAR-string.
                # the [5:] slicing strips beggining 'cg:Z:'
                # the strip() method removes ending '\n'.
                self.cigar = i[5:].strip("\n")
                # analyze the amount of matches and indels in the `CIGAR` string
                # refer to the above function "cig_analysis" to see the keys of
                # the returned dictionary.
                self.cig_summary = self.cig_analysis(self.cigar)

#            # No l'he aconseguit fer funcionar ......
#            elif 'dv:' in i:
#                # divergencia (si es troba)
#                self.divergence = i[5:]
#            elif 'de:' in i:
#                # divergencia gap-compressed (si es troba)
#                self.comp_div = i[5:]

        # Time for the 'pandas' dict definition.
        self.df = {'Qname': self.qname,
                   'Qstart': self.qstart,
                   'Qend': self.qend,
                   'Qlen': (self.qend - self.qstart),
                   'Tname': self.tname,
                   'Tstart': self.tstart,
                   'Tend': self.tend,
                   'Tlen': (self.tend - self.tstart),
                   'Align_len': self.ali_len,
                   'BLAST-id': (self.matches/self.ali_len),
                   'Compressed-id': (self.matches /
                                      (self.cig_summary['M']+self.cig_summary['compressed'])),
                   'MapQ': self.mapQ,
                   'Matches': self.matches,
                   'No_matches': self.no_match,
                   # misses = NoMatch - gaps
                   'Mismatches': (self.no_match -
                                  (self.cig_summary['D']+self.cig_summary['I'])),
                   'Deletions': self.cig_summary['D'],
                   'Insertions': self.cig_summary['I'],
                   'cig_Matches': self.cig_summary['M'],
                   'Ambiguous': self.ambiguous,
                   'Strand': self.strand,
                   'TypeAl': self.typeA,
                   }

def parse_paf_alignment_db (filename: "Path to PAF formatted file"):
    """Returns a list of 'Mapping' classes.

    colnames are very similar to the ones chosen by the manual of minimap2.
    colnames are explained in detail inside the 'Mapping' class.

    Returns a tuple with (wholly-dataframe, unique-list-of-scaffolds)
    Select one or the other with tuple-slicing (i.e. func(input)[1])
    """

    # Parsing paf-file into Mapping() Python3 object:
    paf = []

    # LLista de cromosomes i contigs únics, emparellats amb la seva llargada:
    unique_scaff = {'Q.Scf': 0, 'T.Scf': 0}

    # In theory, could open both gzipped or normal files:
    with gzip.open(filename) if filename.endswith('.gz') else open(filename) as fn:

        for line in fn:
            # Use the mapping class to create an object from 'line'.
            x = Mapping(line)

            # Apend it to 'paf' return list.
            paf.append(x.df)

            # Si és el primer cop que trobem el cromosoma/contig:
            if x.qname == 'Q.Scf':
                unique_scaff['Q.Scf'] += x.qlen
            elif not x.qname in unique_scaff.keys():
                unique_scaff[x.qname] = x.qlen
            if x.tname == 'T.Scf':
                unique_scaff['T.Scf'] += x.tlen
            elif not x.tname in unique_scaff.keys():
                unique_scaff[x.tname] = x.tlen

    return paf, unique_scaff

def chromosomal_coverage(
        intervals):
    """
    INPUT

    * intervals
    A list of sublists where each sublist is composed of two coordinates (begin
    and end) and the chromosome of the mapping. Eg:
    [ [begin1, end1, seq1], [begin2, end2, seq2], ..., [bN, eN, sN] ]

    Submit a single sequid (not multiple sequids) !!!

    OUTPUT

    Sum of non-overlapping mapping lengths. Removes all mapping subsegments
    which are in overlapping and would be counted twice toward the sum.
    """
    # Create a point list from an interval list
    points = []
    for i in range(0, len(intervals)):
        # Left/starting points labeled "0"
        points += [[intervals[i][0], 0, i]]
        # Right/ending points labeled "1"
        points += [[intervals[i][1], 1, i]]
    # Sort points by position (by first value in sublists)
    points.sort()

    # Init variables
    # there should be only one sequid in the whole intervals list
    sequid = intervals[0][2]
    d_return = {sequid: 0}
    size_distrib = {"0-9 bp": [0,0], "10-99 bp": [0,0], "100-999 bp": [0,0],
                    "1-9,999 kb": [0,0], "10-99,999 kb": [0,0],
                    "100-999,999 kb": [0,0], "1-10 Mb": [0,0], "higher": [0,0]}
    # keep track of open and closed intervals
    currentOpen = -1

    # for each point in the list
    for i in range(0, len(points)):

        # If the loop landed on a left point (0)
        # which opens an interval:
        if points[i][1] == 0:
            # if there is no other interval opened:
            if currentOpen == -1:
                # enters interval "i"
                currentOpen = points[i][2]
                currentBegin = int(points[i][0])
                # print("enters", currentOpen) # debug
            # else, there already was an open interval:
            else:
                # from the two mappings that are open, find which has the higher
                # end coordinate and make it the currentOpen interval.
                currentEnd = intervals[currentOpen][1]
                nextEnd = intervals[points[i][2]][1]
                if nextEnd > currentEnd:
                    currentOpen = points[i][2]

        # If the loop landed on a right point (1)
        # which closes an interval:
        else:
            # And it is the right point of the "currentOpen" interval...
            if points[i][2] == currentOpen:
                # compute between begin and end
                end = points[i][0]
                mapping_bps = end - currentBegin
                d_return[sequid] += int(mapping_bps)
                # Size distribution of synteny blocks
                if mapping_bps < 10:
                    size_distrib["0-9 bp"][0] += 1
                    size_distrib["0-9 bp"][1] += mapping_bps
                elif mapping_bps < 100:
                    size_distrib["10-99 bp"][0] += 1
                    size_distrib["10-99 bp"][1] += mapping_bps
                elif mapping_bps < 1000:
                    size_distrib["100-999 bp"][0] += 1
                    size_distrib["100-999 bp"][1] += mapping_bps
                elif mapping_bps < 10000:
                    size_distrib["1-9,999 kb"][0] += 1
                    size_distrib["1-9,999 kb"][1] += mapping_bps
                elif mapping_bps < 100000:
                    size_distrib["10-99,999 kb"][0] += 1
                    size_distrib["10-99,999 kb"][1] += mapping_bps
                elif mapping_bps < 1000000:
                    size_distrib["100-999,999 kb"][0] += 1
                    size_distrib["100-999,999 kb"][1] += mapping_bps
                elif mapping_bps < 10000000:
                    size_distrib["1-10 Mb"][0] += 1
                    size_distrib["1-10 Mb"][1] += mapping_bps
                else:
                    size_distrib["higher"][0] += 1
                    size_distrib["higher"][1] += mapping_bps
                # print("adds", mapping_bps) # debug
                # close currentOpen interval
                currentOpen = -1

    return d_return, size_distrib

def genome_coverage(df):
    """
    Revisa el coverage de cada cromosoma; com de ben coberts es troben els crms.
    També avalua la freqüència de blocs sintènics agrupats dins una certa llargada.

    No és una vertadera funció, és un monstre lleig.
    """
    print("# Printing the size distribution of syntenic blocks")
    print("# The two values in dict lists are num of mappings and bps mapped",
          "in each range\n")
    distrib_list = []
    size_distrib = {"0-9 bp": [0,0], "10-99 bp": [0,0], "100-999 bp": [0,0],
                    "1-9,999 kb": [0,0], "10-99,999 kb": [0,0],
                    "100-999,999 kb": [0,0], "1-10 Mb": [0,0], "higher": [0,0]}
    for sequid in df["Qname"].unique():
        if 'chr' in sequid.casefold():
            dfsel = df.loc[df["Qname"] == sequid]
            begins = list(dfsel["Qstart"])
            ends = list(dfsel["Qend"])
            seqs = [sequid] * len(begins)
            intervals = list(zip(begins, ends, seqs))
            result = chromosomal_coverage(intervals)
            coverage = result[0][sequid]
            distrib_synteny_blocks = result[1]
            distrib_list.append(distrib_synteny_blocks)
            print("*", sequid+": covered as much as", coverage, "unique bp in "+
                  "the chr.")
            print("    ", distrib_synteny_blocks)
    for i in distrib_list:
        for key in i.keys():
            size_distrib[key][0] += i[key][0]
            size_distrib[key][1] += i[key][1]
    print(" ** Sum size distrib. for query:", size_distrib)
    print()
    distrib_list = []
    size_distrib = {"0-9 bp": [0,0], "10-99 bp": [0,0], "100-999 bp": [0,0],
                    "1-9,999 kb": [0,0], "10-99,999 kb": [0,0],
                    "100-999,999 kb": [0,0], "1-10 Mb": [0,0], "higher": [0,0]}
    for sequid in df["Tname"].unique():
        if 'chr' in sequid.casefold():
            dfsel = df.loc[df["Tname"] == sequid]
            begins = list(dfsel["Tstart"])
            ends = list(dfsel["Tend"])
            seqs = [sequid] * len(begins)
            intervals = list(zip(begins, ends, seqs))
            result = chromosomal_coverage(intervals)
            coverage = result[0][sequid]
            distrib_synteny_blocks = result[1]
            distrib_list.append(distrib_synteny_blocks)
            print("*", sequid+": covered as much as", coverage, "unique bp in "+
                  "the chr.")
            print("    ", distrib_synteny_blocks)
    for i in distrib_list:
        for key in i.keys():
            size_distrib[key][0] += i[key][0]
            size_distrib[key][1] += i[key][1]
    print(" ** Sum size distrib. for target:", size_distrib)
    print()

def histogram_indels(
    df, write_to_filepath, title,
    logscale="False",
    ):
    """
    Draw a distribution of indel sizes from a mapping dataset.

    * df: DataFrame of mapping data.
    colname with insertion values: 'Insertions'
    colname with deletion values: 'Deletions'

    * write_to_filepath: save PNG figure to this path

    * title: the title of the plot.

    * logscale: change X axis to logarithmic scale.
    """
    # specify amount of bins (as many as max value?)
    #~ amount_bins = math.floor(df["perc_divg"].max()) + 1
    # Create the matplotlib axes (a row of three consecutive plots)
    fig, (ax1, ax2, ax3,) = plt.subplots(
            nrows=1, ncols=3, figsize=(18, 6))
    # The first two axes are absolute Ins and Del; the third is a relative
    # comparison of the preceding data (stat='percent')
    ax1.set_title('Deletions')
    ax2.set_title('Insertions')
    ax3.set_title('Rel. comparison')
    # Iterate through multiple seaborn axes and give to each of them the
    # adequate ax=axN number.
    # Element "poly" changes the style of the plot (see seaborn doc).
    sns.histplot(data=df, x='Deletions', ax=ax1, element="step")
    sns.histplot(data=df, x='Insertions', ax=ax2, element="step")
    # Prepara dades pel tercer plot relatiu/comparatiu...
    # Per a que pd.concat funcioni, anomena igual les columnes dels dataframes
    # df_ins i df_del (e.g. anomena les dues sèries "size").
    df_del = df['Deletions'].rename("size").to_frame()
    df_ins = df['Insertions'].rename("size").to_frame()
    ## print(df_del) ## DEBUG
    ## print(df_ins) ## DEBUG
    # Etiqueta cada valor per aplicar 'hue'.
    df_ins['hue'] = 'Ins'
    df_del['hue'] = 'Del'
    ## print(df_del) ## DEBUG
    ## print(df_ins) ## DEBUG
    # Empra 'pd.concat()' per a concatenar les dues sèries anteriors (una sola
    # columna amb totes les dades).
    df_concat = pd.concat([df_ins, df_del])
    # Crea el tercer histograma, aquesta vegada relatiu amb totes les dades.
    ## print(df_concat) ## DEBUG
    sns.histplot(data=df_concat, x='size', ax=ax3,
                 multiple='fill',
                 hue='hue', palette=['C1', 'C2'])

    if logscale:
        ax1.set_xscale("log")
        ax2.set_xscale("log")
        ax3.set_xscale("log")
    plt.title(title, fontsize=10)
    # plt.suptitle("Divergence from consensus sequence", fontsize=12)
    plt.xlabel("Categories, in basepairs")
    plt.ylabel("Frequency/Density")
    plt.savefig(write_to_filepath, dpi=300)
    plt.close('all')

    return None

# folder setting    if folder:
# folder setting        if not os.path.exists(folder):
# folder setting            os.makedirs(folder)
# folder setting        write_to_filepath = (folder+"/divhist_"+
# folder setting            species+"_"+"_".join(reptype)+'.png')
# folder setting    else:
# folder setting        write_to_filepath = ("divhist_"+
# folder setting            species+"_"+"_".join(reptype)+'.png')

if __name__ == "__main__":

    # Instruccions respecte els arguments necessaris per cridar l'script:
    if len(sys.argv) < 2:
        sys.exit('\nCrit script: script.py <input-paf-file.paf>\n')

    # Define paf-file
    paf_file = sys.argv[1]

    # List of unique scaffolds and list of 'rows':
    list_paf_rows, unique_scaff = parse_paf_alignment_db(paf_file)
    #print(unique_scaff) ##DEBUG
    #print(list_paf_rows[0]) ##DEBUG
    #print(list_paf_rows[1].df) ##DEBUG

    # Create a mighty dataframe.
    df = pd.DataFrame(list_paf_rows)

    # Options to fully view and print tables:
    pd.set_option('display.max_columns', None)

    # Revisa coverage del genoma:
    genome_coverage(df)

    print("# Printing first five rows of the PAF file...")
    print("----------")
    print(df.head(5), end='\n\n')

    print("# Printing information about the PAF file...")
    print("----------")
    print("Number of rows in file (num. mappings):", df.shape[0])
    for type_alignment in df["TypeAl"].unique():
        val = df.loc[df["TypeAl"] == type_alignment].shape[0]
        print(f"Number of rows with {type_alignment} type:", val)
    print()

    # Descripció d'estadístiques bàsiques.
    print("# Printing a description of the 'dataset' (PAF file).")
    print("----------")
    print(df.describe(), end='\n\n')

    print("# Amount of rows for each chromosome tag...")
    print("----------")
    print(df['Qname'].value_counts(normalize=True), end='\n\n')
    print(df['Tname'].value_counts(normalize=True), end='\n\n')

    ## Intent histograma ##

    histogram_indels(df, "histogram_indels_linearscale.png", "Standard", logscale=False)
    histogram_indels(df, "histogram_indels_logscale.png", "Logscale")

    ## BoxPlotting section ##

    # Create a concatenate of query and target values
    # (mapqual, ali-len, etc.)
    query = df.loc[:, ['Qname', 'Align_len', 'BLAST-id', 'Compressed-id', 'MapQ']]
    query['Indel'] = df['Deletions'] + df['Insertions']
    query.rename(columns={'Qname': 'Sequid'}, inplace=True)
    target = df.loc[:, ['Tname', 'Align_len', 'BLAST-id', 'Compressed-id', 'MapQ']]
    target['Indel'] = df['Deletions'] + df['Insertions']
    target.rename(columns={'Tname': 'Sequid'}, inplace=True)
    df_concat = pd.concat([query, target])

    # natsort df_concat by sequid
    # (removing 'Q.' or 'T.', sorting by 'Dcat' or 'Dtil')
    df_concat = df_concat.sort_values(
                            by=['Sequid'],
                            key=lambda x: np.argsort(index_natsorted(
                            # sort by "Sequid", removing first two
                            # characters from string in column
                            df_concat["Sequid"].str[2:])))

    fig_titles = {
        # 'df_colname': 'Fig. title',
        'Align_len': 'Length of sequences alignment (xlim cut)',
        'BLAST-id': 'Standard/BLAST sequences identity',
        'Compressed-id': 'Gap-compressed sequences identity',
        'Indel': 'Indels/gaps in sequences (xlim cut)',
        'MapQ': 'Mapping quality'}
    fig_xlabel = {
        # 'df_colname': 'Fig. xlabel',
        'Align_len': 'Basepairs',
        'BLAST-id': 'Percent identity',
        'Compressed-id': 'Percent identity',
        'Indel': 'Basepairs',
        'MapQ': 'Mapping quality'}
    for xvar in ['Indel',
            'Align_len', 'BLAST-id', 'Compressed-id', 'MapQ']:
        sns.boxplot(data=df_concat, x=xvar, y='Sequid',
                    showmeans=True,
                    meanprops=dict(color='green',
                                   marker='*',
                                   markeredgecolor='black',
                                   markersize='15'))
        # Els rangs de l'eix X en alguns casos és massa ample
        if xvar == 'Align_len':
            plt.xlim(0, 10000) # Restringeix xlim
        elif xvar == 'Indel':
            plt.xlim(0, 2000) # Restringeix xlim
        plt.title(fig_titles[xvar])
        plt.subplots_adjust(left=0.2, bottom=0.1, right=0.95, top=0.91)
        plt.xlabel(fig_xlabel[xvar])
        print("# Plotting", str(xvar)+'_matplotlib.png')
        print("----------")
        plt.savefig(str(xvar)+'_matplotlib.png', dpi=300)
        plt.close('all')



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
target (end - start). BLAST-id and compressed-id. Number of bases hit (col11).
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
import natsort

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
        if 'scaffold' in self.tname.lower() or 'ctg' in self.qname.lower():
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
                   'Alig. len.': self.ali_len,
                   'BLAST-id.': (self.matches/self.ali_len),
                   'Compressed-id.': (self.matches /
                                      (self.cig_summary['M']+self.cig_summary['compressed'])),
                   'MapQ.': self.mapQ,
                   'Matches': self.matches,
                   'No-matches': self.no_match,
                   # misses = NoMatch - gaps
                   'Mismatches': (self.no_match -
                                  (self.cig_summary['D']+self.cig_summary['I'])),
                   'Deletions': self.cig_summary['D'],
                   'Insertions': self.cig_summary['I'],
                   'cig.Matches': self.cig_summary['M'],
                   'Ambiguous': self.ambiguous,
                   'Strand': self.strand,
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

    # Veure la capçalera.
    print("# Printing first five rows of the PAF file...")
    print("----------")
    print(df.head(5), end='\n\n')

    # Descripció d'estadístiques bàsiques.
    print("# Printing a description of the 'dataset' (PAF file).")
    print("----------")
    print(df.describe(), end='\n\n')

    print("# Amount of rows for each chromosome tag...")
    print("----------")
    print(df['Qname'].value_counts(normalize=True), end='\n\n')
    print(df['Tname'].value_counts(normalize=True), end='\n\n')

    ## Plotting section ##

    # Box plot of mapping quality for each chromosome.
    # Create a dataframe in which each column is the mapping qualities of a
    # single chromosome/scaffold.
    df_mapq = pd.DataFrame()
    for name in ['Tname', 'Qname']:
        for crm in df[name].unique():
            df_crm_mapq = df[df[name]==crm].reset_index()['MapQ.']
            df_crm_mapq = df_crm_mapq.rename(crm)
            df_mapq = pd.concat([df_mapq, df_crm_mapq], axis=1)

    # natsort the df_mapq
    df_mapq = df_mapq.reindex(columns=natsort.humansorted(df_mapq.columns))

    print("# Writing mapqual-boxplot.jpg...")
    print("----------")
    # Plot `df_mapq` with matplotlib.
    df_mapq.plot.box(showmeans=True,
                     #figsize=(20,10),
        )
    plt.xticks(rotation = 30) # Rotates X-Axis Ticks by 30-degrees
    plt.title('Mapping qualities for individual chromosomes')
    plt.savefig("mapqual-boxplot.jpg")



#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# .paf_parser.py
#
# 27 abr 2022  <adria@molevol-OptiPlex-9020>

""" Parse a paf-file. Store its values in a class.
Transforms a paf-file into a more manageable
Python3 object.

Module for paf_analyzer2.py

"""

import sys

# Aids in adding a default value of zero to newly created dict-keys (RELICT).
#from collections import defaultdict

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def round_to_Mbps(number:"Amount of bases", decims:"Nº of trailing decimals"=2):
    """Transforma una mesura qualsevol al prefixe d'unitats `Mega`, que denota
    un factor de 10⁶. Mostra fins a `decims` decimals a l'hora d'arrodonir.
    Facilita la visualització de les dades, tot i que la funció és ximplona i
    tant sols és un embolcall de la funció `round`.
    """
    return round( number/(10**6), decims)


def cig_analysis (cig: "CIGAR string"):
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
    # ...Fes .strip() igualment per eliminar espai buit.

    # Elimina un espai buit al final:
    # c = c[:-1]
    # (Ja es sol fer abans d'entrar a la func.)

    # Prepara retorn:
    answer = {
            'M': 0,
            'I': 0,
            'D': 0,
            'compressed': 0
            }

    # Recompte:
    for i in c:
        if 'M' in i:
            answer['M'] += int(i[:-1])
        elif 'D' in i:
            answer['D'] += int(i[:-1])
        elif 'I' in i:
            answer['I'] += int(i[:-1])

    # Tingues en compte el nombre de gaps
    # i no tant la seva llargada (comprimeix-los).
    answer['compressed'] = cig.count('D') + cig.count('I')

    return answer


class Mapping(object):
    """ Parse fields given a line from a PAF formatted file. Creates an
    instance of the "Mapping" class from a single line of a PAF file.
    """
    def __init__(self, line):
        # Create variables from the given file-line. 
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
            self.n_bases,   # nº of matches+ misses+ gaps
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
        self.n_bases = int(self.n_bases)
        self.mapQ = int(self.mapQ)

        # Detect whether the given scaffold is a chromosome or a minor
        # unassembled contig. If it is the later, change its name to a generic
        # tag (`Scaff`) which groups all of them under a single name: 
        if 'Scaffold' in self.qname:
            self.qname = 'Scaff'
        if 'Scaffold' in self.tname:
            self.tname = 'Scaff'
        # `Scaff` results will be a mean/sum of all minor contigs.

        # add a small tag to visually distinguish query and target when printed to
        # the terminal:
        self.qname = "Q." + self.qname
        self.tname = "T." + self.tname

        # find and add the other lines we're
        # interested in (some measures might not always be present???):
        for i in line.split()[12:]:
            if 'NM:i:' in i:
                # sum of mismatches and gaps in alig.
                self.mismatch = int(i[5:])
            elif 'tp:A:' in i:
                # tipus d'alineament; 1ari o 2ari.
                self.typeA = i
            elif 'nn:' in i:
                # Bases ambigues (NNN)
                self.ambiguous = int(i[5:])
            elif 'cg:' in i:
                # CIGAR-string.
                # the [5:] slicing strips beggining 'cg:Z:'
                # the strip() method removes ending '\n'.
                self.cigar = i[5:].strip("\n")
#            elif 'dv:' in i:
#                # divergencia (si es troba)
#                self.divergence = i[5:]
#            elif 'de:' in i:
#                # divergencia gap-compressed (si es troba)
#                self.comp_div = i[5:]


class Results(object):
    def __init__(self):
        """ Creates a dictionary in which results can be stored.

        The init function stores the results template and creates
        the 'self.r' results-storage dictionary-object.

        These results are stored as the sum for all variables
        inside the chromosome under scrutiny.
        """
        # Template to create many 'results' dicts.
#        self.entry_template = {'ll': 0, # llargada
#                'cgMatch': 0,
#                'cgLost': 0,
#                'cgAdded': 0,
#                'cgCompressed': 0,    # amount of gap occurrences
#                'col10': 0,
#                'col11': 0,
#                'NM': 0,
#                'ambiguous': 0,
#                'amount_maps': 0,
#                'type': defaultdict(lambda: 0), # amount of each type of alignment.
#                # to be a defaultdict means that '+=' can be used with
#                # newly created keys without throwing errors.
#                'default_dv': 0,      # divergence
#                'default_de': 0,      # compressed-dv.
#                'number_additions': 0 # sum one for each parsing pass.
#                }
        self.r = {}  # Where all results will be stored


    def new_entry(self, scaff_name, ll_bases):
        """ Creates a new entry for scaffold -> scff_name,
        with length -> ll_bases.

        Only sets this two parameters; others are added while
        the paf-file is parsed.
        """
        # Create entry
        self.r[scaff_name] = {'ll': ll_bases,
                'cgMatch': 0,
                'cgLost': 0,
                'cgAdded': 0,
                'cgCompressed': 0,    # amount of gap occurrences
                'col10': 0,
                'col11': 0,
                'NM': 0,
                'ambiguous': 0,
                'amount_maps': 0,
                'type': {'tp:A:P': 0,   # 1ary alignments (one hit)
                         'tp:A:S': 0,   # 2ary alignments (multiple hit)
                         'tp:A:i':0,    # inverted alignments
                         'tp:A:I': 0},  # inverted alignments
                # previously employed a defaultdic... not anymore.
                # employing a defaultdict means that '+=' can be used with
                # newly created keys without throwing errors.
                'default_dv': 0,      # divergence
                'default_de': 0,      # compressed-dv.
                'number_additions': 0 # sum one for each parsing pass.
                }


    def update_entry(self,
            scaff: "scaffold's name where results should be added",
            cgm: "sum of cigarMatches" = 0,
            cgl: "sum of cigarLosses (Dels)" = 0,
            cga: "sum of cigarAdditions (Ins)" = 0,
            cgc: "sum of cigarCompressed (amount of gaps)" = 0,
            col10: "amount of perfect matches" = 0,
            col11: "length of alignment (counting gaps)" = 0,
            NM: "non-matching bases (sum of gaps and misses)" = 0,
            tp: "supply string with type of alignment for the line" = None,
            nn: "ambiguous bases (NNNs)" = 0,
            n_ali: "number of mappings calculated" = 0,
            ll: "add length to this particular entry" = 0,
            dv: "default divergence" = 0,
            de: "default compressed-divergence" = 0,
            ):
        # shorthand for results dictionary for a selected chromosome.
        x = self.r[scaff]
        # add cigarMatches.
        x['cgMatch'] += cgm
        # add cigarLosses (deletions in scaff)
        x['cgLost'] += cgl
        # add cigarAdded (insertions in scaff)
        x['cgAdded'] += cga
        # add cigarCompressed (amount of indels)
        x['cgCompressed'] += cgc
        # add total *CORRECTLY MATCHING* bases in mapping.
        x['col10'] += col10
        # add total bases in this mapping (line)
        x['col11'] += col11
        # Add No-Matches (mismatches and gaps in alig.)
        x['NM'] += NM
        # type of alignment?
        # if `tp` is not the 'None' given by default...
        if tp:
            # add one to such type.
            x['type'][tp] += 1
        # ambiguous bases (NNNs)
        x['ambiguous'] += nn
        # number of mappings for each crm:
        x['amount_maps'] += n_ali
        # llargada
        x['ll'] += ll


    def print_header(self):
        """ Print the header in the terminal with the
        print() function. These results are then copy-pasted
        into a new file.
        """
        # Start by printing the headers of the CSV:
        print ("Seq-Name",             # Scaffold name.
            "# Llargada CHR",
            "# Bases Map. (w/gaps)", # col11, number of bases (including gaps) in alignment
            "# Bases Map. (wo/gps)", # cgMatches, number of bases (excluding gaps) in alig.
            "% Covg.",             # Coverage; cgMatches/total_len
            "# matches",
            "# misses",
            "# «Added»",           #added seq in crm
            "# «Lost»",            #lost seq in crm
            "# Indel Sum",
            "% ~BLAST-id",         # col10/col11
            # BLAST-id is the ratio of exact matches versus all other bases in mapping
            # (matches versus matches, mismatches and gaps)
            "% Gap-comprssd-id",   # same but counts each gap as a single mismatch.
            sep=";"
            )

    def print_results(self, crm):
        """ Print the results for a
        requested chromosome with the print()
        function. These results are then copy-pasted
        into a new file.
        """
        # Shortens variable calls:
        x = self.r[crm]

        # Compute percentual variables:

        # matches throughout the chromosome.
        covg = round(
                    (x['cgMatch'] / x['ll'])*100,
                4 )
        # identity: correct vs incorrect matches
        blast_id = round(
                    (x['col10'] / x['col11'])*100,
                    4 )
        # identity: correct vs incorrect matches
        gap_comp = round(
                    (x['col10'] / (x['cgMatch'] + x['cgCompressed']) )*100,
                    4 )
        # Compute sum of gaps and only misses:
        sum_of_gaps = int( x['cgAdded'] + x['cgLost'] )
        # Misses is 'NM'-'tot_gaps' or 'cg_Matches' - 'col10'
        misses = int( x['NM'] ) - sum_of_gaps


        print(crm,                              # crm name
                round_to_Mbps(x['ll']),         # llargada crm
                round_to_Mbps(x['col11']),      # col11
                round_to_Mbps(x['cgMatch']),    # cgMatches
                f"{covg} %",                    # covered scaffold by the others
                round_to_Mbps(x['col10']),      # only matches
                round_to_Mbps(misses),# only misses
                round_to_Mbps(x['cgAdded']),# extra seqs in crm
                round_to_Mbps(x['cgLost']),# missing seqs in crm
                round_to_Mbps(sum_of_gaps),# sum of indels
                f"{blast_id} %",# blast-id
                f"{gap_comp} %",# Gap-comprssd-id
                sep=";")


def parse_paf_alignment_db (filename: "Path to PAF formatted file"):
    """Returns a list of 'Mapping' classes.

    To call the newly created list of Mappings, use:
    map-object = parse_paf_alignment_db(filename)

    In order to print the value of a given column created by the class:
    for i in map:
        print( i.colname )

    colnames are the same as the manual for minimap2.
    colnames are explained in detail inside the
    'Mapping' class.

    Returns a tuple with (parseable-paf, unique-list-of-scaffolds)
    Select one or the other with tuple-slicing (i.e. func(input)[1])
    """

    # Parsing paf-file into Mapping() Python3 object:
    paf = []

    # LLista de cromosomes i contigs únics:
    unique_scaff = []

    # In theory, could open both gzipped or normal files:
    with gzip.open(filename) if filename.endswith('.gz') else open(filename) as fn:

        for line in fn:
            # Use the mapping class to create
            # a 'line' object. Append it to 'paf' list.
            paf.append(Mapping(line))

            # Last Class/line entry is:
            x = paf[-1]

            # Si és el primer cop que trobem el cromosoma/contig:
            if not [x.qname, int(x.qlen)] in unique_scaff:
                print(x.qname, unique_scaff) ##DEBUG
                unique_scaff += [ [x.qname, int(x.qlen)] ]
            if not [x.tname, int(x.tlen)] in unique_scaff:
                unique_scaff += [ [x.tname, int(x.tlen)] ]

    return paf, unique_scaff


if __name__ == "__main__":

    # Instruccions respecte els arguments necessaris per cridar l'script:
    if len(sys.argv) < 2:
        sys.exit('\nCrit script: script.py <input-paf-file.paf>\n')

    # Define paf-file
    paf_file = sys.argv[1]

    # List of unique scaffolds:
    unique_scaff = parse_paf_alignment_db(paf_file)[1]
    #print(unique_scaff) ##DEBUG

    # Create a dict~class in which to store definitive results.
    results = Results()

    # Crea una entrada per a resultats generals...
    # (total alineament, resum alineament)
    # amb llargada inicial de zero.
    results.new_entry('general', 0)

#    # Inclou ja les entrades per a la
#    # mescla d'Scaffolds petits:
#    # (( SI EXISTEIXEN SCAFFOLDS MINORITARIS! ))
#    # Fa una passada pels scaffolds únics i si troba
#    # algun de minoritari, crea una entrada on embutir-los.
#    for crm in unique_scaff:
#        # If you find a Scaffold:
#        if 'Scaffold' in crm[0]:
#            # and it is from the query:
#            if crm[0].startswith('Q.'):
#                results.new_entry('Q.Scaff', 0)
#            # otherwise it should be from query:
#            elif crm[0].startswith('T.'):
#                results.new_entry('T.Scaff', 0)

    # Populate the results dict with chromosomes names:
    for crm in unique_scaff:
        # If the crm is a scaffold/contig
        # (not a main chromosome):
        # Put them all in a single pool of Scaffolds...

        # crm[0]: scaffold's name
        # crm[1]: scaffold's length

        if 'Scaff' in crm[0]:
            if crm[0].startswith('Q.'): # query scaffold
                # suma la llargada al conjunt d'scaffolds
                results.update_entry('Q.Scaff', ll=int(crm[1]) )
            elif crm[0].startswith('T.'): # target scaffold
                # suma la llargada al conjunt d'scaffolds
                results.update_entry('T.Scaff', ll=int(crm[1]) )
        else:
            results.new_entry( crm[0], crm[1] )

    # Calcula la llargada completa de l'alineament, la suma
    # de totes les bases de cada cromosoma, sigui Q o T.
    ll_total = 0
    for crm in results.r.keys():
        ll_total += results.r[crm]['ll']
    # Afegeix-ho a general:
    results.update_entry('general', ll=ll_total)

    # Parse and extract results.
    # Add them to the above populated results dictionary.
    for line in parse_paf_alignment_db(paf_file)[0]:
        # Analyze the cig with a function (creates dictionary).
        c = cig_analysis(line.cigar)

        # Change qname or tname if they're generic Scaffolds:
        if 'Scaff' in line.qname:
            qname = 'Q.Scaff'
        else:
            qname = line.qname

        if 'Scaff' in line.tname:
            tname = 'T.Scaff'
        else:
            tname = line.tname

        # Update 'general':
        results.update_entry('general',
                cgm=c['M'],
                cgl=c['D'],
                cga=c['I'],
                cgc=c['compressed'],
                col10=line.matches,
                col11=line.n_bases,
                NM=line.mismatch,
                tp=line.typeA,
                nn=line.ambiguous,
                n_ali=1
                )

        # Update 'Target' Scaffold:
        results.update_entry(tname,
                cgm=c['M'],
                cgl=c['I'],
                cga=c['D'],
                cgc=c['compressed'],
                col10=line.matches,
                col11=line.n_bases,
                NM=line.mismatch,
                tp=line.typeA,
                nn=line.ambiguous,
                n_ali=1
                )

        # Update 'Query' Scaffold:
        results.update_entry(qname,
                cgm=c['M'],
                cgl=c['D'],
                cga=c['I'],
                cgc=c['compressed'],
                col10=line.matches,
                col11=line.n_bases,
                NM=line.mismatch,
                tp=line.typeA,
                nn=line.ambiguous,
                n_ali=1
                )


    # Create a CSV to easily export unto a spreadsheet:
    print(paf_file) # print the name of the analyzed *.paf

    # Print the headers of the CSV:
    results.print_header()

    # For all chromosomes:
    for crm in results.r.keys():
    # Print the results:
        results.print_results( crm )

    print() # Espai estètic per separar múltiples resultats.




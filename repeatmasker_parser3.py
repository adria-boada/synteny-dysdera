#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# repeatmasker2.py
#
# 09 may 2022  <adria@molevol-OptiPlex-9020>

""" Parses a repeatmasker output and returns a table
with the amount, length and average length of each
group of classified repetitions.

For example, returns a row for unspecified DNA TEs,
another one for unspecified LTRs, another one for 
TEs classified in the Gypsy superfamily...

To parse only for a given scaffold, create a new
sub-file with grep "Scaffold" < mainfile_reps.out?
"""

""" IMPROVEMENTS/ADDITIONS:
"""

import sys

# Create Markdown-like tables, effortlessly, from lists of lists.
import tabulate

import gzip

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def unique_class_list(fn):
    """Funció per extreure les famílies úniques d'un fitxer de
    repeticions de RepeatMasker. 
    """

    out=[]
    
    with gzip.open(fn) if fn.endswith('.gz') else open(fn) as fh:
        for line in fh:
            # Create a Repeat obj.
            obj = Repeat(line)
        
            # Si repout (la classe extreta) és nova:
            if not obj.repclass in out:
                out += [obj.repclass]

    out.sort()
    return out


# Simple class that parses a line from these files and stores
# the individual values in its fields:
class Repeat(object):
    def __init__(self, ln):
        """ Parse fields given a line.split() object
        as input and store them.
        """
        # Split the line by columns
        ln = ln.split()

        (self.swsc,  # score
                self.pctdiv,
                self.pctdel, 
                self.pctins, 
                self.refid, # query sequence
                self.begin, # begin in query
                self.end, # end in query
                self.ref_remain, 
                self.orient, # strand
                self.rep_nm,
                self.original_cl, # original assigned classification
                self.rep_i, # begin in repeat
                self.rep_f, # end in repeat
                self.rep_prior, 
                self.id, # arbitrary id, from 1 to infinity.
                self.unknown, # unknown type of data.
                self.quality, # quality tag, either '.' or '?'
                self.Class,  # retrotrans., dnatrans., tandem repeat...
                self.Subclass,  # DNA subclass one or two.
                self.Order,  # LTR, DIRS, etc.
                self.Super_family,
                self.notes # MITE or nMITE?
        ) = ln
        # int-ize the reference coordinates
        self.begin, self.end = int(self.begin), int(self.end)
        # Calculate the repetitive element's length:
        self.replen = abs(self.end - self.begin) + 1

        # Get a unique tag for each repeat, excluding
        # the last 'N/A's.
        repclasses = (self.Class, self.Subclass, self.Order, self.Super_family)
        # [0]: retrotrans., dnatrans. or others. 
        # [1]: DNA subclass or tandem repeat subclass.
        # [2]: Transposable Element Order.
        # [3]: Transposable Element Super-family.
        (a, b, c, d) = repclasses

        if d != 'N/A':
            repout = repclasses
        elif c != 'N/A':
            repout = repclasses[:-1] 
        elif b != 'N/A':
            repout = repclasses[:-2] 
        else:
            repout = repclasses[:-3]

        self.repclass = repout


def parse_repeats_with_results(fn, crm='all'):
    """ Parse a RepeatMasker file, and for each line extract
    results. Do not store all the line in a list, do the
    analysis line by line, and afterwards void it.

    In this way we ensure that the computer does not
    freeze all the time.

    Parse the results only from crm. argument...
    only includes line if crm. arg. matches regularly
    with self.refid column.
    The default value parses all scaffolds, not taking
    into account if it is either main or secondary.
    """

    # Needs the preceding `unique_class_list` function.
    unique_list = unique_class_list(fn)
    # Create a results dict.
    results = {}

    """ Initialize the results list, in a JSON format.
    """
    # Init class (first column)
    for i in unique_list:
        if len(i) == 1:
            results[i[0]] = {'amount': 0, 'length': 0}
        elif len(i) > 1:
            results[i[0]] = {}
            results[i[0]]['amount'] = 0
            results[i[0]]['length'] = 0
    # Init subclass (second column)
    for i in unique_list:
        if len(i) == 2:
            results[i[0]][i[1]] = {'amount': 0, 'length': 0}
        elif len(i) > 2:
            results[i[0]][i[1]] = {}
            results[i[0]][i[1]]['amount'] = 0
            results[i[0]][i[1]]['length'] = 0
    # Init Order (third col.)
    for i in unique_list:
        if len(i) == 3:
            results[i[0]][i[1]][i[2]] = {'amount': 0, 'length': 0}
        elif len(i) > 3:
            results[i[0]][i[1]][i[2]] = {}
            results[i[0]][i[1]][i[2]]['amount'] = 0
            results[i[0]][i[1]][i[2]]['length'] = 0
    # Init Super-Family (fourth col.)
    for i in unique_list:
        if len(i) == 4:
            results[i[0]][i[1]][i[2]][i[3]] = {'amount': 0, 'length': 0}

    # DEBUG:
    # print(results)

    # List of filtered in scaffolds.
    scaffolds = []
    
    if crm=='all':
        """ Parse the repeatmasker file. For each line, append
        the amount and length of each classified element.
    
        This block takes in all scaffolds.
        """
        # Can open both gzipped or normal files...?
        with gzip.open(fn) if fn.endswith('.gz') else open(fn) as fh:
            for line in fh:
                # Create a Repeat obj.
                obj = Repeat(line)
                # Get the repclass from them.
                x = obj.repclass
                if len(x) == 1:
                    # DEBUG:
                    # print("1.", x)
                    results[x[0]]['amount'] += 1
                    results[x[0]]['length'] += obj.replen
                elif len(x) == 2:
                    # DEBUG:
                    # print("2.", x)
                    results[x[0]][x[1]]['amount'] += 1
                    results[x[0]][x[1]]['length'] += obj.replen
                elif len(x) == 3:
                    # DEBUG:
                    # print("3.", x)
                    results[x[0]][x[1]][x[2]]['amount'] += 1
                    results[x[0]][x[1]][x[2]]['length'] += obj.replen
                elif len(x) == 4:
                    # DEBUG:
                    # print("4.", x)
                    results[x[0]][x[1]][x[2]][x[3]]['amount'] += 1
                    results[x[0]][x[1]][x[2]][x[3]]['length'] += obj.replen

    else:
        """ Parse the repeatmasker file. For each line, append
        the amount and length of each classified element.

        This block takes in only specified scaffold regexp.
        """
        # Can open both gzipped or normal files...?
        with gzip.open(fn) if fn.endswith('.gz') else open(fn) as fh:
            for line in fh:
                # Create a Repeat obj.
                obj = Repeat(line)
                
                # If the regexp is in the line's REF ID (scaffold name):
                # .casefold() method is recommended in case-insensitive comparison.
                if crm.casefold() in obj.refid.casefold():
                    # Get the unique and allowed REFIDs into a list.
                    if not obj.refid in scaffolds:
                        scaffolds += [ obj.refid ]
                    # Get the repclass from them.
                    x = obj.repclass
                    if len(x) == 1:
                        # DEBUG:
                        # print("1.", x)
                        results[x[0]]['amount'] += 1
                        results[x[0]]['length'] += obj.replen
                    elif len(x) == 2:
                        # DEBUG:
                        # print("2.", x)
                        results[x[0]][x[1]]['amount'] += 1
                        results[x[0]][x[1]]['length'] += obj.replen
                    elif len(x) == 3:
                        # DEBUG:
                        # print("3.", x)
                        results[x[0]][x[1]][x[2]]['amount'] += 1
                        results[x[0]][x[1]][x[2]]['length'] += obj.replen
                    elif len(x) == 4:
                        # DEBUG:
                        # print("4.", x)
                        results[x[0]][x[1]][x[2]][x[3]]['amount'] += 1
                        results[x[0]][x[1]][x[2]][x[3]]['length'] += obj.replen

    # DEBUG:
    # print(results)

    # Return the dictionary with the computed results.
    return (results, scaffolds)


def walkresults(r, i=0):
    """ Funció recursiva que entra dins un diccionari
    amb múltiples nivells (fins a `cols` nivells), 
    busca les branques apicals amb dos resultats,
    `amount` i `length`,
    els imprimeix amb format tabular, i torna enrere
    per imprimir els resultats dels ítems antecedents.

    S'ha de cridar amb for k, v in results: walkdict(v)
    """ 
    out = []

    for k, v in r.items():
        if ((k=="amount") or (k=="length")):
            # skip the amount or length values...
            continue
        # ln = v['length']
        # am = v['amount']
        # Si el dict conté més de cols (dos) ítems:
        # Si no és l'última branca apical, indaga més profundament:
        if len(v) > 2:
            out += walkresults(v, i+1)

        TE_name = ["."]*i + [str(k)] + ["."]*(3-i)
        # Si no es divideix per zero:
        if int(v['amount']) != 0:
            avg_length = int(v['length']) / int(v['amount'])
        else:
            avg_length = 0

        out += [ TE_name + [v['amount'], v['length'], round(avg_length, 1)] ]
    return out
        

if __name__ == '__main__':

    # Instruccions respecte els arguments necessaris per cridar l'script:
    if len(sys.argv) < 2:
        sys.exit('\nCrit script: script.py <repeatmasker.parseable.out> [<Match-this-scaffold-regex>]\n')

    # Load repeatmasker file as a list of "line" classes.
    # repeats = parse_repeat_masker_db(sys.argv[1])
    # YOU DO NOT NEED ANYMORE TO LOAD A FILE (OTHERWISE
    # FREEZES AND CRASHES THE PC).

    # Si hi han més de dos arguments a la terminal:
    if len(sys.argv) > 2:
        # Usa el 2n argv. per buscar REFID, doncs:
        (results, scaffolds) = parse_repeats_with_results( sys.argv[1], sys.argv[2] )
        # Comença a printar header de resultats.
        print()
        print("Analysed File:", sys.argv[1])
        print("Filtered case-insensitively by '", sys.argv[2], "' in the refid column. Allowed scaffolds list:", sep="")
        [print(x) for x in scaffolds]
        print()
    else:
        # Parse the given file in the first argv.
        (results, scaffolds) = parse_repeats_with_results( sys.argv[1] )
        # Comença a printar header de resultats.
        print()
        print("Analysed File:", sys.argv[1])
        print()
    
    # Crea capçalera:
    capçalera = ["Class", "Sub-Class", "Order", "Super-family",
            "Amount", "Base-pair-len", "Avg. Length"]
    # print the actual results?
    # tabulated thanks to tabulate() module (tabulates a given list).
    data_listed = walkresults(results)

    # Add, at the end of the list, a sum of each column:
    # First, loop through the result list and compute it:
    total = {'amount':0, 'BP-Total':0, 'avg':0, 'avg.counter':0}
    for i in data_listed:
        # Suma la columna quatre (nombre)
        total['amount'] += int(i[4])
        # Suma la columna cinc (bp-len)
        total['BP-Total'] += int(i[5])
        # Suma la columna sis (average)
        if int(i[6]) != 0:
            total['avg'] += int(i[6])
            # Tingues un comptador a mà:
            # o utilitza len(data_listed) com a comptador...
            total['avg.counter'] += 1


    # Separate it from the rest of the table:
    data_listed += [ ('----------', '----------', '----------', '----------', '----------', '----------', '----------') ]
    data_listed += [ ('', '', '', '', '', '', '') ]

    # DEBUG
    if total['avg.counter']!=0:
        data_listed += [ ('Total', '.', '.', '.', total['amount'], total['BP-Total'], round(total['avg']/total['avg.counter'], 1) ) ]
    else:
        data_listed += [ ('Total', '.', '.', '.', total['amount'], total['BP-Total'], 0) ]

    print( tabulate.tabulate (data_listed, headers=capçalera, tablefmt='pipe') )




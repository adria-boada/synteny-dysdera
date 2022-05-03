#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# paf_parser.py
#
# 2022-02-21  <adria@molevol-OptiPlex-9020>

""" parse paf-file to measure coordinates covered for query and target 
sequences.

Help on solving this problem: 

https://www.baeldung.com/cs/finding-all-overlapping-intervals

"""

import sys

# Funcions per evaluar solapament d'intervals:
import overlap_paf_parse as papaf


# Instruccions respecte els arguments necessaris per cridar l'script:
if len(sys.argv) < 2:
    sys.exit('\nCrit script: script.py <paf-file.paf>\n')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Define paf-file.
paf_file = sys.argv[1]

# Creates a dict in which to store intervals
# of alignments for each chr. 
seqdic = {}

# Create a dict in which to store definitive results.
results = {}

# Later on both preceding dicts will get
# populated with alignment sequences/chromosomes.

# Populate dictionaries with chromosomes.
with open(paf_file) as paf:
    for line in paf:
        # Split the lines by tabs and cut columns
        # we're interested in:

        # Query and Target name+len for this line
        # (line = one alignment/hit)
        qname = "Q." + line.split("\t")[0]
        qlen = int( line.split("\t")[1] )

        tname = "T." + line.split("\t")[5]
        tlen = int( line.split("\t")[6] )

        # Create an entry for each chromosome.
        # Col·lapsing length of intervals:
        seqdic[qname] = []
        seqdic[tname] = []
        # Length of supplied Q and T.

        results[qname] = {
                'll': qlen,
                'overlap': 0,
                'cgMatch': 0,
                'cgIndel': 0,
                'cgOpp': 0,
                'cgCompr': 0,
                'col10': 0,
                'col11': 0,
                'non-Ovlap': 0,
                'NM': 0
                }
        results[tname] = {
                'll': tlen,
                'overlap': 0,
                'cgMatch': 0,
                'cgIndel': 0,
                'cgOpp': 0,
                'cgCompr': 0,
                'col10': 0,
                'col11': 0,
                'non-Ovlap': 0,
                'NM': 0
                }
        # Other info (?):



# Prepare a couple of vars
# to store general stats...

# Total sum of added to Target
# (not present in Query):
# sum(cgIndel) from only Target
total_sum_added_to_query = 0
# Total sum of lost from Target
# (only present in Target):
# sum(cgOpp) from only Target
total_sum_lost_from_query = 0

# Append the data once the dictionaries
# have been populated:
with open(paf_file) as paf:
    for line in paf:

        # Query and Target namefor this line
        # (line = one alignment/hit)
        qname = "Q." + line.split("\t")[0]
        tname = "T." + line.split("\t")[5]

        # Extract coordinate range (start, end)  covered 
        # by the alignment. Also, bases mapped.
        # Don't forget about being gap inclusive...
        # 'chr' can either be query or target.
        # seqdic[chr][0 to 1]: start to end.
        # seqdic[chr][2]: bases mapped excluding indels.
        # seqdic[chr][3]: b.    m.     including indels.
        seqdic[qname] += [ [ int(x) for x in line.split("\t")[2:4] ] + [int( line.split("\t")[9])] + [int( line.split("\t")[10])] ]

        # Do the same for target names:
        seqdic[tname] += [ [ int(x) for x in line.split("\t")[7:9] ] + [int( line.split("\t")[9])] + [int( line.split("\t")[10])]]


        # cigar string (for both Q and T):
        # [-1]: cuts the last column.
        # [5:-1]: cuts beggining 'cg:Z:' and ending '\n'.
        cig = line.split("\t")[-1][5:-1]

        # papaf.cig_analysis(cig) # how to call the f(x).
        # returns a dict with the format {'M':m,'I':i,'D':d,'compressed':g}
        c = papaf.cig_analysis(cig)

        # Sum CIGAR matches...
        results[qname]['cgMatch'] += c['M']
        results[tname]['cgMatch'] += c['M']

        # Sum CIGAR indels and opposite gaps.
        # 'I' means gap only present in query.
        # 'D' means gap only present in target.
        results[qname]['cgIndel'] += c['I']
        results[tname]['cgIndel'] += c['D']

        results[qname]['cgOpp'] += c['D']
        results[tname]['cgOpp'] += c['I']

        total_sum_lost_from_query += c['D']
        total_sum_added_to_query += c['I']

        results[qname]['cgCompr'] += c['compressed']
        results[tname]['cgCompr'] += c['compressed']
       
        # matches and gaps:
        # col10: matching bases in mapping:
        col_ten = int( line.split("\t")[9] )

        results[qname]['col10'] += col_ten
        results[tname]['col10'] += col_ten

        # col11: matching bases including gaps (indels).
        col_eleven = int( line.split("\t")[10] )
        
        results[qname]['col11'] += col_eleven
        results[tname]['col11'] += col_eleven

        # Number of MisMatches (NM, col13)
        # Sum of missmatches and gaps in ali.
        nmiss = int( line.split("\t")[12][5:] )

        results[qname]['NM'] += nmiss
        results[tname]['NM'] += nmiss



# Overlapping algorithm:
# SWEEP LINE APPROACH (ref.: blog-post in description).

# Create an ordered list with all the points 
# (both Left and Right) for every interval range.

# look at module overlaping_paf_parser.py.

# Sort the Q-T-dictionary:
for crm, l in seqdic.items():
    # Sort each list of intervals:
    seqdic[crm].sort()

for crm, l in seqdic.items():

    # Create couples of overlapping intervals.
    solapats = papaf.couples_overlapping(seqdic[crm])

    # Create thorough overlapping interval indices.
    thorough_idx = papaf.thorough_indices(seqdic[crm])

    # Create thorough idx's complimentary list
    # (all non-overlapped intervals)
    # Afterward, measure mapped bases from this list.
    integre = {'tot':0, 'gap':0}

    # Per a cada index de la llista completa:
    for i in range(0, len(seqdic[crm])):
        # Suma'n les bases si no és solapat:
        if i not in thorough_idx:
            integre['tot'] += l[i][2]
            integre['gap'] += l[i][3]

    results[crm]['non-Ovlap'] = integre

    # Start overlapping-interval-fusion:
    while solapats:

        extensos = []

        for i in solapats:
            # Extreu una parella solapada de la llista (i[0] i i[1])
            overlap = [ seqdic[crm][i[0]] , seqdic[crm][i[1]] ]

            # Crea nous intervals extensos (suma intervals parella):
            # Agafa el mínim de l'esquerra (x[0]) i el màxim de la dreta (x[1]) 
            extensos += [[ min( [x[0] for x in overlap] ), max( [x[1] for x in overlap] ) ]]

        # Afageix un complementari de la llista d'indices solapats
        # (és a dir, afageix objectes que no solapen a 'extensos'):
		
        # Necessita una llista on pugui avaluar si index pertany
        # (les subllistes de 'solapats' molesten en aquest pas).
        des_solapats = [x[0] for x in solapats] + [x[1] for x in solapats]
        des_solapats.sort()

        # DEBUGGGGGGGGG:
        #print(des_solapats, extensos)
        #print("LEN:", len(seqdic[crm]) )
		
        # per a cada possible índex:
        for i in range(0, len(seqdic[crm])):
			
            # Guarda'l si no es troba dins dels solapats:
            if i not in des_solapats:
                extensos += [ seqdic[crm][i] ]
            
        # Ordena els intervals segons la 1era coordenada.
        extensos.sort()
		
        # Guarda els nous intervals extensos al dict.
        seqdic[crm] = extensos
		
        # DEBUG
        #print(extensos, "\n"*5, seqdic[crm])
		
        # Crea una altra llista de punts per tornar a iterar per
        # la funció que extreu parelles de solapats:
        point_list = []
	
        # Create a list with all the points (both left and right)
        # of every interval.
        for i in range(0, len(extensos)):
            # Left point == 0.
            point_list += [ [extensos[i][0], 0, i] ]
            # Right point == 1.
            point_list += [ [extensos[i][1], 1, i] ]
		
        # Ordena la llista de punts per trobar solapaments:
        point_list.sort()
		
        # Buida solapats i torna a iterar per veure si troba
        # solapaments addicionals:
        solapats = papaf.couples_overlapping (seqdic[crm]) 

        # Si en troba torna a dalt del while;
        # Alternativament, continua.

# Final algorisme que extreu intervals solapats.

# Computa distància (start, end) i suma-la
# a 'overlap'.
for crm, l in seqdic.items():
	
    # Per a cada interval [inici, final]:
    for interval in l:
        results[crm]['overlap'] += int(interval[1] - interval[0] + 1)



# Let us compute general over-arching stats:

# Length of the sum of all chromosomes:
total_sum_alibases = 0
# Length of the sum of mapped bases (col11):
total_sum_map_w_indels = 0
# Length of the sum of 'correctly' mapped bases (col10):
total_sum_mapped = 0
# Length of indels + mismatches:
total_sum_mindles = 0
# Length of only mismatches.
total_sum_misses = 0
# Length of matches and misses:
total_sum_overall_hit = 0
# Compressed gaps sum:
total_sum_compressed = 0

for crm, dict_val in results.items():
    total_sum_alibases += dict_val['ll']
    total_sum_map_w_indels += dict_val['col11']
    total_sum_mapped += dict_val['col10']
    # mindels is col11-col10.
    total_sum_mindles += dict_val['col11'] - dict_val['col10']
    # sense encertar = NM - gaps = errors + gaps - gaps
    total_sum_misses += dict_val['NM'] - dict_val['cgOpp'] - dict_val['cgIndel']    
    # sum of cigar M (overall hit lengths)
    total_sum_overall_hit += dict_val['cgMatch']
    # Sum the amount of compressed gaps:
    total_sum_compressed += dict_val['cgCompr']

# BLAST-identity-like parameter:
blast_id = total_sum_mapped / total_sum_map_w_indels

# Average coverage of genome by alignment:
avg_cov = total_sum_map_w_indels / total_sum_alibases

# ~ # # ~ # # ~ # # ~ #  # ~ #  # ~ #  # ~ #  # ~ #  # ~ #
# PRINTING() RESULTS...

# Transforma un nombre ( x*10**6 ) bases a x Megabasepairs.
# fins a 'decims' decimals mostrats.
def round_to_Mbps(number:"Amount of bases", decims:"Nº of trailing decimals"=2):
    return round( number/(10**6), decims)


#print() # Espai estètic
#
#print("Suma de la ll. de tots els CRM.s: ", round_to_Mbps(total_sum_alibases), "Mbps" )
#print("Average coverage of genomes:      ", round(avg_cov*100, 3), "%" )
#
## The following stats are divided by 2 because I summed Q+T (oops).
#print("Total mapped bases (w/ indels):   ", round_to_Mbps(total_sum_map_w_indels/2), "Mbps" )
#print("Mapped bases which match:         ", round_to_Mbps(total_sum_mapped/2), "Mbps" )
#print("Mapped bases which mismatch:      ", round_to_Mbps(total_sum_misses/2), "Mbps" )
## La quantitat d'indels/gaps és la suma de gaps i misses menys els misses...:
#print("Mapped bases which are gaps:      ", round_to_Mbps((total_sum_mindles-total_sum_misses)/2), "Mbps" )
#
#print("Ratio indels + misses / total:    ", round((total_sum_mindles/total_sum_map_w_indels)*100, 3), "%" )
#print("Ratio Matches/total == BLAST-id:  ", round(blast_id*100, 3), "%" )
#print("Gap-compressed identity:          ", round((total_sum_mapped*100)/(total_sum_overall_hit+total_sum_compressed), 3), "%" )
#
#print() # Espai estètic.
#
## Specific chromosome stats:
#for crm in results.keys():
#    # Shorter variable calling:
#    r = results[crm]
#
#    print("*"*70, f"\nFor chromosome {crm}:" )
#
#    print("Llargada cromosoma:               ", round_to_Mbps(r['ll']), 'Mbps')
#    
#    # col11 length: adding INDELS to col10 (sum of cigar's M+D+I).
#    # Doesn't discriminate insertions from deletions and viceversa. 
#    print("Bases hit/Mapped bases w/ indels: ", round_to_Mbps(r['col11']), 'Mbps')
#    print("Ratio Match/Base hit == BLAST-id: ", round( r['col10']/r['col11'], 3) )
#
#    # col10 length: matches and mismatches
#    # llargada de regions homòlogues.
#    print("Mapped bases which match:         ", round_to_Mbps(r['col10']), 'Mbps')
#    print("Mapped bases which mismatch:      ", round_to_Mbps(r['NM']-r['cgOpp']-r['cgIndel']), 'Mbps')
#    print("Percent matching bases:           ", round( r['col10']/r['ll'], 3) )
#
#    print("CIGAR hit length (match+miss):    ", round_to_Mbps(r['cgMatch']), 'Mbps')
#    print("CIGAR Added bases to crm:         ", round_to_Mbps(r['cgIndel']), 'Mbps')
#    print("CIGAR Lost bases from crm:        ", round_to_Mbps(r['cgOpp']), 'Mbps')
#    # 'NM' includes gaps; rest them to find mismatches.  
#
#    # solapats hauria de ser M+I o M+D, no?
#    print("Ll. alineament col·lapsat:          ", round_to_Mbps(r['overlap']), 'Mbps')
#    print("Percentatge col·lapsat:             ", round( r['overlap']/r['ll'], 3) )
#
#    # Gràcies al blog de lh3 
#    # https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
#    # Gap-compressed identity:
#    print("Gap excluded id.:                 ", round(r['col10']/r['cgMatch'], 3))
#    print("Gap compressed id.:               ", round(r['col10']/(r['cgMatch']+r['cgCompr']), 3))
#
#    print() # Espai estètic.


# Crea un final estil CSV per exportar fàcilment DATA a un full de càlcul:
print("^CSV\n - General - ")

print("File", 
        "# Sum(every CHR len)", 
        "# Bases with hit", 
        " ", # blank formatting line (follows colnames from 'specific')
        " ", # blank formatting line (follows colnames from 'specific')
        "% CRM all gaps", 
        "# matches", 
        "# misses", 
        "# «Added» to Qry.", #added to query == query insertions
        "# «Lost» from Qry.", #lost from query == query deletions
        "# sum indels", 
        " ", # blank formatting line
        "% ~BLAST-id", 
        "% Gap-compressed-id", 
        sep=";")

print( sys.argv[1], #file
        round_to_Mbps( total_sum_alibases), #Sum(every CHR len)
        round_to_Mbps( total_sum_map_w_indels/2 ), #Bases with hit.
        # Divided by two because I summed bases from both target and query
        # 30 aligned from query + 30 aligned from target => 30 aligned bases overall.
        " ", # blank formatting line
        " ", # blank formatting line
        f"{round( avg_cov*100, 4)} %", # %(Genomes len.) covered by alignment.
        round_to_Mbps( total_sum_mapped/2), #sum(Matches)
        round_to_Mbps( total_sum_misses/2), #sum(Misses)
        round_to_Mbps( total_sum_added_to_query), #sum(Added)
        round_to_Mbps( total_sum_lost_from_query), #sum(Lost)
        round_to_Mbps( (total_sum_mindles - total_sum_misses)/2), #sum(gaps)
        " ", #blank line for formatting (stay in line with 'specific')
        f"{round(blast_id*100, 4)} %", #blast-id
        f"{round((total_sum_mapped/(total_sum_overall_hit+total_sum_compressed))*100, 4)} %", #gap-compressed
        sep=";")

# Print the stats for each chromosome
print() # Espai estètic
print(" - Especific per cada CRM - ")
print("DNA sequence",  #seqname
        "# CHR len",  #chrlen
        "# Bases with hit",  #col11
        "% CRM gapless", #(match+ miss) /chrlen
        "% CRM partial gaps",  #(match+ miss+ addedgaps) /chrlen
        "% CRM all gaps", #col11 / chrlen
        "# matches",  #matches
        "# misses",  #misses
        "# «Added»",  #added seq in crm
        "# «Lost»",  #lost seq in crm
        "# sum indels",  #sum of added and lost 
        "# bases ali col·lapsat",  #colapsed compute val
        "% ~BLAST-id",  #matches/(col11)
        "% Gap-compressed-id",  #gap-compressed==each gap one mismatch
        sep=";")

for crm in results.keys():

    # Shortens variable calls.
    r = results[crm]

    # Compute variables.
    gapless_cov = round( r['col10'] / r['ll']*100 , 4 ) #matches
    partial_cov = round( ((r['cgMatch']+r['cgIndel']) / r['ll'])*100, 4 ) #matches with misses and "added" gaps
    complete_cov = round( (r['col11'] / r['ll'])*100, 4 ) #both gaps, added and lost.
    blast_id = round( (r['col10'] / r['col11'])*100, 4 ) #identity
    gap_comp = round( (r['col10'] / (r['cgMatch']+r['cgCompr']))*100, 4 ) #identity

    print(crm, # crm name
            round_to_Mbps(r['ll']), #llargada chr
            round_to_Mbps(r['col11']), #bases amb hit
            f"{gapless_cov} %", #discounting all gaps
            f"{partial_cov} %", #coverage discounting "lost" gaps
            f"{complete_cov} %", #coverage counting all gaps
            round_to_Mbps(r['col10']), #only matches.
            round_to_Mbps(r['NM']-r['cgOpp']-r['cgIndel']), #only misses.
            round_to_Mbps(r['cgIndel']), #added
            round_to_Mbps(r['cgOpp']), #lost
            round_to_Mbps( r['cgOpp'] + r['cgIndel'] ), #sum of gaps
            round_to_Mbps( r['overlap'] ), #sum of overlapping bases
            f"{blast_id} %", 
            f"{gap_comp} %", 
            sep=";")




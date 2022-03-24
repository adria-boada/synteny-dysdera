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

# Create a dict. in which to store chr coverage.
seqdic = {} 

# Create a dict in which to store chr total len
chr_total_len = {}

# Create a list in which to store gaps (to calc. gap percent)
gaps_len = []

# Get all chr-names (query and target):
with open(paf_file) as paf:
	for line in paf:
		# Split the lines by tabs, and cut column one [0].
		qname=line.split("\t")[0]
		# Create an entry in the dict for each query name
		seqdic[f"Q.{qname}"]=[]
		
		# Split the line by tabs and cut qlen column:
		qlen = line.split("\t")[1]
		chr_total_len[f"Q.{qname}"] = qlen
		
		# Do the same for target names:
		tname=line.split("\t")[5]
		seqdic[f"T.{tname}"]=[]
		
		tlen = line.split("\t")[6]
		chr_total_len[f"T.{tname}"] = tlen

		# cut out the gaps column (col 10 and 11):
		# col 10: matching bases in mapping;
		# col 11: including gaps, mismatches... etc.
		temp_gaps = [ int( line.split("\t")[9] ), int( line.split("\t")[10] ) ]
		# amb els nombres de la línia anterior,
		# fés un càlcul de bases ocupades per gaps:
		temp_gaps += [ (temp_gaps[1] - temp_gaps[0]) ]
		
		# temp_gaps = [ col10, col11, col11 - col10 ]
		
		# Store the computed vals in a sub-list, for posterity;
		gaps_len += [ temp_gaps ]
		

# Compute the comparison of gaps versus mapped bases:	
total_gap_len = 0
total_bases_mapped = 0
total_missand_matches_mapped = 0

for x in gaps_len:
	total_gap_len += x[2]
	total_bases_mapped += x[1]
	total_missand_matches_mapped += x[0]
	
# Sum of total length of alignment:
total_sum_bases_ali = 0
for chr, x in chr_total_len.items():
	total_sum_bases_ali += int(x)

# DEBUG:
# Rounds most values to 2 places and displays Megabasepairs.
print("Sum of alignment bases (sum of all chr len):", round(total_sum_bases_ali/(10**6), 2), "Mbps" )
print("Total mapped bases:", round(total_bases_mapped/(10**6), 2), "Mbps" )
print("Total 'correctly' mapped bases:", round(total_missand_matches_mapped/(10**6), 2), "Mbps" )
# ~ print("SUM(Mapped bases)/SUM(whole_aln_len)/ 2:", round(total_bases_mapped/(total_sum_bases_ali/2), 4) )
print("Total gaps:", round(total_gap_len/(10**6), 2), "Mbps" )
# ~ print("Percent gaps:", round(total_gap_len/total_sum_bases_ali, 4) )
print("Ratio Gaps/Mapped:", round(total_gap_len/total_bases_mapped, 4) )

# View the created dictionaries:
# ~ print(seqdic)


# Get all coordinate ranges covered by the alignment:
with open(paf_file) as paf:
	for line in paf:
		# Split the lines by tabs, and cut column one [0].
		qname=line.split("\t")[0]
		# Extract coordinate range (start, end) and bases mapped.
		# Don't forget about being gap inclusive...
		# seqdic[chr][0 to 1]: start to end.
		# seqdic[chr][2]: bases mapped excluding gaps 
		# seqdic[chr][3]: b.    m.     including gaps. 
		seqdic[f"Q.{qname}"] += [ [ int(x) for x in line.split("\t")[2:4] ] + [int( line.split("\t")[9])] + [int( line.split("\t")[10])] ]
		
		# Do the same for target names:
		tname=line.split("\t")[5]
		seqdic[f"T.{tname}"] += [ [ int(x) for x in line.split("\t")[7:9] ] + [int( line.split("\t")[9])] + [int( line.split("\t")[10])]]
		
		
		
for kaho, x in seqdic.items(): 
	
	total_map_perchr = 0
	for i in x:
		total_map_perchr += i[2]
		
	# DEBUG:
	# Rounds most values to 2 places and displays Megabasepairs.
	print(f"Mapatge de {kaho}:", round(total_map_perchr/(10**6), 2), "Mbps" )
	print(f"Llargada de {kaho}:", round(int(chr_total_len[kaho])/(10**6), 2), "Mbps" )
	print(f"Percentatge de {kaho}:", round((total_map_perchr/int(chr_total_len[kaho])), 3) )
		
		

# Sort the Q-T-dictionary: 
for key, val in seqdic.items():
	# Sort each list of values.
	val.sort()



# View the created dictionaries:		
# ~ print(seqdic)


# Create a dictionary to store lengths for every chr.
lendic = {}

# For each { key: [list of intervals] }, store { key: nº of intervals }
# in the lendic dictionary.
for key, x in seqdic.items(): lendic[key] = [len(x)]

# DEBUG
# ~ print(lendic)


 ##### NAÏVE ALGORITHM
 
# ~ # Per a cada cromosoma de la llista:
# ~ for chr, length in lendic.items():
	
	# ~ # Create a list to get the amount of overlapping intervals
	# ~ amount = []
	
	# ~ # Faci l'algorisme que troba intervals solapats.
	# ~ for i in range(0, length-1):
		# ~ for j in range(0, length-1):
			# ~ # No miris interval i==j; com que son el mateix
			# ~ # segur que solapen un amb l'altre.
			# ~ if i==j: continue
			# ~ if (max(seqdic[chr][i][0], seqdic[chr][j][0] )) <= (min(seqdic[chr][i][1], seqdic[chr][j][1] )):
				# ~ print("{}:".format(chr), Qseqdic[chr][i])
				# ~ amount += Qseqdic[chr][i]
				# ~ break

# ~ print("Total amount of overlapping genes is", len(amount) )



 ##### SWEEP LINE APPROACH


def point_list_overlapping_indices (p_list: "A list of sorted points from *point_list*"):
	
	""" Acquire and return indices of overlapping intervals from 
	*point_list*. Returns couples of overlapping intervals in the form 
	of sublists. Only does one sweep to find pairs; to find multiple 
	overlaps an iterative approach is required. """
	
	# State variables for algorithm:
	currentOpen = -1
	answer = []
	
	# for each point in the point_list:
	for i in range(0, len(p_list)):
		
		# Si el punt actual és 'Left; opens an interval'
		if p_list[i][1] == 0:
			
			# Si no hi ha cap interval actualment obert:
			if currentOpen == -1:
				# 'Entra' a l'interval 'i'.
				currentOpen = p_list[i][2]
			
			# Si es trobava dins un interval anterior:
			else:
				# Index del nou interval:
				index = p_list[i][2]
				# Afageix la parella a la resposta:
				answer += [ [currentOpen, index] ]
				# Tanca l'interval previ
				currentOpen = -1

		# Si troba un punt que és 'Right'; que tanca interval:
		else:
			# Surt del currentOpen.
			currentOpen = -1
	
	return answer




# Per a cada chr, sigui query o target, dins el paf-file:
for chr, length in lendic.items():

	point_list = []

# Create a list with all the points (both left and right)
# of every interval.
# lendic is now a list so must be called with len[0]. 
	for i in range(0, length[0]):
		# Left point == 0.
		point_list += [ [seqdic[chr][i][0], 0, i] ]
		# Right point == 1.
		point_list += [ [seqdic[chr][i][1], 1, i] ]
	
	# Sort the list by position (first val in sublists):
	point_list.sort()
	
	# ~ print("Imprimeix la llista de punts ordenada:")
	# ~ print(f"{chr} : {point_list}") # Sembla que si que ordena correctament.

	# Create couples of overlapping intervals indices. 
	solapats = point_list_overlapping_indices (point_list)
	
	# Create thorough overlapping intervals indices
	thorough_idx = papaf.thorough_indices(seqdic[chr])
	
	# Create the latter list's complimentary list (non-overlapped):
	# Afterwards, measure the bases mapped of this list (var. 'integre')
	integre = 0
	integre_gaps = 0
	
	# ~ # Necessita(va) una llista on pugui avaluar si index pertany
	# ~ # (les subllistes de 'solapats' molesten en aquest pas).
	# ~ des_solapats = [x[0] for x in solapats] + [x[1] for x in solapats]
	# ~ des_solapats.sort()
	
	# per a cada possible índex de la llista completa:
	for i in range(0, len(seqdic[chr])):
		
		# Suma'n les bases si no és solapat:
		if i not in thorough_idx:
			integre += seqdic[chr][i][2]
			integre_gaps += seqdic[chr][i][3]
	
	# DEBUG:
	# ~ print(f"Non-overlapping mapped b. in chr {chr}:", round(integre/(10**6), 2), "Mbps" )
	# ~ print(f"Non-overlapping gapped b. in chr {chr}:", round(integre_gaps/(10**6), 2), "Mbps" )
	
	# Guarda 'íntegre' a dins de lendic[chr] [1]:
	lendic[chr] += [integre]		

	# Mentre hi hagi objectes dins de solapats:
	while solapats:
		
		# DEBUG:
		# ~ print(solapats)
		
		# Prepara variable per algorisme:
		extensos = []
		
		# Traspassa la fusió de parelles dins solapats a la llista 'extensos'. 
		for i in solapats:
			
			# Extreu una parella de la llista (i[0] i i[1])
			overlap = [ seqdic[chr][i[0]] , seqdic[chr][i[1]] ]
			
			# Crea nous intervals extensos (suma intervals parella):
			# Agafa el mínim de l'esquerra (x[0]) i el màxim de la dreta (x[1]) 
			extensos += [[ min( [x[0] for x in overlap] ), max( [x[1] for x in overlap] ) ]]
			
			# DEBUG:
			# ~ print(overlap, extensos)
	
		# Afageix un complementari de la llista d'indices solapats
		# (és a dir, afageix objectes que no solapen a 'extensos'):
		
		# Necessita una llista on pugui avaluar si index pertany
		# (les subllistes de 'solapats' molesten en aquest pas).
		des_solapats = [x[0] for x in solapats] + [x[1] for x in solapats]
		des_solapats.sort()
		
		# per a cada possible índex:
		for i in range(0, len(seqdic[chr])):
			
			# Guarda'l si no es troba dins dels solapats:
			if i not in des_solapats:
				extensos += [ seqdic[chr][i] ]
		
		# DEBUG
		# ~ print(extensos)
		# ~ print("\n"*20)
		
		# Ordena els intervals segons la 1era coordenada.
		extensos.sort()
		
		# DEBUG
		# ~ print(len(extensos), len(seqdic[chr]), len(solapats))
		
		# Guarda els nous intervals extensos al dict.
		seqdic[chr] = extensos
		
		# DEBUG
		# ~ print(extensos, "\n"*5, seqdic[chr])
		
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
		
		# DEBUG
		# ~ print(point_list)
		
		# Ordena la llista de punts per trobar solapaments:
		point_list.sort()
		
		# DEBUG
		# ~ print(f"{chr}\n", "\n"*4, point_list)
		
		# Buida solapats i torna a iterar per veure si troba
		# solapaments addicionals:
		solapats = point_list_overlapping_indices (point_list)
		
		# Si en troba torna a dalt del while;
		# Alternativament, continua.
		
		
# DEBUG:
# ~ for key, x in seqdic.items():
	# ~ print("**"*38, f"\n{key}", f"\n\n{x}\n")	

# Final algorisme que extreu intervals solapats.

# Calcula la longitud total d'alineaments de cada chr:
llargada_intervals_sumats = {}

for chr, x in seqdic.items():
	llargada_intervals_sumats [chr] = 0
	
	# Per a cada interval [inici, final]:
	for interval in x:
		llargada_intervals_sumats [chr] += int(interval[1] - interval[0] + 1)


# DEBUG:
for chr, x in llargada_intervals_sumats.items():
	print ("*"*78, f"\nFor chromosome {chr}:")
	print ("Total len:", round(int(chr_total_len[chr])/(10**6), 2), "Mbps")
	print ("Aligned len:", round(int(llargada_intervals_sumats[chr])/(10**6), 2), "Mbps")
	print ("Aligned percent:", round(int(llargada_intervals_sumats[chr])/int(chr_total_len[chr]), 4) )
	
	print() # Espai estètic


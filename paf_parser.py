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




# View the created dictionaries:		
print(seqdic)




# Get all coordinate ranges covered by the alignment:
with open(paf_file) as paf:
	for line in paf:
		# Split the lines by tabs, and cut column one [0].
		qname=line.split("\t")[0]
		# Extract coordinate range (start, end).
		seqdic[f"Q.{qname}"] += [ [ int(x) for x in line.split("\t")[2:4] ] ]
		
		# Do the same for target names:
		tname=line.split("\t")[5]
		seqdic[f"T.{tname}"] += [ [ int(x) for x in line.split("\t")[7:9] ] ]
		
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
for key, x in seqdic.items(): lendic[key] = len(x)
print(lendic)

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


# Per a cada chr, sigui query o target, dins el paf-file:
for chr, length in lendic.items():
 
	point_list = []

# Create a list with all the points (both left and right)
# of every interval.
	for i in range(0, length):
		# Left point == 0.
		point_list += [ [seqdic[chr][i][0], 0, i] ]
		# Right point == 1.
		point_list += [ [seqdic[chr][i][1], 1, i] ]
	
	# Sort the list:
	point_list.sort()
	
	print("Imprimeix la llista de punts ordenada:")
	print(f"{chr} : {point_list}") # Sembla que si que ordena correctament.

	# State variables for algorithm:
	currentOpen = -1
	added = False
	answer = []
	
	# for each point in the point_list:
	for i in range(0, len(point_list)):
		
		# Si el punt actual és 'Left; opens an interval'
		if point_list[i][1] == 0:
			# Si no hi ha cap interval actualment obert:
			if currentOpen == -1:
				# 'Entra' a l'interval 'i'.
				currentOpen = point_list[i][2]
				added = False
			# Si es trobava dins un interval anterior:
			else:
				# Index del nou interval:
				index = point_list[i][2]
				# Aquest forma part de la resposta (intervals solapats).
				answer += [index]
				
				# DEBUG: prints each found answer.
				print( [ [chr, index] ] )
				
				# Si l'interval currentOpen no ha sigut encara afegit:
				if (not added):
					# Fes-ho:
					answer += [currentOpen]
					
					# DEBUG: prints each found answer.
					print( [ [chr, currentOpen] ] )
					
					added = True
				# Mira quin dels dos intervals que s'estan comparant
				# te la cua més llarga, i per tant solaparà amb
				# més intervals (cua = Right point).
				if (seqdic[chr][currentOpen][1] < seqdic[chr][index][1]):
					currentOpen = index
					added = True
					
		# Si troba un punt que és 'Right'; que tanca interval:
		else:
			# I és el punt esquerra de l'interval actual:
			if point_list[i][2] == currentOpen:
				# Surt de l'interval
				currentOpen = -1
				added = False

	# Traspassa la answer a intervals dins del cromosoma chr:
	for i in range(0, length):
		# Si 'i' es troba a 'answer', afageix a solapats:
		if i in answer:
			solapats[chr] += [[seqdic[chr][i] ]]
		# else, a no solapats:
		else:
			no_solapats[chr] += [[ seqdic[chr][i] ]]

# LLargada dels solapats:
llargades_solapament = {}
for chr, x in solapats.items(): llargades_solapament[f"S.{chr}"] = len(x)
for chr, x in no_solapats.items(): llargades_solapament[f"NS.{chr}"] = len(x)

print(f"SOLAPATS: ")
print(solapats)

print(f"NO SOLAPATS: ")
print(no_solapats)
				
print(llargades_solapament)



# Anem a mirar coverage per sols els intevals NO solapats:

for chr, value in no_solapats.items():
	# var to count length
	total_len = 0	
	
	# iterate through list and get len of intervals
	for i in value:
		total_len += i[0][1] - i[0][0] + 1
	
	# calculate coverage
	cov = (total_len/ int(chr_total_len[chr]) )
	print(f"For chr {chr}, coverage is {cov}.")	


# for each pair of intervals in the list, remove overlapping:


def point_list_overlapping_indices (p_list: "A list of sorted points from *point_list*"):
	
	""" Acquire and return indices of overlapping intervals from 
	*point_list*. Returns pairs of overlapping intervals in the form of 
	sublists. """
	
	# State variables for algorithm:
	currentOpen = -1
	added = False
	answer = []
	
	# for each point in the point_list:
	for i in range(0, len(p_list)):
		
		# Si el punt actual és 'Left; opens an interval'
		if p_list[i][1] == 0:
			# Si no hi ha cap interval actualment obert:
			if currentOpen == -1:
				# 'Entra' a l'interval 'i'.
				currentOpen = p_list[i][2]
				added = False
			# Si es trobava dins un interval anterior:
			else:
				# Index del nou interval:
				index = p_list[i][2]
				# Aquest forma part de la resposta (intervals solapats).
				answer += [index]
				
				# DEBUG: prints each found answer.
				# ~ print( [ [chr, index] ] )
				
				# Si l'interval currentOpen no ha sigut encara afegit:
				if (not added):
					# Fes-ho:
					answer += [currentOpen]
					
					# DEBUG: prints each found answer.
					# ~ print( [ [chr, currentOpen] ] )
					
					added = True
				# Mira quin dels dos intervals que s'estan comparant
				# te la cua més llarga, i per tant solaparà amb
				# més intervals (cua = Right point).
				if (seqdic[chr][currentOpen][1] < seqdic[chr][index][1]):
					currentOpen = index
					added = True
					
		# Si troba un punt que és 'Right'; que tanca interval:
		else:
			# I és el punt esquerra de l'interval actual:
			if p_list[i][2] == currentOpen:
				# Surt de l'interval
				currentOpen = -1
				added = False

	# Traspassa la answer a intervals dins del cromosoma chr:
	for i in range(0, length):
		# Si 'i' es troba a 'answer', afageix a solapats:
		if i in answer:
			solapats[chr] += [[seqdic[chr][i] ]]
		# else, a no solapats:
		else:
			no_solapats[chr] += [[ seqdic[chr][i] ]]


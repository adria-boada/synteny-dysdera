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
Qseqdic = {} ; Tseqdic = {} 

# Get all query-names:
with open(paf_file) as paf:
	for line in paf:
		# Split the lines by tabs, and cut column one [0].
		qname=line.split("\t")[0]
		# Create an entry in the dict for each query name
		Qseqdic[qname]=[]
		
		# Do the same for target names:
		tname=line.split("\t")[5]
		Tseqdic[tname]=[]




# View the created dictionaries:		
# ~ print(Qseqdic)
# ~ print("\n"+"="*80+"\n")
# ~ print(Tseqdic)




# Get all coordinate ranges covered by the alignment:
with open(paf_file) as paf:
	for line in paf:
		# Split the lines by tabs, and cut column one [0].
		qname=line.split("\t")[0]
		# Extract coordinate range (start, end).
		Qseqdic[qname]+= [ [ int(x) for x in line.split("\t")[2:4] ] ]
		
		# Do the same for target names:
		tname=line.split("\t")[5]
		Tseqdic[tname]= [ [ int(x) for x in line.split("\t")[7:9] ] ]
		
# Sort the Qdictionary: 
for key, val in Qseqdic.items():
	# Sort each list of values.
	val.sort()
	
# Same for Tdict. :
for key, val in Tseqdic.items():
	val.sort()



# View the created dictionaries:		
print(Qseqdic)
# ~ print("\n"+"="*80+"\n")
# ~ print(Tseqdic)








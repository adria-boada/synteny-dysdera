#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# overlap_paf_parse.py 
#
# 2022-03-18  <adria@molevol-OptiPlex-9020>

""" Functions used by the paf_parser.py script.

"""

import sys

# Instruccions respecte els arguments necessaris per cridar l'script:
#~ if len(sys.argv) < 2:
    #~ sys.exit('\nCrit script: script.py <input> <output>\n')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 # Funció per extreure parelles de solapats que es poden concatenar.
def couples_overlapping (interval_list: "A list of sorted points from *point_list*"):
	
    """ Acquire and return indices of overlapping intervals from 
	*point_list*. Returns couples of overlapping intervals in the form 
	of sublists. Only does one sweep to find pairs; to find multiple 
	overlaps an iterative approach is required. """

    # Create point_list from interval_list:
    point_list = []

    for i in range(0, len(interval_list)):
        # Left point == 0.
        point_list += [ [interval_list[i][0], 0, i] ]
        # Right point == 1.
        point_list += [ [interval_list[i][1], 1, i] ]
	
    # Sort the list by position (first val in sublists):
    point_list.sort()


    # State variables for algorithm:
    currentOpen = -1
    answer = []
	
    # for each point in the point_list:
    for i in range(0, len(point_list)):
	
        # Si el punt actual és 'Left; opens an interval'
        if point_list[i][1] == 0:
			
            # Si no hi ha cap interval actualment obert:
            if currentOpen == -1:
                # 'Entra' a l'interval 'i'.
                currentOpen = point_list[i][2]
		
            # Si es trobava dins un interval anterior:
            else:
                # Index del nou interval:
                index = point_list[i][2]
                # Afageix la parella a la resposta:
                answer += [ [currentOpen, index] ]
                # Tanca l'interval previ
                currentOpen = -1

        # Si troba un punt que és 'Right'; que tanca interval:
        else:
            # Surt del currentOpen.
            currentOpen = -1

    return answer


def toy_examples():

    """ Returns a list wich can be used as examples to test parsing functions.

    """
    readme = "This list contains an example in each item.\n0: Explanation.\n1: Intervals list.\n2: Point list."
    toy_paf_plist = [ [10, 0, 0], [15, 0, 1], [20, 1, 0], [25, 1, 1], [30, 0, 2], [35, 1, 2], [40, 0, 4], [40, 0, 3], [45, 1, 3], [50, 1, 4] ]
    toy_paf_intervals = [ [10, 20], [15, 25], [30, 35], [40, 45], [40,50] ]

    return ([ readme, toy_paf_intervals, toy_paf_plist ])


 # Funció per extreure la totalitat dels índexos dels intervals solapats.
def thorough_indices (interval_list: "A list of sublists with paired points, which represent intervals."):

    """ Return all overlapping indices in a simple list

    point_list -> a list with the format:
    [ [position, Left(0)/Right(1), Interval_Index] ]

    """
    
    # Create point_list from interval_list:
    point_list = []

    for i in range(0, len(interval_list)):
        # Left point == 0.
        point_list += [ [interval_list[i][0], 0, i] ]
        # Right point == 1.
        point_list += [ [interval_list[i][1], 1, i] ]
	
    # Sort the list by position (first val in sublists):
    point_list.sort()

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
    
                # Si l'interval currentOpen no ha sigut encara afegit:
                if (not added):
                    # Fes-ho:
                    answer += [currentOpen]
                    added = True

                    # Mira quin dels dos intervals que s'estan comparant
                    # te la cua més llarga, i per tant solaparà amb
                    # més intervals (cua = Right point).
                    if (interval_list[currentOpen][1] < interval_list[index][1]):
                        currentOpen = index
                        added = True
    
        # Si troba un punt que és 'Right'; que tanca interval:
        else:
            # I és el punt esquerra de l'interval actual:
            if point_list[i][2] == currentOpen:
                # Tanca l'interval previ
                currentOpen = -1
                added = False
  
    return answer



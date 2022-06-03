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
    readme = "This list contains an example in each item.\n0: Explanation.\n1: Intervals list.\n2: Point list.\n3: real CIGAR string. \n4: short artificial CIGAR string"
    toy_paf_plist = [ [10, 0, 0], [15, 0, 1], [20, 1, 0], [25, 1, 1], [30, 0, 2], [35, 1, 2], [40, 0, 4], [40, 0, 3], [45, 1, 3], [50, 1, 4] ]
    toy_paf_intervals = [ [10, 20], [15, 25], [30, 35], [40, 45], [40,50] ]
    cigarro = "28M3D86M1I64M1I181M3I69M33D113M2I82M6D20M3D99M7I43M231D43M1I103M1D220M3I92M2D252M1I145M382D42M6I114M6I66M3I37M2D149M3D13M4D287M1D18M3I87M2D55M2D3M29D106M7D19M2D9M6D24M1D24M1I154M1D8M10I98M2D45M1D51M3D107M54D56M7I101M4I24M10D745M3D17M1I240M1D189M2D18M4I80M3D775M6I58M10D10M9D53M1D79M13I196M1I441M1I108M1D62M1D27M1D131M2D33M1I361M1D522M2D123M1D228M1I17M4I54M1D100M1I189M1D257M1I56M1I441M7D122M2I55M3D50M4I4M11I54M1I22M38I4M2I299M1D130M2I220M4D212M1D280M1839D1012M1I71M6I19M3I117M8D63M1I134M4D11M1I197M8D23M13D67M1D89M8I222M6I41M1I42M15I47M1D84M6I15M8D238M8I154M6D25M2D12M1000D194M8I70M2I533M1I96M1D79M4D27M7I24M1D100M4D28M3I20M3I78M26D51M2I191M4I160M1D358M3D84M176D123M8I203M161D437M4D39M58D33M1D75M2I31M1I328M6D136M2I202M1I136M2D34M79D336M2I121M1D126M1I63M185D204M1D34M1D49M11I130M2I11M1D21M2D65M132D17M9I145M3I129M10I36M7I70M2D19M6I619M1D57M1D20M2I44M2D7M12D138M1I12M1D5M6D8M5D2M8D55M1I344M1D177M7D96M3530D66M1I28M3D33M1I95M367D101M1D19M1D30M1D37M3D14M1I170M8D2M5D73M1D53M34D82M6D8M1D537M1D39M32D91M6D68M2I125M2D52M19I151M1D14M4I9M6D24M2D39M1D27M1D60M5I66M2D150M1I45M1I74M4I150M1D8M10D220M1I554M3I109M2I36M2I27M2I266M5D176M1D47M5I20M1D75M5D70M2I179M"
    artificial_cig = "100M10I100M10D100M"

    return ([ readme, toy_paf_intervals, toy_paf_plist, cigarro, artificial_cig ])


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


# Accepts a CIGAR string and calculates total and percent M, I and D. 
def cig_analysis (cig: "CIGAR string"):
    """ Returns a dictionary with 'tM', 'tI','tD'.
    Accepts a string with Matches, Deletions and Insertions 
    (and no other type).

    """
    # Input:
    c = cig

    # Split string by adding ' ' to the end of any letter:
    for i in 'MID':
        c = c.replace(i, i+' ')

    # Transform spaces into list separations:
    c = c.split(' ')

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


if __name__ == '__main__':
    print (cig_analysis(toy_examples()[3]))



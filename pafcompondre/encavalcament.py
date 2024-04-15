#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# encavalcament.py
#
# 13 d’abr. 2024  <adria@molevol-OptiPlex-9020>

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def remove_overlapping_with_score(intervals: list):
    """
    Take a list of intervals (pairs of begin and end coordinates), coupled with
    a scoring index (mapping quality, Smith-Waterman score, etc.) and remove
    overlapping regions, prioritising the removal of lower scoring intervals.

    Parameters
    ----------

    intervals : list
    A list of intervals which, in turn, are sublists of index,  begin and end
    coordinates, and score. E.g.

    [ (idx1, beg1, end1, score1), ... (idxN, begN, endN, scoreN) ]

    Returns
    -------

    A list of sublists. The first item in the sublist is the provided index for
    intervals, whereas the second item is the length of the nonoverlapping
    annotation, in basepairs. E.g.

    [ (idx1, len1), ... (idxN, lenN) ]
    """
    # Inicialitza un diccionari per guardar-hi els resultats.
    answer = dict()
    # Crea una llista de punts a partir d'una llista d'intervals. Exemple:
    # Llista de punts → ( [beg1, 0, idx1, score1], [end1, 1, idx1, score1],
    #                     [beg2, 0, idx2, score2], [end2, 1, idx2, score2], ...
    #                     [begN, 0, idxN, scoreN], [endN, 1, idxN, scoreN] )
    llista_punts = []
    for i in range(0, len(intervals)):
        # Esquerra o «beg» punts etiquetats amb "0".
        llista_punts.append(list([intervals[i][1],
                                  0,
                                  intervals[i][0],
                                  intervals[i][3], ]))
        # Dreta o «end» punts etiquetats amb "1".
        llista_punts.append(list([intervals[i][2],
                                  1,
                                  intervals[i][0],
                                  intervals[i][3], ]))
        # Crea una entrada al diccionari on sumar-hi els resultats.
        answer[intervals[i][0]] = 0
    # Ordena la llista de punts segons la posició (primer item de subllistes).
    llista_punts = sorted(llista_punts, key=lambda x: x[0])
    # Inicialitza altres variables de l'algorisme.
    currentBest = -1
    intervals_oberts = list()

    # Comença iterant a través de la llista de punts. Descompta bp encavalcats
    # de les anotacions amb pitjor puntuació.
    for punt in llista_punts:
        # Si el loop arriba a un punt esquerre "0", obre interval.
        if punt[1] == 0:

            # Si no hi ha altres intervals oberts:
            if currentBest == -1:
                # El millor interval és aquest mateix.
                currentBest = int(punt[2])
                currentBegin = int(punt[0])
                currentScore = int(punt[3])

            # Altrament,ja hi havia un interval encetat.
            elif currentScore > punt[3]:
                # No tanquis currentBest, perquè el pròxim interval és d'una
                # puntuació inferior a l'actual.
                intervals_oberts.append(list([
                    int(punt[2]),       # ID
                    int(punt[3]), ]))   # Puntuació.
                # Ordena intervals oberts segons puntuació (primer camp).
                intervals_oberts = sorted(intervals_oberts,
                                          key=lambda x: x[1],
                                          reverse=True)
            else:
                # El próxim interval és d'una puntuació superior. Computa
                # llargada 'currentBest' fins ara.
                llargada_interval = int(punt[0] - currentBegin)
                answer[currentBest] += llargada_interval
                # Afageix "currentBest" als intervals prèviament encetats.
                intervals_oberts.append(list([
                    int(currentBest),       # ID
                    int(currentScore), ]))  # Puntuació.
                # Ordena intervals oberts segons puntuació (primer camp).
                intervals_oberts = sorted(intervals_oberts,
                                          key=lambda x: x[1],
                                          reverse=True)
                # Actualitza currentBest, enceta pròxim interval.
                currentBest = int(punt[2])
                currentBegin = int(punt[0])
                currentScore = int(punt[3])

        # Altrament, el bucle es troba en un punt a la dreta "1", que tanca
        # un interval actualment encetat/obert.
        elif punt[2] == currentBest:
            # Computa la llargada de l'interval currentBest fins a tancar-se.
            # RepeatMasker i minimap2 diria que són "0-based" (suma 1 a la
            # diferència entre principi i final per obtenir llargada).
            llargada_interval = int(punt[0] - currentBegin +1)
            answer[punt[2]] += llargada_interval
            # Tanca "currentBest". Revisa si hi ha intervals oberts que puguin
            # prendre-li el lloc.
            if len(intervals_oberts) == 0:
                currentBest = -1
            else:
                currentBest, currentScore = intervals_oberts.pop(0)
                # Calcula principi a partir de l'interval que acabem de tancar.
                currentBegin = punt[0] +1

        else:
            # El punt que volem tancar no és actualment obert. S'ha d'eliminar
            # aquest punt dels intervals oberts.
            intervals_oberts.remove(list([
                int(punt[2]),      # ID
                int(punt[3]), ]))  # Puntuació

    return list(zip(answer.keys(), answer.values()))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


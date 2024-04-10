#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# remescla.py
#
# 10 d’abr. 2024  <adria@molevol-OptiPlex-9020>

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def mescla_i_redueix_fitxer_enorme(
    cami_entrada: str, cami_sortida: str,
    n_linies: int=2500 ):
    """
    Redueix la mida d'un fitxer enorme tot recollint línies aleatoriament del
    fitxer original. Se'n recolliràn tantes com el paràmetre `n_linies`.

    Aquesta funció es pot usar per crear una mostra d'un fitxer enorme, amb
    l'intenció de facilitar-ne les proves d'anàlisis computacionalment
    intensius. La idea prové de <https://stackoverflow.com/q/60623645/21787735>.
    Aquesta funció no exporta el nombre exacte de línies especificades, sinó una
    quantitat aproximada de línies (degut al seu comportament aleatori).

    Input
    =====

    cami_entrada: Un camí que indiqui quin fitxer llegir com a entrada.

    cami_sortida: Un camí que indiqui a on guardar el resultat.

    n_linies [opcional]: Número de línies aproximades del fitxer resultant.

    Output
    ======

    Escriu a `cami_sortida` un subconjunt de línies provinents del fitxer
    localitzat a `cami_entrada`.
    """
    # Computa el nombre total de línies del fitxer d'entrada.
    with open(cami_entrada, "r") as fitxer:
        for t_linies, linia in enumerate(fitxer):
            pass
    # Reobre el fitxer altra vegada, retornant el punter al principi. A més a
    # més, obre el cami de sortida en mode "a"; d'aquesta manera evita eliminar
    # fitxers preexistents i afageix resultat al final.
    with open(cami_entrada, "r") as fitxer_entrada, \
            open(cami_sortida, "a") as fitxer_sortida:
        # Computa la probabilitat d'imprimir una línia de manera que s'obtinguin
        # aproximadament `n_linies` resultants.
        prob = n_linies / t_linies
        [fitxer_sortida.write(linia) for linia in fitxer_entrada
            if random.random() < prob]

    return fitxer_sortida

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


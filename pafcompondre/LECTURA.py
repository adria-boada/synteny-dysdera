#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# lectura.py
#
# 11 d’abr. 2024  <adria@molevol-OptiPlex-9020>

# Printing colourful messages.
from . import missatges
# Accounting for overlapping intervals.
from . import encavalcament as encav

# Dataframes manipulation, akin to the R project
import pandas as pd
# CIGAR pattern matching
import re

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class Mapping(object):
    """
    Create a Python object from a PAF mapping file. A documentation of the
    columns in a PAF file can be found in the addendum of the manual of
    `minimap2`:

    <https://lh3.github.io/minimap2/minimap2.html#10>
    # accessed 16-april-2024

    In summary, we have used the following column names:

    "matches" + "mismatches_and_gaps" == "ali_len" == \
        "cigar_M" + "cigar_D" + "cigar_I"

    "mismatches" + "matches" == "cigar_M"

    Parameters
    ----------

    path_to_paf : str
    The path to a PAF file. The file will be read and transformed into a
    `pd.Dataframe`.

    memory_efficient : bool, default False
    Whether memory expensive columns should be loaded or left discarded. For
    instance, CIGAR and difference strings will not be loaded when True. Might
    prevent Python to run out of memory. Set to `True` if the Python process is
    being killed.
    """
    def __init__(self, filename: str, df: pd.DataFrame,
                 genome_subsets_patterns: dict=None):
        # Revisa si s'han entregat patrons dins el mètode de fàbrica!
        if not genome_subsets_patterns:
            # Patrons predeterminats, genèrics, d'exemple.
            genome_subsets_patterns = {
                "genome": ".",
                "autosomes": "chr\d", }

        return None

    @classmethod
    def from_PAF(cls, path_to_paf: str, memory_efficient: bool=False):
        """
        `Factory method'. It accepts a path to a PAF file, and returns an
        instance of the `Mapping' class.

        Parameters
        ----------

        path_to_paf : str
        The path to a PAF file. The file will be read and transformed into a
        `pd.Dataframe`.

        memory_efficient : bool, default False
        Whether memory expensive columns should be loaded or left discarded. For
        instance, CIGAR and difference strings will not be loaded when True. Might
        prevent Python to run out of memory. Set to `True` if the Python process is
        being killed.
        """
        # Revisa que el camí fins al fitxer PAF existeixi.
        try:
            with open(path_to_paf):
                pass
            missatges.Estat("Opening the file " + str(path_to_paf) + ".")
        except FileNotFoundError:
            missatges.Error("The file " + str(path_to_paf) +
                            " was not found.")
            return None

        # Emmagatzema el nom del fitxer. Elimina els directoris que precedeixen
        # el camí `path_to_paf`.
        paf_filename = path_to_paf.strip("\n").split("/")[-1]
        # Elimina l'extensió PAF del fitxer.
        fitxer_separat = paf_filename.split(".")
        if fitxer_separat[-1].lower() == "paf":
            paf_filename = ".".join(fitxer_separat[:-1])
        else:
            paf_filename = ".".join(fitxer_separat)

        # Llegeix el fitxer PAF, incorporant-lo com a `pandas.DataFrame()`
        missatges.Estat("Reading the first 12 columns of the file "+
                        "(always present).")
        df_main = cls._read_main_paf_file(path_to_paf)
        # Si `memory_efficient` està desactivat (False), llegeix la resta de
        # columnes opcionals que ocupen considerablement més memòria.
        if not memory_efficient:
            missatges.Estat("Reading optional columns beyond the 12th field.")
            df_main = cls._read_optional_paf_file(df_main, path_to_paf)

        return cls(filename=paf_filename, df=df_main)

    @classmethod
    def _read_main_paf_file(cls, path_to_paf: str):
        """
        """
        # PAF may optionally have additional fields, i.e. an irregular number of
        # columns (ragged TSV).
        # This function initially reads columns which are always present.
        colnames = ["Qname",    # Name of a query scaffold/chr.
                    "Qlen",     # Length of above sequence.
                    "Qstart",   # Start of alignment in Qname.
                    "Qend",     # End of alignment in Qname.
                    "strand",   # Complement DNA seq.
                    # Same variables for target genome/assembly.
                    "Tname", "Tlen", "Tstart", "Tend",
                    "matches",  # Number of matching bases.
                    "ali_len",  # Sum of matches, mismatches and gaps.
                    "mapQ"]     # Mapping quality.
        coltypes = {"Qname": "category", "Qlen": "int", "Qstart": "int",
                    "Qend": "int", "strand": "category", "Tname": "category",
                    "Tlen": "int", "Tstart": "int", "Tend": "int",
                    "matches": "int", "ali_len": "int", "mapQ": "int"}

        df = pd.read_table(path_to_paf, header=None,
            # The regex "\s+" stands for
            # "one or more blank spaces" (includes spaces and TABs).
            sep="\s+",
            # Load the appropiate (always present) columns and no more.
            usecols=list(range(12)), # [0, 1, 2, etc., 10, 11]
            # Use the following column names (from first to last column).
            names=colnames,
            # Specify the dtypes of each column.
            dtype=coltypes)

        return df

    @classmethod
    def _read_optional_paf_file(cls, df: pd.DataFrame, path_to_paf: str):
        """
        """
        # Crea una llista amb els camps opcionals. Són camps que podrien
        # trobar-se en totes, alguna o cap filera. Veieu el manual de referència
        # per una explicació de cada etiqueta/columna:
        # <https://lh3.github.io/minimap2/minimap2.html#10>
        # Recordeu que ordeno les columnes manualment segons "importància"
        # subjectiva i personal.

        # Drecera fins al dataframe principal, evitant paramatritzar-lo.
        df = self.df

        etiqueta_a_columna = {"tp": [], "NM": [], "nn": [], "dv": [], "de": [],
                              "cg": [], "cs": [], "SA": [], "ts": [], "cm": [],
                              "s1": [], "s2": [], "MD": [], "AS": [], "ms": [],
                              "rl": [], "zd": [], }
        # A més, interpreta la cadena CIGAR (si existeix).
        cigar_columnes = {"cig_deletions": [], "cig_insertions": [],
                          "cig_matches": [], "cig_compressed": [],
                          "indel_list": [], }
        # Desplaça't a través de les línies del fitxer PAF i omple les llistes
        # dels diccionaris anteriors.
        with open(path_to_paf) as fitxer_paf:
            for count, line in enumerate(fitxer_paf):
                # Obté una llista amb tots els camps a partir del 12é.
                line = line.strip("\n").split("\t")[12:]
                while line:
                    camp = line.pop(0)
                    # L'etiqueta dels camps opcionals es troba a
                    # la primera parella de caràcters.
                    etiqueta = camp[0:2]
                    # El valor del camp es troba a partir del 5é caràcter.
                    valor = camp[5:]
                    etiqueta_a_columna[etiqueta].append(valor)
                # Omple amb "None" els camps que no s'han trobat dins la filera:
                for key, val in etiqueta_a_columna.items():
                    if len(val) < count+1:
                        etiqueta_a_columna[key].append(None)
                # Si aquesta línia contenia una cadena CIGAR, interpreta-la:
                if etiqueta_a_columna["cg"][-1] != None:
                    last_cigar = etiqueta_a_columna["cg"][-1]
                    resum_cigar = interpret_cigar_string(last_cigar)
                    cigar_columnes["cig_deletions"].append(resum_cigar["D"])
                    cigar_columnes["cig_insertions"].append(resum_cigar["I"])
                    cigar_columnes["cig_matches"].append(resum_cigar["M"])
                    cigar_columnes["cig_compressed"].append(resum_cigar["compressed"])
                    cigar_columnes["indel_list"].append(resum_cigar["indels"])
                # Si no contenia CIGAR, afageix "None" per representar buit a
                # les dades finals.
                else:
                    for key in cigar_columnes.keys():
                        cigar_columnes[key].append(None)

        # Tradueix aquestes curtes etiquetes a noms de columna més descriptius:
        reanomenar_etiquetes = {"tp": "type_aln", "cm": "num_minimizers",
                          "s1": "chaining_score", "s2": "second_chain_score",
                          "NM": "mismatches_and_gaps",
                          "MD": "regenerate_refseq", "AS": "DP_ali_score",
                          "SA": "list_supp_ali", "ms": "DP_max_score",
                          "nn": "ambiguous_bases", "ts": "transcript_strand",
                          "cg": "CIGAR_string", "cs": "diff_string",
                          "dv": "divergence", "de": "gap_compr_diverg",
                          "rl": "length_repetitive_seeds",
                          "zd": "zd_unknown", }
        # Finalment, afageix els camps opcionals al marc de dades.
        for etiqueta, llista_valors in etiqueta_a_columna.items():
            df[reanomenar_etiquetes[etiqueta]] = llista_valors
        # Afageix l'interpretació / anàlisi del CIGAR:
        for nom_columna, llista_valors in cigar_columnes.items():
            df[nom_columna] = llista_valors

        # Assegura que els «dtypes» d'aquestes columnes són correctes.
        columnes_enters = ["mismatches_and_gaps", "ambiguous_bases",
                           "num_minimizers", "chaining_score", "DP_ali_score",
                           "DP_max_score", "length_repetitive_seeds",
                           "cig_deletions", "cig_insertions", "cig_matches",
                           "cig_compressed"]
        df[columnes_enters] = df[columnes_enters].apply(pd.to_numeric,
                                                        downcast="integer")
        columnes_decimals = ["gap_compr_diverg", "second_chain_score",
                             "zd_unknown"]
        df[columnes_decimals] = df[columnes_decimals].astype("float")
        df["type_aln"] = df["type_aln"].astype("category")

        # Elimina columnes buides (totes les entrades amb valors nuls).
        df = df.dropna(axis="columns", how="all")

        # Per evitar un error estrany de "SettingWithCopyWarning", faig una
        # còpia abans de computar noves columnes del marc de dades???
        df = df.copy()
        # Havent assignat el «dtype» numèric, computa altres paràmetres i
        # afageix-los al marc de dades:
        df["blast_identity"] = df["matches"] / df["ali_len"]
        # Computa altres valors a partir del CIGAR, si és present.
        if ("cig_matches" in df.columns) and \
           ("cig_compressed" in df.columns) and \
           ("cig_deletions" in df.columns) and \
           ("cig_insertions" in df.columns) and \
           ("mismatches_and_gaps" in df.columns):
            df["gap_compr_identity"] = (
                df["matches"] / (df["cig_matches"] + df["cig_compressed"]))
            df["mismatches"] = \
                df["mismatches_and_gaps"] - \
                (df["cig_insertions"] + df["cig_deletions"])

        return df

    def _list_intervals():
        """
        """

    def 

def interpreting_cigar_string(cigar_string: str):
    """
    Interpret a CIGAR string; compute its total amount of "M" (matches), "I"
    (insertions to query) and "D" (deletions to query).

    Parameters
    ----------

    cigar_string : str
    A string of letters, each paired with a number. They symbolize the gap
    composition of an alignment/mapping. Accepts a string with the letters
    "M", "I" and "D". Each letter must be preceded by a number, which is the
    length of the gap/match. Take a look at the specifications of the "SAM"
    format for a more in-depth explanation:

    <https://samtools.github.io/hts-specs/SAMv1.pdf>
    # accessed 12-april-2024

    Returns
    -------

    A dictionary with the sum of "M", "I" and "D" fragments. Moreover, it
    returns the amount of insertions and deletions instead of their lengths
    in basepairs. It is useful to compute "compressed gap identity". For a
    discussion on "compressed gap identity", see:

    <https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity>
    # accessed 12-april-2024
    """
    # Crea una llista de llargada (bp) d'indels
    # identificant el patró "\d+[ID]".
    indels = re.findall(r"\d+[ID]", cigar_string)
    indels = re.findall(r"\d+", "".join(indels))
    indels = [int(i) for i in indels]
    # Separa la cadena CIGAR per mitjà de
    # concatenar espais a les lletres "MID":
    for i in "MID":
        cigar_string = cigar_string.replace(i, i+" ")
    cigar_list = cigar_string.strip("\n").split(" ")
    # Inicialitza un diccionari de retorn:
    answer = {"M": 0, "I": 0, "D": 0, "indels": indels}
    # Compta i emmagatzema resultats a `answer`.
    for i in cigar_list:
        if "M" in i:
            answer["M"] += int(i[:-1])
        elif "D" in i:
            answer["D"] += int(i[:-1])
        elif "I" in i:
            answer["I"] += int(i[:-1])
    # Prén en consideració la quantitat de gaps en comptes de la seva
    # llargada (és a dir, comprimeix els gaps a la llargada igual a 1).
    answer["compressed"] = cigar_string.count("D") + cigar_string.count("I")

    return answer

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# lectura.py
#
# 11 d’abr. 2024  <adria@molevol-OptiPlex-9020>

help_msg = """
Describe here the goal, input, output, etc. of the script
as a multi-line block of text.
"""

import sys
from . import missatges

import pandas as pd

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class Mapping:
    """
    Create a Python object from a PAF mapping file. A documentation of the
    columns in a PAF file can be found in the addendum of the manual of
    `minimap2`:

    <https://lh3.github.io/minimap2/minimap2.html#10>

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
    def __init__(self, path_to_paf: str, memory_efficient: bool=False):

        self.df = self.read_paf_file(path_to_paf)

        return None

    def read_paf_file(self, path_to_paf: str):
        # PAF may optionally have additional fields, i.e. an irregular number of
        # columns (ragged TSV).
        # First, read the columns that are always present.
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

        ##> missatges.Estat("Reading the first 12 columns of the file "+
        ##>                "(always present)")
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

    def interpret_cigar_string(self, cigar_string: str):
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
        # accessed 12-04-2024

        Returns
        -------

        A dictionary with the sum of "M", "I" and "D" fragments. Moreover, it
        returns the amount of insertions and deletions instead of their lengths
        in basepairs. It is useful to compute "compressed gap identity". For a
        discussion on "compressed gap identity", see:

        <https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity>
        # accessed 12-04-2024
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

# Si es crida com a script:
if __name__ == '__main__':

    ### PARSE ARGUMENTS ###
    import argparse
    parser = argparse.ArgumentParser(description=help_msg,
        # make sure the 'help_msg' is not automatically
        # wrapped at 80 characters (manually assign newlines).
        formatter_class=argparse.RawTextHelpFormatter)
    # file-name: positional arg.
#    parser.add_argument('filename', type=str, help='Path to ... file-name')
    # integer argument
#    parser.add_argument('-a', '--numero_a', type=int, help='Paràmetre "a"')
    # choices argument
#    parser.add_argument('-o', '--operacio', 
#            type=str, choices=['suma', 'resta', 'multiplicacio'],
#            default='suma', required=False,
#            help='')

    args = parser.parse_args()
    # call a value: args.operacio or args.filename.



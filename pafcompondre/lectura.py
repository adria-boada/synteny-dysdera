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
import numpy as np
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

    filename : str
    The name of the mapping file.

    df : pd.DataFrame
    The contents of the mapping file formatted as a pandas DataFrame. See the
    factory method `cls.from_PAF', which creates a pandas df from a PAF file.

    genome_subsets_patterns : dict
    A dictionary where keys are names of genome subsets and values are
    corresponding patterns which can be used to filter the main dataframe.
    For example:

        genome_subsets_patterns = dict({
            "genome": ".",
            "autosomes": "chr\\d",
            "sex_chr": "chrx",
            "minor_scaff": "scaff", }
    """
    def __init__(self, filename: str, df: pd.DataFrame,
                 genome_subsets_patterns: dict=None):
        # Emmagatzema `df' i `filename'.
        self.filename = filename
        self.df = df
        # Revisa si s'han entregat patrons dins el mètode de fàbrica!
        if not genome_subsets_patterns:
            # Patrons predeterminats, genèrics, d'exemple.
            self.genome_subsets_patterns = {
                "genome": r".",
                "autosomes": r"chr\d",
                "sex_chr": r"chrx",
                "minor_scaff": r"scaff|ctg", }

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
        prevent Python to run out of memory. Set this parameter to `True` if the
        Python process ends up killed.
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

        missatges.Estat("Finished reading the PAF file.")
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
            sep=r"\s+",
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

    def _list_intervals(self, df: pd.DataFrame, column_sequid: str,
                        column_start: str, column_end: str):
        """
        Create a dictionary with all the intervals found in each sequid of
        either the queried or targeted genome.

        Parameters
        ----------

        df : pd.DataFrame
        The main dataframe created from the PAF file.

        column_sequid : str
        Either "Qname" or "Tname", where sequence IDs can be found.

        column_start : str
        Either "Qstart" or "Tstart", where interval openings can be found.

        column_end : str
        Either "Qend" or "Tend", where interval endings can be found.
        """
        # Inicialita un diccionari de retorn. Emmagatzemarà intervals indicant
        # regions incloses al mapatge. Una llista d'intervals per cromosoma.
        # P. ex. {"sequid": ((beg1, end1), (beg2, end2), ... (begN, endN)) }
        intervals_per_seq = dict()
        for seq in df[column_sequid].unique():
            df_seq = df.loc[df[column_sequid] == seq]
            intervals = list( zip(
                list(df_seq.index),
                list(df_seq[column_start]),
                list(df_seq[column_end]),
                list(df_seq["mapQ"]), ))
            # Ordena intervals segons la coordenada inicial.
            intervals = sorted(intervals, key=lambda x: x[0])
            intervals_per_seq[seq] = intervals

        return intervals_per_seq

    def compute_nonoverlap_len(self, df: pd.DataFrame):
        """
        Compute how many basepairs of each mapping are not overlapping with
        other mappings. Resolve conflicts by retaining the basepairs of higher
        mapping quality alignments.

        Parameters
        ----------

        df : pd.DataFrame
        The main dataframe created from the PAF file. Does not have the columns
        `Qalgorbp' nor `Talgorbp'.

        Returns
        -------

        A dict with two values: (1) the inputted `df' with two newly added
        columns, "Qalgorbp" and "Talgorbp", and (2) a dictionary which relates each
        sequence to a list of their nonoverlapping, aligned regions.

        The newly added columns convey the length of the nonoverlapping aligned
        regions, in basepairs. They're obtained from `Qstart' or `Tstart' and
        `Qend' and `Tend' columns.
        """
        # No sobrescriguis el "df". Crea'n una drecera.
        df = df.copy()
        # Inicialita una llista on guardar-hi els resultats de l'algorisme.
        answer = {"Qalgorbp": list(), "Talgorbp": list()}
        index = {}
        intervals_per_seq_sans_encavalcament = {"Qinterval": dict(),
                                                "Tinterval": dict(), }
        # Repeteix el procediment per a cada base de dades "query" i "target".
        for db in ["Q", "T"]:
            llista_extensible = list()
            # Obté un diccionari amb totes les seqüències de la base de dades,
            # juntament a la seva llista de coordenades d'intervals mapats.
            intervals_per_seq = self._list_intervals(
                df=df,
                column_sequid=str(db) + "name",
                column_start=str(db) + "start",
                column_end=str(db) + "end", )
            # Travessa el diccionari de seqüències.
            for seq, intervals in intervals_per_seq.items():
                # Per cada seqüencia, elimina nucleòtids encavalcats.
                sans_encavalcament = \
                    encav.remove_overlapping_with_score(intervals)
                # Extén la resposta afegint-hi parelles "(index, llargada)".
                llista_extensible.extend(sans_encavalcament["algorbp"])
                # Amplia el diccionari d'intervals sense encavalcar.
                intervals_per_seq_sans_encavalcament \
                    [str(db+"interval")][seq] = \
                    sans_encavalcament["intervals"]
            answer[str(db) + "algorbp"] = [x[1] for x in llista_extensible]
            index[db] = [x[0] for x in llista_extensible]

        # Crea dues `Series' amb les llargades sense encavalcar (en bp.)
        column_Qalgorbp = pd.Series(answer["Qalgorbp"],
                               name="Qalgorbp", index=index["Q"],
                               dtype=int).sort_index()
        column_Talgorbp = pd.Series(answer["Talgorbp"],
                               name="Talgorbp", index=index["T"],
                               dtype=int).sort_index()
        # Enganxa les dues columnes "Qalgorbp" i "Talgorbp" al df.
        df = df.join([column_Qalgorbp, column_Talgorbp, ])

        # Finalment, crea dues pd.Series a partir de les parelles anteriors.
        return {"df": df,
                "intervals_sans_encav": \
                    intervals_per_seq_sans_encavalcament, }

    def genome_subset_length(self, df: pd.DataFrame, column_sequid: str,
                             column_len: str, ):
        """
        """
        # Obtén una pd.Series amb les llargades de totes les seqüències. Evita els
        # duplicats que podrien sorgir d'una simple i directa operació "unique()"
        # a la columna "Qlen" o "Tlen".
        series_seqlen = df[[column_sequid, column_len]] \
            .groupby(column_sequid, observed=True) \
            .agg("min") \
            [column_len]

        return series_seqlen

    def evaluate_overlapping_alignment(self, df: pd.DataFrame,
                                       genome_subset_patterns: dict):
        """
        Parameters
        ----------

        df : pd.DataFrame
        The main dataframe created from the PAF file. Does have the columns
        `Qalgorbp' and `Talgorbp' obtained by running the function
        `self.compute_nonoverlap_len()'.

        genome_subset_patterns : dict
        The keys are genome subsets, whereas the values are patterns that will
        be matched against sequence names to filter the main dataframe. For
        instance,

        {"genome": ".", "autosomes": "chr\\d", "sex_chr": "chrx",
         "minor_scaff": "scaff|ctg", }

        The pattern for "genome" encompasses any sequence name ("." matches
        everything). The pattern "chr\\d" matches "chr" followed by amy digit.
        The pattern "scaff|ctg" matches either "scaff" or "ctg".
        """
        # Revisa que `df' conté ambdós columnes "Qalgorbp" i "Talgorbp". Si no,
        # demanar que computi `self.compute_nonoverlap_len'.
        if ("Qalgorbp" not in df.columns) or \
           ("Talgorbp" not in df.columns):
            missatges.Error("The function "+
                            "`self.evaluate_overlapping_alignment()' "+
                            "takes as the `df' parameter the output of "+
                            "`self.compute_nonoverlap_len()['df'].")
            return None

        answer = {"Subset": [], "Kind": [], "Length (#bp)": [],
                  "Length (%naive)": [], "Length (%subset)": [], }

        for db in ["Q", "T"]:
            for subset, pat in genome_subset_patterns.items():
                # The filtered df is the main df where the column
                # "str(db+'name')" contains the pattern, case-insensitively.
                df_filtered = df.loc[df[str(db+"name")].str.contains(
                    pat, case=False)]
                # Compute the length of alignments naively, end - start +1
                series_naive = df_filtered[str(db+"end")] - \
                    df_filtered[str(db+"start")] +1
                series_algor = df_filtered[str(db+"algorbp")]
                # Recull les llargades de totes les seqüències (llargada del
                # subset genòmic). Paràgraf una mica complicat per
                # aconseguir-ho.
                series_seqlen = self.genome_subset_length(
                    df=df_filtered,
                    column_sequid=str(db+"name"),
                    column_len=str(db+"len"), )

                # Emmagatzema aquests resultats per a crear-ne un pd.DataFrame.
                answer["Subset"].extend([str(db + "." + subset)] *2)
                answer["Kind"].extend(["naive", "algor"])
                answer["Length (#bp)"].extend([sum(x) \
                                        for x in (series_naive, series_algor)])
                answer["Length (%naive)"].extend([sum(x)/sum(series_naive) \
                                        if sum(series_naive) != 0 else 0 \
                                        for x in (series_naive, series_algor)])
                answer["Length (%subset)"].extend([sum(x)/sum(series_seqlen) \
                                        if sum(series_naive) != 0 else 0 \
                                        for x in (series_naive, series_algor)])

        return pd.DataFrame(answer)

    def len_mapped_and_unmapped(self, intervals: list, chr_end: int,
                                chr_begin: int=0, ):
        """
        Parameters
        ----------

        intervals : list
        A list that contains sublists [begin, end], representing the intervals
        of mapped regions within a single sequence/chromosome.

        chr_end : int
        The ending coordinate of the sequence/chromosome; sequence length.

        chr_begin : int
        The starting coordinate of the sequence/chromosome; zero.
        """
        # Crea una llista de punts a partir de la llista d'intervals.
        punts = []
        for i in range(0, len(intervals)):
            # Extén la llista afegint l'interval.
            punts.extend(list([intervals[i][0], intervals[i][1], ]))

        # Inicialitza variables.
        interdist, alig_lens = list(), list()
        # Calcula interdistància entre principi cromosoma i el 1r mapatge.
        principi = int(punts[0] - chr_begin)
        if principi > 0:
            interdist.append(principi)
        # Calcula 1r interval mapat (while loop necessita parelles "unmapped" i
        # "mapped", el 1r interval queda desaparellat).
        coord_beg, coord_end = [punts.pop(0) for i in range(2)]
        alig_lens.append(int(coord_end - coord_beg))
        # Itera a través de la llista de punts.
        while len(punts) > 0:
            # Computa distància entre alineaments, inter-alineaments.
            coord_beg = punts.pop(0)
            interdist.append(int(coord_beg - coord_end))
            # Computa llargada dels alineaments, intra-alineaments.
            coord_end = punts.pop(0)
            alig_lens.append(int(coord_end - coord_beg))
        # Finalment, calcula interdistància entre l'últim alineament i el final
        # del cromosoma.
        final = int(chr_end - coord_end)
        if final > 0:
            interdist.append(final)

        # Converteix les llistes a pd.Series. Haurien de ser més fàcils de
        # manipular.
        answer = {
            "mapped_coords": pd.Series(alig_lens, dtype=int),
            "unmapped_coords": pd.Series(interdist, dtype=int), }

        return answer

    def delimit_mapped_regions(self, intervals_per_seq: list,
                               genome_subset_patterns: dict):
        """
        """
        # Inicialitza un diccionari que emmagatzema llargada d'alineament i
        # distàncies entre mapatges (intervals de regions no mapades,
        # complementari als intervals de regions mapades).
        answer = dict()
        # Drecera fins al df principal de la classe.
        df = self.df
        # Computa pd.Series amb una llista de valors observats.
        for db in ("Q", "T"):
            answer[db] = dict()
            for seqtype, patt in genome_subset_patterns.items():
                # Filtra el df segons el patró del bucle.
                mask_seqtype = df[str(db + "name")] \
                    .str.contains(patt, case=False)
                df_filtered = df.loc[mask_seqtype].copy()
                if df_filtered.empty:
                    missatges.Avis("The pattern "+str(patt)+" does not "+
                                   "match any sequence in the column "+
                                   str(db)+"name.")
                    continue
                # Inicialitza un diccionari per a guardar-hi regions.
                answer[db][seqtype] = {
                    "mapped_coords": pd.Series(dtype="int"),
                    "unmapped_coords": pd.Series(dtype="int"), }
                # Crea una drecera fins a "answer".
                a = answer[db][seqtype]
                # Emmagatzema la llista de llargades dels alineaments, per
                # separat dels intervals mapats sense encavalcaments.
                a["ali_len"] = df_filtered["ali_len"].astype(int)
                # Emmagatzema "matches" i "mismatches" dels alineaments (la part
                # dels alineaments que no són "gaps").
                a["matches"] = df_filtered["matches"].astype(int)
                a["mismatches"] = df_filtered["mismatches"].astype(int)
                # Obtén un diccionari on es relacioni cada "seqtype" amb una
                # llista dels seus intervals mapats, eliminant encavalcaments.
                intervals = intervals_per_seq[str(db) + "interval"]
                for seq, intrv in intervals.items():
                    # Crea dues pd.Series amb les llargades de les regions
                    # mapades i no mapades.
                    seq_matches_pat = re.search(
                        pattern=patt, string=seq, flags=re.IGNORECASE)
                    if not seq_matches_pat:
                        continue
                    final_chr = df_filtered.loc[
                        df_filtered[str(db + "name")] == seq, str(db + "len")] \
                        .iloc[0]
                    series_regions = self.len_mapped_and_unmapped(
                        intervals=intrv, chr_end=final_chr)
                    a["mapped_coords"] = pd.concat([
                        a["mapped_coords"],
                        series_regions["mapped_coords"]], ignore_index=True)
                    a["unmapped_coords"] = pd.concat([
                        a["unmapped_coords"],
                        series_regions["unmapped_coords"]], ignore_index=True)

        return answer

    def indel_list_per_genome_subset(self, df: pd.DataFrame,
                                     patterns: dict, ):
        """
        Create a list of indel lengths per genome subset found in `patterns`.
        """
        indels_per_region = dict()
        ali_len_per_indel = dict()
        for db in ("Q", "T"):
            indels_per_region[db] = dict()
            ali_len_per_indel[db] = dict()
            for seqtype, patt in patterns.items():
                # Filtra df segons "db + name" (columna sequid) utilitzant el
                # patró adient.
                mask_seqtype = df[str(db + "name")].str.contains(patt, case=False)
                df_seqtype = df.loc[mask_seqtype]
                # Inicialitza una entrada al dict. i usa "list comprehension"
                # per a allisar una llista de subllistes.
                indels_per_region[db][seqtype] = list()
                [indels_per_region[db][seqtype].extend(i)
                 for i in df_seqtype["indel_list"].to_list()]
                # Computa la llargada d'alineament per a cada mida indel.
                num_indels = df_seqtype["indel_list"].map(len)
                ali_len_per_indel[db][seqtype] = \
                    np.repeat(df_seqtype["ali_len"], num_indels).to_list()

        return {"indels_per_regions": indels_per_region,
                "ali_len_per_indel": ali_len_per_indel, }

    def prepend_sequid_with_QT(self, series: pd.Series, label: str,
                               separador: str="."):
        """
        Prepend all strings within a pd.Series (either Qname or Tname) of
        categorical or object «dtype» with a label.

        Parameters
        ----------

        series : pd.Series
        A series with sequid values; either "Qname" or "Tname" columns of the
        main dataframe.

        label : str
        A label which will be added to the beginning of all sequids; for
        instance, "Q" for "Qname" and "T" for "Tname".

        Returns
        -------

        pd.Series with the aforementioned changes.
        """
        # Afageix etiquetes al principi.
        series = series.copy()
        tipus = series.dtype
        series = pd.Series([label] * len(series)).str.cat(series, sep=separador)
        series = series.astype(tipus)

        return series

    def prova(self):
        """
        """
        import seaborn as sns, matplotlib.pyplot as plt
        self.df, self.intervals_per_seq = self.compute_nonoverlap_len(self.df).values()
        m = self.delimit_mapped_regions(
            intervals_per_seq=self.intervals_per_seq,
            genome_subset_patterns=self.genome_subsets_patterns)
        # ["T"]["autosomes"]["unmapped_coords"]
        dataset = pd.DataFrame({"x": [], "y": [], "hue": []})
        for db, d in m.items():
            for subset, dd in d.items():
                dd = dd["unmapped_coords"].value_counts().sort_index()
                new_rows = pd.DataFrame({
                    "x": dd.index, "y": dd.values,
                    "hue": [str(db + " - " + subset)] * len(dd.index), })
                dataset = pd.concat([dataset, new_rows])

        sns.kdeplot(data=dataset, x="x", weights="y", hue="hue")
        plt.xscale("log")
        #> g.set(
        #> g.set_xlim(left=0)
        plt.show()

        return None

def interpret_cigar_string(cigar_string: str):
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


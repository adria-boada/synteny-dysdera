#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# paranome_families_representation.py
#
# 24 de març 2023  <adria@molevol-OptiPlex-9020>

"""
Converts a couple of `file.gff3.gene_filtered` and `file.mcl` into a more
adequate file for representation with CIRCOS.

The newly created file is a TSV with 5 columns: gene ID, beginning and ending
coordinates, in which chromosome is it located, and of which paralogous family
is it (numbered arbitrarilly from 1 to N).

This format allows sorting the file by paralogous family, facilitating the
representation efforts.

"""

import sys, tabulate

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class Paranome:

    def __init__(self, file1, file2):
        """
        Define which two files will create the Paranome class.
        Assign each file to the pertinent variable by detecting the *.gff3 extension
        """
        print("WARNING: remove all blank lines from the TSV.MCL file "+
              "with, for example, `sed '/^$/d' bad.tsv.mcl > good.tsv.mcl`")
        if ".gff3" in file1:
            self.file_gff = file1
            self.file_paralogy = file2
        elif ".gff3" in file2:
            self.file_gff = file2
            self.file_paralogy = file1
        else:
            sys.exit("ERROR: no file with *.gff3 "+
                  "extension has been provided")
        # verifica que el camí fins als fitxers és vàlid
        try:
            with open(self.file_gff) as fg, open(self.file_paralogy) as fp:
                pass
        # si no és vàlid, exit()
        except:
            sys.exit('ERROR: The path to the provided files appears to be '+
                     'incorrect')
        self.parsed_genes = self.parsing()

    # define all fields of a given GFF3 line thanks to a nested class
    # (requires to be called as 'self.Classname' inside the class)
    class GFF3_line:
        def __init__(self, line):
            (self.sequid,           # scaffold or chromosome ID or name
             self.source,           # source of annotation
             self.feature_type,     # type of feature
             self.start,
             self.end,
             self.score,
             self.strand,
             self.phase,
             self.attributes,
             ) = line.split()

            # int-ize numeric values:
            self.start = int(self.start)
            self.end = int(self.end)

    def parsing(self):
        """
        For each family, find all of their paralogous genes inside the gff3 and
        pull out information about them.
        """
        family = 0
        stored_genes = []
        with open(self.file_paralogy) as fp:
            for line_par in fp:
                family += 1
                print(f"Status: processing family number {family}")
                for gene_par in line_par.split():
                    with open(self.file_gff) as fg:
                        for line_gff in fg:
                            # "gene_gff3" stores the gene ID from a given gff3 line.
                            gene_gff = line_gff.split()[8].split('=')[1].split(';')[0]
                            # if both genes are concordant and feature type is "gene"
                            if str(gene_par) == str(gene_gff):
                                # create a "Paranome" object from the gff3 line:
                                par = self.GFF3_line(line_gff)
                                # store values of interest in a list
                                # each object from the list is another list of fields
                                stored_genes += [[par.sequid,   # chromosome/scaff
                                                  gene_par,     # gene ID
                                                  family,       # arbitrary family number
                                                  par.start,    # gff3 start
                                                  par.end,      # gff3 end
                                                  par.strand,   # gff3 strand
                                                  par.end - par.start # gene length
                                            ]]
                                break
                        if str(gene_par) != str(gene_gff):
                            print(f"WARNING: the paralogous gene {gene_par} "+
                                  "has not found its coordinates in the GFF3 "+
                                  "file")

        # return the list of fields of interest for each gene
        # once the whole file with paralogous families has been parsed.
        return stored_genes

    def families_dict(self):
        """
        Get a dictionary with arbitrary family numbers as keys and lists of
        genes in the family as items.
        """
        fam = 0
        dresult = {}
        for gene in self.parsed_genes:
            if fam == gene[2]:
                dresult[fam] += [list(gene)]
            else:
                fam = gene[2]
                dresult[fam] = [list(gene)]

        return dresult

    def linkage(self):
        """
        Parse the items of the dictionary generated in the previous function.
        Decipher linkage inside and between chromosomes. Return a file with the
        linkages of interest.
        """
        # Define intrachromosomal ranges to detect tandem paralogs
        # (e.g. less than 10 kbs, less than 100 kbs, less than 1 Mbs, beyond 1 Mbs)
        intrachr_len_filter = [10**4, 10**5, 10**6]

        # List the unique scaffolds/chromosomes which will be analyzed
        unique_chr = []
        for gene in self.parsed_genes:
            if gene[0] not in unique_chr:
                unique_chr += [gene[0]]
        print("Status: List of unique chr queried for analyses:")
        [print(f"Status:  + {x}") for x in unique_chr]
        print() # newline to separate previous list

        # List the unique possible connections between chromosomes
        # (chr1-chr1->1Mb ; chr1-chr2 ; etc.)
        # Assign to each connection a dictionary with the variables of interest
        # ('l' for sum of both paralogs' bp, 'n' for number of connections)
        connection = {}
        # Intrachromosomal connections:
        for c in unique_chr:
            for s in intrachr_len_filter:
                connection[f"{c}--dup<{s}"] = {'n':0, 'l':0}
            connection[f"{c}--dup>={s}"] = {'n':0, 'l':0}
        # Interchromosomal connections:
        reduce_chr = unique_chr[:]
        while len(reduce_chr) != 1:
            first_chr = reduce_chr.pop(0)
            for c in reduce_chr:
                connection[f"{first_chr}--{c}"] = {'n': 0, 'l': 0}
            print("Status: finished interchromosomal connections for", first_chr)
        #print(connection)#DEBUG

        # Find the information of linked paralogs and fill the previous list of
        # unique connections.
        total_connection_num = 0
        total_connection_len = 0
        for par_gen_fam in self.families_dict().values():
            print(f"Status: processing family {[x[1] for x in par_gen_fam]}")
            # If there is one gene remaining, break
            while len(par_gen_fam) != 1:
                # Take the first gene in the family
                first_gene = par_gen_fam.pop(0)
                print(f"Status: comparing gene {first_gene[1]} with the rest of "+
                      f"the family {[x[1] for x in par_gen_fam]}")
                # Compute connections with the remaining genes
                for gp in par_gen_fam:
                    # if different sequid, interchromosomal connection:
                    if first_gene[0] != gp[0]:
                        try:
                            connection[f"{first_gene[0]}--{gp[0]}"]['n'] += 1
                            connection[f"{first_gene[0]}--{gp[0]}"]['l'] += first_gene[6] + gp[6]
                            total_connection_num += 1
                            total_connection_len += first_gene[6] + gp[6]
                        except:
                            connection[f"{gp[0]}--{first_gene[0]}"]['n'] += 1
                            connection[f"{gp[0]}--{first_gene[0]}"]['l'] += first_gene[6] + gp[6]
                            total_connection_num += 1
                            total_connection_len += first_gene[6] + gp[6]
                    # if the same sequid, intrachromosomal connection:
                    else:
                        if any([
                        first_gene[3] <= gp[3] <= first_gene[4],
                        gp[3] <= first_gene[3] <= gp[4]
                        ]):
                            # if they overlap, distance is zero
                            distance = 0
                        else:
                            # otherwise...
                            distance = int(
                                max(first_gene[3], gp[3]) -
                                min(first_gene[4], gp[4])
                            )
                        # check in which range does the computed distance fall:
                        for s in intrachr_len_filter:
                            if distance < s:
                                connection[f"{first_gene[0]}--dup<{s}"]['n'] += 1
                                connection[f"{first_gene[0]}--dup<{s}"]['l'] += first_gene[6] + gp[6]
                                total_connection_num += 1
                                total_connection_len += first_gene[6] + gp[6]
                                break
                        if distance >= s:
                            connection[f"{first_gene[0]}--dup>={s}"]['n'] += 1
                            connection[f"{first_gene[0]}--dup>={s}"]['l'] += first_gene[6] + gp[6]
                            total_connection_num += 1
                            total_connection_len += first_gene[6] + gp[6]

        print(f"Status: found as many as {total_connection_num} connections, "+
              f"with total len {total_connection_len}")
        return connection

    def basic_parfam_properties(self):
        """
        Print a table with the number of paralogous genes,
        number of families, proportion of genes/bp per family...
        """
        genes_in_families = []
        families = 0
        with open(self.file_paralogy) as fp:
            for line in fp:
                # recupera el nombre de paràlegs dins la família d'aquesta línia
                genes_in_families += [len(line.split())]
                # recupera el nombre de línies al fitxer (núm. famílies + orfes)
                families += 1

        # recompta famílies amb "membres >= 2" (sols famílies de paralogia)
        paralogous_families = 0
        with open(self.file_paralogy) as fp:
            for line in fp:
                if len(line.split()) == 1:
                    break
                else:
                    paralogous_families += 1

        # recompta el nombre total de gens paràlegs
        # (suma els membres de cada família)
        paralogous_genes = 0
        for i in range(2, max(genes_in_families)+1):
            paralogous_genes += i*genes_in_families.count(i)

        return {
            "Màxim de membres en una sola família":
            max(genes_in_families),
            "Núm. de gens paràlegs":
            paralogous_genes,
            "Núm. de gens paràlegs + no paràlegs":
            paralogous_genes + genes_in_families.count(1),
            "Núm. de famílies paralogues + gens orfes":
            families,
            "Núm. de famílies paralogues":
            paralogous_families,
        }

    """
    The following two table generators parse the files independently.
    If both are called, the parsing will be done twice (slow)
    """
    def tsv_tbl(self):
        # print the list parsed previously in parsing()
        for gene_fields in self.parsed_genes:
            # imprimeix fins al penúltim camp:
            for field in gene_fields[0:-1]:
                print(field, end='\t')
            # imprimeix l'últim camp amb nova línia:
            print(gene_fields[-1])

    def pipe_tbl(self):
        # print a formatted list parsed previously in parsing()
        capçalera = ["Sequence ID", "Gene ID", "Arbitrary paralogous family",
                  "Start", "End", "Strand"]
        print(tabulate.tabulate(self.parsed_genes, tablefmt='pipe',
                            headers=capçalera))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Si es crida com a script:
if __name__ == '__main__':

    ### PARSE ARGUMENTS ###
    import argparse

    parser = argparse.ArgumentParser(
        description="Provide two files: (a) GFF3 annotations with *.gff3 "+
        "extension and (b) paranome file (in any order)")
    # file-name: positional arg.
    parser.add_argument('file1', type=str,)
    parser.add_argument('file2', type=str,)
    # integer argument
#    parser.add_argument('-a', '--numero_a', type=int, help='Paràmetre "a"')
    # choices argument
#    parser.add_argument('-o', '--operacio', 
#            type=str, choices=['suma', 'resta', 'multiplicacio'],
#            default='suma', required=False,
#            help='')

    args = parser.parse_args()
    # call a value: args.operacio or args.filename.

    # parse both the GFF3 and TSV.MCL thanks to the "Paranome" class

    # file1 and file2 have to be GFF3 and MCL (in any order):
    paranome = Paranome(args.file1, args.file2)
    # create a dictionary where each key-value pair are an arbitrary family of
    # paralogy and a list of their genes
#    print("LIST OF GENES FOUND BY PARSING()")#DEBUG
#    [print("+ ", gene) for gene in paranome.parsed_genes]#DEBUG
#    print(paranome.linkage())
    for key, val in paranome.linkage().items():
        print(key, val)
    # print table 'pipe' formatted:
    paranome.pipe_tbl()
    # print table 'tsv' formatted:
#    paranome.tsv_tbl()
    # print basic stats from the MCL file:
#    for key, val in paranome.basic_parfam_properties().items():
#        print(key, ':', val)


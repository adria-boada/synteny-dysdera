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
        Accepts two files and defines which one contains the *.gff3 extension
        """
        if ".gff3" in file1:
            self.file_gff = file1
            self.file_paralogy = file2
        elif ".gff3" in file2:
            self.file_gff = file2
            self.file_paralogy = file1
        else:
            print("Error in creating the Paranome class: no file with gff3",
                  "extension found")
        try:
            with open(self.file_gff) as fg, open(self.file_paralogy) as fp:
                pass
        except:
            sys.exit('The files have not been correctly provided')

    # define all fields of a given GFF3 line thanks to a nested class
    # (requires to be called as 'self.Classname' inside the class)
    class GFF3_line(object):
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
                            if str(gene_par) == str(gene_gff) and line_gff.split()[2] == "mRNA":
                                # create a "Paranome" object from the gff3 line:
                                paralog = self.GFF3_line(line_gff)
                                # store values of interest in a list
                                # each object from the list is another list of fields
                                stored_genes += [[paralog.sequid,   # chromosome/scaff
                                                  gene_par,        # gene ID
                                                  family,           # arbitrary family number
                                                  paralog.start,    # gff3 start
                                                  paralog.end,      # gff3 end
                                                  paralog.strand,   # gff3 strand
                                            ]]
                                break

        # return the list of fields of interest for each gene
        # once the whole file with paralogous families has been parsed.
        return stored_genes

    def basic_parfam_properties(self):
        # print a table with number of paralogous genes, etc.

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Si es crida com a script:
if __name__ == '__main__':

    ### PARSE ARGUMENTS ###
    import argparse

    parser = argparse.ArgumentParser(
        description="Provide two files: (1) GFF3 annotations with *.gff3 extension and (2) paranome file (in any order)")
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

#    for key, val in Paranome.basic_properties(args.file_parg).items():
#        print(key, ':', val)
    # parse both the GFF3 and TSV.MCL with the "Paranome" class
    paranome = Paranome(args.file1, args.file2)
    formatted_paranome_list = paranome.parsing()
    # print the formatted list with "tabulate" module
    capçalera = ["Sequence ID", "Gene ID", "Arbitrary paralogous family",
                  "Start", "End", "Strand"]
    print(tabulate.tabulate(formatted_paranome_list, tablefmt='pipe',
                            headers=capçalera))
    for gene_fields in formatted_paranome_list:
        for field in gene_fields:
            print(field, end='\t')
        print()



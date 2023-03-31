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

class Paranome(object):
    """ Parse fields given a GFF3 and MCL formatted files. Create an instance of
    the "Paranome" class from a single line of the GFF3 file, and store it in a
    list. Return this list with the parsing() function.
    """
    def __init__(self, gff3_line):
        # define all fields of the given GFF3 line
        (self.sequid,           # scaffold or chromosome ID or name
         self.source,           # source of annotation
         self.feature_type,     # type of feature
         self.start,
         self.end,
         self.score,
         self.strand,
         self.phase,
         self.attributes,
         ) = gff3_line.split()

        # int-ize numeric values:
        self.start = int(self.start)
        self.end = int(self.end)

    def parsing(
        file_gff3: "Path to a GFF3 formatted file which has been filtered by 'gene' features",
        file_parg: "Path to an MCL-like file in which each line is a family of paralogous genes",
                ):
        # for each family, find all of their paralogous genes inside
        # the gff3 and pull out information about them.
        family = 0
        stored_genes = []
        with open(file_parg) as fparg:
            for line_parg in fparg:
                family += 1
                print(f"Status: processing family number {family}")
                for gene_parg in line_parg.split():
                    with open(file_gff3) as fgff3:
                        for line_gff3 in fgff3:
                            # "gene_gff3" stores the gene ID from a given gff3 line.
                            gene_gff3 = line_gff3.split()[8].split('=')[1].split(';')[0]
                            # if both genes are concordant and feature type is "gene"
                            if str(gene_parg) == str(gene_gff3) and line_gff3.split()[2] == "mRNA":
                                # create a "Paranome" object from the gff3 line:
                                paralog = Paranome(line_gff3)
                                # store values of interest in a list
                                # each object from the list is another list of fields
                                stored_genes += [[paralog.sequid,   # chromosome/scaff
                                                  gene_parg,        # gene ID
                                                  family,           # arbitrary family number
                                                  paralog.start,    # gff3 start
                                                  paralog.end,      # gff3 end
                                                  paralog.strand,   # gff3 strand
                                            ]]
                                break

        # return the list of fields of interest for each gene
        # once the whole file with paralogous families has been parsed.
        return stored_genes

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Si es crida com a script:
if __name__ == '__main__':

    ### PARSE ARGUMENTS ###
    import argparse

    parser = argparse.ArgumentParser(description='')
    # file-name: positional arg.
    parser.add_argument('file_gff3', type=str,
                        help='Path to GFF3 annotations file-name')
    parser.add_argument('file_parg', type=str,
                        help='Path to paranome file-name')
    # integer argument
#    parser.add_argument('-a', '--numero_a', type=int, help='Paràmetre "a"')
    # choices argument
#    parser.add_argument('-o', '--operacio', 
#            type=str, choices=['suma', 'resta', 'multiplicacio'],
#            default='suma', required=False,
#            help='')

    args = parser.parse_args()
    # call a value: args.operacio or args.filename.

    # parse both the GFF3 and TSV.MCL with the "Paranome" class
    formatted_paranome_list = Paranome.parsing(args.file_gff3, args.file_parg)
    # print the formatted list with "tabulate" module
    capçalera = ["Sequence ID", "Gene ID", "Arbitrary paralogous family",
                  "Start", "End", "Strand"]
    print(tabulate.tabulate(formatted_paranome_list, tablefmt='pipe',
                            headers=capçalera))
#    for gene_fields in formatted_paranome_list:
#        for field in gene_fields:
#            print(field, end='\t')
#        print()



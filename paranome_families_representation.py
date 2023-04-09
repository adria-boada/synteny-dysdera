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

import sys, tabulate, pandas as pd
import matplotlib.pyplot as plt, numpy as np

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class Paranome:

    def __init__(self, file1, file2):
        """
        Define which two files will create the Paranome class.
        Assign each file to the pertinent variable by detecting the *.gff3 extension
        """
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

        # finally, parse and withdraw paralogous genes info...
        self.parsed_genes = self.parsing()
        self.df_tracks, self.df_links = self.linkage()

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
        # read the amount of families in the file (=number of lines)
        with open(self.file_paralogy) as fp:
            for i, l in enumerate(fp):
                if len(l) == 1:
                    print("WARNING: remove all blank lines from the TSV.MCL file "+
                    "with, for example, `sed '/^$/d' bad.tsv.mcl > good.tsv.mcl`")
            tot_fam = i+1
        # start searching information through the GFF3
        family = 0
        stored_genes = []
        with open(self.file_paralogy) as fp:
            for line_par in fp:
                family += 1
                print(f"Status: Withdrawing gene information from the GFF3 "+
                      f"for family {family}/{tot_fam}")
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
            # gene[0] is sequid
            if gene[0] not in unique_chr:
                unique_chr += [gene[0]]
        # Sort sequids alphabetically
        unique_chr.sort()
        print("Status: List of unique chr queried for analyses:")
        [print(f"Status:  + {x}") for x in unique_chr]

        # Create pandas DataFrame storing the results
        list_colnames = [
            'Chromosome length',
            "Sum of paralogs' length",
            'Number of paralogs',
            'Number of interchr. conns.',
        ]
        for s in intrachr_len_filter:
            list_colnames += [f'Dups < {s} bp']
        list_colnames += [f'Dups >= {s} bp']
        # Connections for individual sequids
        blocks = pd.DataFrame()
        for coln in list_colnames:
            blocks[coln] = 0
        for crm in unique_chr:
            blocks.loc[crm] = [0]*len(list_colnames)
        # Connections between sequids
        intchr_links = pd.DataFrame()
        for crm in unique_chr:
            intchr_links[crm] = 0
        for crm in unique_chr:
            intchr_links.loc[crm] = [0]*len(unique_chr)
        print("Status: Created the blank pandas DataFrames")

        # Fill the information of the pandas DataFrames:
        for gene in self.parsed_genes:
            # Sum of paralogs' length for each sequid:
            blocks.loc[gene[0], "Sum of paralogs' length"] += gene[6]
            # Number of paralogs for each sequid:
            blocks.loc[gene[0], "Number of paralogs"] += 1

        # Find the information of linked paralogs and keep filling the previous
        # list with unique connections between pairs of sequids.
        tot = len(self.families_dict().keys())
        for fam_key, par_gen_fam in self.families_dict().items():
            if len(par_gen_fam) == 1:
                print(f"Status: Skipping pairwise connections between genes "+
                      f"from the family {fam_key}/{tot} (only one member)")
            else:
                print(f"Status: Calculating pairwise connections between genes "+
                      f"from family {fam_key}/{tot}")
            # If there is one gene remaining, break
            while len(par_gen_fam) != 1:
                # Take the first gene in the family
                first_gene = par_gen_fam.pop(0)
#DBUG                print(f"Status: comparing gene {first_gene[1]} with the rest of "+
#DBUG                      f"the family ({len(par_gen_fam)-1} connections left)")
                # Compute connections with the remaining genes
                for gp in par_gen_fam:
                    # if different sequid, interchromosomal connection:
                    if first_gene[0] != gp[0]:
                        s = [first_gene[0], gp[0]]
                        s.sort()
                        # add a connection in the matrix between chr[0] and
                        # chr[1] and add one connection for each to the 'blocks'
                        # pd.DataFrame()
                        intchr_links.loc[s[1], s[0]] += 1
                        intchr_links.loc[s[0], s[1]] += 1
                        blocks.loc[s[0], 'Number of interchr. conns.'] += 1
                        blocks.loc[s[1], 'Number of interchr. conns.'] += 1
                        continue
                    # else, they are in the same chromosome (intrachr):
                    # add in the diagonal of the matrix (connection chr1 to chr1)
                    intchr_links.loc[first_gene[0], gp[0]] += 1
                    # find the distance of genes in the chromosome
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
                    # check in which range does the computed distance fall,
                    # and add one connection to that distance range.
                    for s in intrachr_len_filter:
                        if distance < s:
                            blocks.loc[first_gene[0], f'Dups < {s} bp'] += 1
                            break
                    if distance >= s:
                        blocks.loc[first_gene[0], f'Dups >= {s} bp'] += 1

        # Compute the total of each column:
        for coln in list_colnames:
            blocks.loc['*Total*', coln] = blocks[coln].sum()

        return (blocks, intchr_links)

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

    def heatmap(self, filename: "Output filename of the heatmap" ='linkage_heatmap.png'):
        """
        Create a heatmap from the pairwise connections between chr
        """
        df = self.df_links
        plt.pcolormesh(df)
        plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
        plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns)
        plt.savefig(filename)

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
                  "Start", "End", "Strand", "Gene length"]
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

    args = parser.parse_args()
    # call a value: args.operacio or args.filename.

    # parse both the GFF3 and TSV.MCL thanks to the "Paranome" class
    # file1 and file2 have to be GFF3 and MCL (vars can go in any order):
    paranome = Paranome(args.file1, args.file2)
    # print retrieved pandas Dataframes
    print(paranome.df_tracks)
    print(paranome.df_links)
    paranome.heatmap()
    # print table 'pipe' formatted:
#    paranome.pipe_tbl()
    # print table 'tsv' formatted:
#    paranome.tsv_tbl()
    # print basic stats from the MCL file:
    for key, val in paranome.basic_parfam_properties().items():
        print(key, ':', val)


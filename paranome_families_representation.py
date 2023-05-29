#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# paranome_families_representation.py
#
# 24 de març 2023  <adria@molevol-OptiPlex-9020>

"""
Converts a couple of `file_gene_filtered.gff3` and `file.mcl` into a more
adequate file for representation with CIRCOS. Leaving only geneIDs should
improve performance.

The newly created file is a TSV with 5 columns: gene ID, beginning and ending
coordinates, in which chromosome is it located, and of which paralogous family
is it (numbered arbitrarilly from 1 to N).

This format allows sorting the file by paralogous family, facilitating the
representation efforts.

Es podria filtrar els sequids amb menys de N paralegs i evitar renombrar.
"""

import sys, tabulate
import pandas as pd
import matplotlib.pyplot as plt, numpy as np
import math
import natsort

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class Paranome:

    def __init__(self, file1, file2,
                 idx_file: 'Optional index file with chr. lengths'=0):
        """
        Define which two files will create the Paranome class.
        Assign each file to the pertinent variable by detecting the *.gff3
        extension.

        One file is a GFF3 with a description of all annotated genes.
        The second, an MCL file which delimits groups of paralogs.
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
        # verify that the provided paths are valid:
        try:
            with open(self.file_gff) as fg, open(self.file_paralogy) as fp:
                pass
        # if paths are not valid, exit the script.
        except:
            sys.exit('ERROR: The path to the provided files appears to be '+
                     'incorrect')
        if idx_file:
            try:
                with open(idx_file) as fx:
                    pass
            except:
                sys.exit('ERROR: The path to the provided IDX file appears to be '+
                     'incorrect')

        # parse and withdraw paralogous genes info...
        self.df_parsed_genes = pd.DataFrame(self.parsing(),
                                  columns=['sequid', 'geneID',
                                         'arbitrary_fam_number',
                                         'start', 'end', 'strand',
                                         'gene_length'])
        # afterwards, find presence of families in sequids...
        self.df_family_presence = self.fam_presence()
        # compute matrix of linkage strengths between sequids
        # (keep in mind that there are multiple func. with diff. weighs)
        self.df_linkage = self.linkage_v3()
        # remove rows and cols with only zeroes from the linkage dataframe
        self.df_clean_linkage =(
            self.df_linkage.loc[self.df_linkage.any(axis=1),
                                self.df_linkage.any(axis=0)])
        # reorder index and colnames with `humansorted` from `natsort`:
        # Scaffold_2 should be placed between Scaffold_1 and Scaffold_11
        self.df_clean_linkage =(
            self.df_clean_linkage.reindex(columns=natsort.humansorted(self.df_clean_linkage.columns)
                                ).reindex(index=natsort.humansorted(self.df_clean_linkage.index)))

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
                # print status every so often (remainder of division by 25 == 1)
                if int(family) % 25 == 1:
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

    def fam_presence(self):
        """
        Parse the families of paralogs. Store the presence of particular
        families in each sequid.

        Example pandas df output:
                chrA    chrB
        fam_1   4       0
        fam_2   2       1
        """
        # Define intrachromosomal ranges to detect tandem paralogs
        # It'd be interesting to report at multiple levels (compare with
        # i-adhore tandem gene tags)
        #
        # The paralogs df is stored in the var 'self.df_parsed_genes'
        # Shorten the calls to this specific df:
        df = self.df_parsed_genes
        # Init. a matrix of zeroes. Each row, a family of paralogy. Each column,
        # the amount of paralogs from the index family present in column
        # chromosome.
        return_families_presence = pd.DataFrame(0,
                                            columns=df['sequid'].unique(),
                                            index=df['arbitrary_fam_number'].unique())
        counter_end = self.df_parsed_genes['arbitrary_fam_number'].max()
        # Compute the amount of paralogs for given families in every sequid.
        # For all families of paralogy in the dataframe:
        for family in df['arbitrary_fam_number'].unique():
            if int(family) % 25 == 1:
                print('Status: searching for presence of family '+
                     f'{family}/{counter_end}')
            # For all sequids containing members of the selected family:
            for crm in df[df['arbitrary_fam_number']==family]['sequid'].unique():
                subset = df.query(f'arbitrary_fam_number == {family} & '+
                                  f'sequid == "{crm}"')
                return_families_presence.loc[family, crm] = len(subset)
        # maybe it would be better to return a dataframe with the amount of
        # paralog genes present in each sequid for each family
        return return_families_presence

    def linkage_v3(self,
                   mode:"Which function to use when computing sequids relatedness"=1):
        """
        Given a dataframe of families presence across sequids, calculate
        relatedness/linkage in and between sequids
        """
        if mode == 1:
            funct = self.analytical_sharedness_pairwise
        elif mode == 2:
            pass
        # The paralogs df is stored in the var 'self.df_parsed_genes'
        # Shorten the calls to this specific df:
        df = self.df_parsed_genes
        # Init. a square matrix of zeroes. Each colname and rowname are unique
        # chromosome names, representing strength of connections between
        # every pair of chromosomes (even with themselves; diagonal)
        return_paralogy_sharedness = pd.DataFrame(0,
                                         index=df['sequid'].unique(),
                                         columns=df['sequid'].unique())
        for row in self.df_family_presence.iterrows():
            r = row[1][row[1] >0]
            # iterate until all the values from the series have been popped()
            while len(r) >0:
                # pop the first value of family presence in `r` series
                # at the same time, store its index in another var (crm.)
                curr_crm, curr_val = (r.index[0], r.pop(r.index[0]))
                # compare popped() value to itself
                return_paralogy_sharedness.loc[curr_crm, curr_crm] +=(
                    funct(curr_val) )
                # compare popped() value to the remaining values in `r` series
                # both `curr_crm` and `i` are sequids
                for i in r.index:
                    return_paralogy_sharedness.loc[curr_crm, i] +=(
                        funct(curr_val, r[i]) )
                    # in order to avoid triangular dataframes (dataframe with
                    # zeroes) add the computed value above and below the
                    # diagonal:
                    return_paralogy_sharedness.loc[i, curr_crm] +=(
                        funct(curr_val, r[i]) )
        return return_paralogy_sharedness

    def analytical_sharedness_pairwise(self, val1, val2=0):
        """
        Compute sharedness of chromosome by calculating the number of all possible
        pairwise connections between paralogs of this family.
        """
        if val1==1 and not val2:
            return 0
        elif val1>1 and not val2:
            return (val1*(val1-1))/2
        else: # val1 and val2
            return math.ceil( (val1*val2)/2 )

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

    def linkage(self,
                crm_idx_len: 'Index file'):
        """
        Parse the items of the dictionary generated in the previous function.
        Decipher linkage inside and between chromosomes. Return a file with the
        linkages of interest.
        """
        # Define intrachromosomal ranges to detect tandem paralogs
        # (e.g. less than 10 kbs, less than 100 kbs, less than 1 Mbs, beyond 1 Mbs)
        # from now on, intra=outer
        outer_chr_len_filter = [10**4, 10**5, 10**6]

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
            'chr_len',
            'num_paralogs',
            'sum_paralogs_len',
            'num_outer_conn',
            'num_interior_conn',
        ]
        for s in outer_chr_len_filter:
            # Outer connections less than {s} -> putatively tandem dups.
            list_colnames += [f'interior_conn_lt_{s}']
        # Outer connections greater (or equal) than {s} -> probably not tandem
        list_colnames += [f'interior_conn_gt_{s}']
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
            blocks.loc[gene[0], "sum_paralogs_len"] += gene[6]
            # Number of paralogs for each sequid:
            blocks.loc[gene[0], "num_paralogs"] += 1

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
                        blocks.loc[s[0], 'num_outer_conn'] += 1
                        blocks.loc[s[1], 'num_outer_conn'] += 1
                        continue
                    # else, they are in the same chromosome (outer_chr):
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
                    for s in outer_chr_len_filter:
                        if distance < s:
                            blocks.loc[first_gene[0], f'outer_conn_lt_{s}'] += 1
                            blocks.loc[first_gene[0], 'num_interior_conn'] += 1
                            break
                    if distance >= s:
                        blocks.loc[first_gene[0], f'outer_conn_gt_{s}'] += 1
                        blocks.loc[first_gene[0], 'num_interior_conn'] += 1

        # If an idx file has been given, specify length of each chr in the
        # DataFrame...
        if crm_idx_len:
            with open(crm_idx_len) as idx:
                for line in idx:
                    crm, length = line.split()
                    if crm in unique_chr:
                        blocks.loc[crm, 'chr_len'] += float(length)

        # Translate amount of connections to percentages of total
        # Then, the CIRCOS track will be divided by these percentages:
        perc_conn = blocks.loc[:, 'num_outer_conn':].drop('num_interior_conn', axis=1)
        circos_blocks = pd.DataFrame()
        circos_blocks['chr_len'] = blocks['chr_len']
        tot_con = blocks['num_outer_conn'] + blocks['num_interior_conn']
        for col in perc_conn:
            circos_blocks[col] = blocks['chr_len']*(perc_conn[col]/tot_con)

        # Compute the total of each column:
        tot = pd.DataFrame()
        for coln in blocks:
            tot.loc['*Total*', coln] = blocks[coln].sum()
        for coln in circos_blocks:
            tot.loc['*Total*', coln] = circos_blocks[coln].sum()

        return (blocks, intchr_links, circos_blocks, tot)

    def linkage_v2(self,
                crm_idx_len: 'Index file'=0):
        """
        Second attempt at representing intrachromosomal and interchromosomal
        (from now on 'outer' and 'interior' connections) of each family.
        Allocate how many genes could originate from tandem dups and how many
        from recent WGD (not a rearrangement?)
        """
        # Define intrachromosomal ranges to detect tandem paralogs
        # (e.g. less than 10 kbs, less than 100 kbs, less than 1 Mbs, beyond 1 Mbs)
        # from now on, intra=outer
        tandem_genes_dist_filter = [10**4, 10**5, 10**6]

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

        df_colnames = [
            'chr_len',
            'sum_paralogs_len',
            'sum_paralogs_dist',
            'trials_pargs_dist',
            'avg_paralogs_dist',
            'num_paralogs',
            'num_paralogs_outer',
        ]
        for s in tandem_genes_dist_filter:
            # Outer connections less than {s} -> putatively tandem dups.
           df_colnames += [f'parg_lt_{s}']
        # Outer connections greater (or equal) than {s} -> probably not tandem
        df_colnames += [f'parg_gt_{s}']
        # Overall stats for each chromosome
        blocks = pd.DataFrame(columns=df_colnames, index=unique_chr)
        blocks.loc[:]=0
        # Shared families between chromosomes
        # (in how many chromosomes is this family present,
        # and in which gene quantity?)
        intchr_links = pd.DataFrame(columns=unique_chr, index=unique_chr)
        intchr_links.loc[:]=0
        print("Status: Created the blank pandas DataFrames")

        # Fill the information of the pandas DataFrames:
        for gene in self.parsed_genes:
            # Sum of paralogs' length for each sequid:
            blocks.loc[gene[0], "sum_paralogs_len"] += gene[6]
            # Number of paralogs for each sequid:
            blocks.loc[gene[0], "num_paralogs"] += 1
        # If an idx file has been given, specify length of each chr in the
        # DataFrame...
        if crm_idx_len:
            with open(crm_idx_len) as idx:
                for line in idx:
                    crm, length = line.split()
                    if crm in unique_chr: #DEBUG
                        blocks.loc[crm, 'chr_len'] += float(length)

        # Compute the 'outer' column
        # (how much are genes of the same family split between chromosomes)
        tot = len(self.families_dict().keys())
        for fam_key, par_gen_fam in self.families_dict().items():
            if len(par_gen_fam) == 1:
                print(f"Status: Skipping outer connections and distances "+
                      f"for family {fam_key}/{tot} (only one member)")
            else:
                print(f"Status: Calculating intrachromosomal distances and "+
                      f"the distribution of paralogs across chromosomes "+
                      f"for family {fam_key}/{tot}")
                outer = pd.Series(index=unique_chr, dtype=float)
                outer.loc[:]=0
            # If there is one gene remaining, break
            while len(par_gen_fam) != 1:
                # Take the first gene in the family
                gene = par_gen_fam.pop(0)
#DBUG                print(f"Status: comparing gene {first_gene[1]} with the rest of "+
#DBUG                      f"the family ({len(par_gen_fam)-1} connections left)")
                # Compute distance to the other genes if they're in the same sequid
                for gp in par_gen_fam:
                    if gene[0] == gp[0]: # zero index is sequid
                        dst = self.gene_distance(gene, gp)
                        blocks.loc[gene[0], 'sum_paralogs_dist']+=dst
                        blocks.loc[gene[0], 'trials_pargs_dist']+=1
                        # check in which range does the computed distance fall,
                        # and add one connection to that distance range.
                        for s in tandem_genes_dist_filter:
                            if dst < s:
                                blocks.loc[gene[0], f'parg_lt_{s}'] += 1
                                blocks.loc[gene[0], 'num_paralogs'] += 1
                                break
                        if dst >= s:
                            blocks.loc[gene[0], f'parg_gt_{s}'] += 1
                            blocks.loc[gene[0], 'num_paralogs'] += 1
                # keep track of in which chr. can we find this family's paralogs
                outer.loc[gene[0]]+=1
            outer.loc[par_gen_fam.pop(0)[0]]+=1
            # take only presence of family, not the amount of genes, dividing by
            # itself (+1 present, +0 absent family in 'column' chromosome)
            f_outer = outer[outer!=0]/outer[outer!=0]

            # Compute how many outer genes from this family each chr has...
            for crm in unique_chr:
                # suma tots els gens a `crm` menys la columna propia:
                outsiders = outer.sum()-outer.loc[crm]
                blocks.loc[crm, 'num_paralogs_outer']+=outsiders
            # Compute how strong are the links between chr...
            intchr_links.loc[f_outer.index, f_outer.index]+=1

        # Compute average distances between paralogs of the same chr:
        blocks['avg_paralogs_dist']=(
            blocks.loc[blocks['sum_paralogs_dist']!=0, 'sum_paralogs_dist'] /
            blocks.loc[blocks['sum_paralogs_dist']!=0, 'trials_pargs_dist'])

        return (blocks, intchr_links)

    def gene_distance(self, stored_gene1, stored_gene2):
        """
        Compute distance between given pair of genes parsed by the above
        function
        """
        # if they overlap, distance is zero
        if any([
        stored_gene1[3] <= stored_gene2[3] <= stored_gene1[4],
        stored_gene2[3] <= stored_gene1[3] <= stored_gene2[4]
        ]):
            distance = 0
        # otherwise, the difference between max beginning and min end
        else:
            distance = int(
                max(stored_gene1[3], stored_gene2[3]) -
                min(stored_gene1[4], stored_gene2[4])
            )
        return distance

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
        # count the total number of paralog genes
        # (sum the paralogs across all families)
        paralogous_genes = 0
        for i in range(2, max(genes_in_families)+1):
            paralogous_genes += i*genes_in_families.count(i)

        # create the dictionary which will be printed and returned by the funct.
        return_measures = {
            "Màxim de membres en una sola família":
            max(genes_in_families),
            "Núm. de gens paràlegs":
            paralogous_genes,
            # orphan genes are defined by not having identifiable paralogs
            "Núm. de gens 'orfes'":
            genes_in_families.count(1),
            "Núm. de gens paràlegs + 'orfes'":
            paralogous_genes + genes_in_families.count(1),
            "Núm. de famílies paralogues":
            paralogous_families,
            "Núm. de 'famílies orfes'":
            families-paralogous_families,
            "Núm. de famílies paralogues + 'orfes'":
            families,}
        # when the funct. is called, print the values to screen:
        for key, val in return_measures.items():
            print(key, ': ', val, sep='')
        return return_measures

    def heatmap(self, filename: "Output filename of the heatmap" ='linkage_heatmap.png'):
        """
        Create a heatmap from the pairwise connections between chr
        """
        df = self.df_clean_linkage
        plt.pcolormesh(df)
        plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
        plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns,
                   # rotate 45 degrees xticks
                   rotation=45)
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
        "extension and (b) paranome file as in MCL's output (files can be in any order)")
    # file-name: positional arg.
    parser.add_argument('file1', type=str,)
    parser.add_argument('file2', type=str,)
    parser.add_argument('-i', '--idx', type=str,
                        help="Chr. index file (with lengths in basepairs)",
                        required=False, default=0)

    args = parser.parse_args()
    # call a value: args.operacio or args.filename.

    # parse both the GFF3 and TSV.MCL thanks to the "Paranome" class
    # file1 and file2 have to be GFF3 and MCL (vars can go in any order):
    paranome = Paranome(args.file1, args.file2)

    # print retrieved pandas Dataframes
    t = "WRITING BOTH DATAFRAMES TO TSV"
    print("-"*len(t), t, sep='\n')
    paranome.df_linkage.to_csv('tmp_raw_linkage.tsv', sep='\t')
    paranome.df_clean_linkage.to_csv('tmp_clean_linkage.tsv', sep='\t')
    # also add a heatmap
    paranome.heatmap()
    # and return paranome stats
    x=paranome.basic_parfam_properties()

    # print table 'pipe' formatted:
#    paranome.pipe_tbl()
    # print table 'tsv' formatted:
#    paranome.tsv_tbl()
    # print basic stats from the MCL file:


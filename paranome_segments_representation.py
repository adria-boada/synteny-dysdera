#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# paranome_segments_representation.py
#
# 05 May 2023  <adria@molevol-OptiPlex-9020>

"""
Working with i-adhore's `segments.txt` output
"""

import sys
import pandas as pd

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class Segments:

    def __init__(self, file1, file2):
        """
        Catches a pair of files (GFF3 and segments.txt)
        in any order.
        Checks that they are valid file paths.
        Creates pandas dataframes and stores as class vars.

        Adds columns of zeros to segments.txt dataframe;
        these initialized columns will later be filled with data
        from the GFF3.
        """
        if ".gff3" in file1 and "segments" in file2:
            path_to_gff3 = file1
            path_to_segments = file2
        elif ".gff3" in file2 and "segments" in file1:
            path_to_gff3 = file2
            path_to_segments = file1
        else:
            sys.exit('ERROR: the path to the provided files '+
                     'appears to be incorrect')
        # verify that provided paths are valid:
        try:
            with open(path_to_segments) as f1, open(path_to_gff3) as f2:
                pass
        # if paths are not valid, exit the script with exception
        except:
            sys.exit('ERROR: the path to the provided files '+
                     'appears to be incorrect')
        # create pandas dataframes from inputted files
        df = pd.read_table(path_to_gff3, header=None)
        df = df.drop(columns=1)
        df.columns = ['sequid', 'feat', 'start', 'end', 'score', 'strand',
                      'phase', 'attributes']
        self.df_input_genes = df
        df = pd.read_table(path_to_segments)
        df = df.drop(columns=['id', 'genome'])
        df.columns = ['multiplicon', 'sequid', 'first', 'last', 'order']
        self.df_input_segments = df
        del df
        print('Status: populating segment dataframe with gene info from GFF3')
        self.df_segments = self.segments_with_gene_info()
        print('Status: linking pairs of segments from the same multiplicon')
        self.df_links, self.parelles = self.paired_segments()
        print('Status: computing ratios of each miscellaneous variables')
        # create a copy from df_links
        df = self.df_links.drop(columns=[
            'gn_from', 'gL_from', 'gD_from', 'gn_to', 'gL_to', 'gD_to',])
        # repeat for all misc. variables (mean gene number, length, distance)
        for i in ['gn', 'gL', 'gD']:
            # compute the ratio of gene length between pair of segments
            df[i] = self.df_links[f'{i}_from']/self.df_links[f'{i}_to']
            # if the ratio is less than 1, substitute by inverse
            df.loc[df[i]<1, i] = 1 / df.loc[df[i]<1, i]
        self.df_links_ratio = df
        del df

    def segments_with_gene_info(self):
        """
        For each row of segments,
            find the row that matches first gene's ID and 'mRNA'
            '   '   '   '   '   '    second gene's   '   '   '  '
            compute the coordinate range for this segment
            store in begin and end

            create a column with a dict of genes found in the
            previously specified coordinate range.
            dict keys are gene IDs, while values are gene coords
            (stored in tuple)
        """
        df_return = self.df_input_segments
        # init columns for later filling of info
        df_return[ ['start', 'end', 'n_genes', 'mean_gene_len', 'mean_gene_dist'] ]=0
        for row in self.df_input_segments.iterrows():
            # select the 'series' object from row tuple
            row = row[1]
            # track progress
            if int(row.name) % 25 == 1:
                print(f'Status: segment {row.name} out of',
                      len(self.df_input_segments))
            # find the genes in gff3 with matching ID in attributes
            gene_first = self.df_input_genes[
                           (self.df_input_genes['attributes'].str.contains(row['first'])) &
                           # we're interested in coordinates (mRNA row)
                           (self.df_input_genes['feat'] == 'mRNA')].squeeze()
            # repeat for 'last' gene in the segment
            gene_last = self.df_input_genes[
                           (self.df_input_genes['attributes'].str.contains(row['last'])) &
                           (self.df_input_genes['feat'] == 'mRNA')].squeeze()
            # the last gene should have higher coords than first gene
            if not (gene_first['start']<gene_first['end']<gene_last['start']<gene_last['end']):
                print('WARNING: segment', row.name,
                      'from multiplicon', row['multiplicon'],
                      'has unsorted coordinates '+
                      '(ending coordinates precede starting coordinates). '+
                      'It will be removed from analysis.')
                continue
            # compute the coordinate range for the current segment
            df_return.loc[row.name, 'start'] = gene_first['start']
            df_return.loc[row.name, 'end'] = gene_last['end']
            # filter the GFF3 dataframe in order to find which genes does the
            # segment contain inside their coordinate range
            df_contains = self.df_input_genes[
                self.df_input_genes['start'].between(gene_first['start'],
                                                     gene_last['end']) |
                self.df_input_genes['end'].between(gene_first['start'],
                                                     gene_last['end'])]
            df_contains = df_contains[
                (df_contains['feat'] == 'gene') &
                (df_contains['sequid'] == gene_last['sequid'])]
            # sort the genes of the segment by starting coord
            df_contains = df_contains.sort_values(by='start').reset_index()
            # compute mean distance between the genes of this colinear segment
            mean_dist = 0
            mean_lens = 0
            for i in range(0, len(df_contains)-1):
                mean_dist += df_contains['start'][i+1] - df_contains['end'][i]
                mean_lens += df_contains['end'][i] - df_contains['start'][i]
            mean_lens += df_contains['end'][i+1] - df_contains['start'][i+1]
            # to get the mean, divide by total events (i+1)
            mean_dist = int(mean_dist/(i+1))
            mean_lens = int(mean_lens/(i+2))
            df_return.loc[row.name, 'mean_gene_len'] = mean_lens
            df_return.loc[row.name, 'mean_gene_dist'] = mean_dist
            df_return.loc[row.name, 'n_genes'] = i+2
# are you interested in storing the IDs/coords from the genes in this specific segment?
#            for rowi in df_contains.iterrows():
#                rowi=rowi[1]
#                gene_id = rowi['attributes'].split('=')[1].split(';')[0]
#                print({gene_id: (rowi['start'], rowi['end'])})
        return df_return

    def paired_segments(self):
        """
        Create a PAF like file, where two regions of the genome are linked to
        each other. For our case, it links multiplicons that are composed by no
        more than a pair of segments. The columns are:

        sequid_from, start_from, end_from, misc_var_from,
        sequid_to,   start_to,   end_to,   misc_var_to
        """
        # init a dataframe to create a PAF-structured object
        df_return = pd.DataFrame(columns=[
            # columns to be filled by first segment
            'sequid_from', 'start_from', 'end_from', 'gn_from', 'gL_from', 'gD_from',
            # columns to be filled by second segment
            'sequid_to', 'start_to', 'end_to', 'gn_to', 'gL_to', 'gD_to'])
        # keep track of how many multiplicons are size 2 or >2:
        size_multiplicons = {'2':0, '>2':0}
        for multiplicon in self.df_segments.groupby('multiplicon'):
            multiplicon=multiplicon[1]
            if len(multiplicon) == 2:
                size_multiplicons['2'] += 1
                f = [multiplicon.iloc[0,[1,5,6,7,8,9]]]
                t = [multiplicon]
                # append the couple of segments to the previous dataframe:
                # index with the multiplicon number
                df_return.loc[multiplicon['multiplicon'].unique()[0]] = (
                    list(multiplicon.iloc[0,[1,5,6,7,8,9]])+
                    list(multiplicon.iloc[1,[1,5,6,7,8,9]]))
            else:
                size_multiplicons['>2'] += 1
                print('Status: multiplicon',
                      multiplicon['multiplicon'].unique()[0],
                      'has more than one segment (will not plot)')

        return df_return, size_multiplicons

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Si es crida com a script:
if __name__ == '__main__':

    ### PARSE ARGUMENTS ###
    import argparse

    parser = argparse.ArgumentParser(description='Add information from GFF3 to'+
                    'segments.txt. Stipulate these files in any order.')
    # file-name: positional arg.
    parser.add_argument('file1', type=str, help='Path to required file-name')
    parser.add_argument('file2', type=str, help='Path to required file-name')
    # integer argument
#    parser.add_argument('-a', '--numero_a', type=int, help='Paràmetre "a"')
    # choices argument
#    parser.add_argument('-o', '--operacio', 
#            type=str, choices=['suma', 'resta', 'multiplicacio'],
#            default='suma', required=False,
#            help='')

    args = parser.parse_args()
    # call a value: args.operacio or args.filename.

    seg = Segments(args.file1, args.file2)
    pd.set_option('display.min_rows', 30)
    print(seg.df_segments)
    print(seg.df_links)
    print(seg.df_links_ratio)
    seg.df_links_ratio.to_csv('segments_linked_pairs.tsv',
                              sep='\t',)


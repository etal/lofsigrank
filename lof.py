#!/usr/bin/env python

from __future__ import print_function

import collections
import random
import sys

import pandas


# ____________________________________________________________________________
# Step 0: Permutation

def permute_table(dtable):
    """Permute a mutation data table's gene, sample and NMAF columns."""
    shuffle_field(dtable, 'gene')
    shuffle_field(dtable, 'sample')
    shuffle_field(dtable, 'Normalized')
    if 'Filler' in dtable:
        del dtable['Filler']


def shuffle_field(dframe, field):
    """Shuffle a column of a pandas DataFrame in-place."""
    column = list(dframe[field])
    random.shuffle(column)
    dframe[field] = column


# _____________________________________________________________________________
# Step 1: Calculate gene-level mutational statistics

def read_list(fname):
    """Parse a "list" file of one string per line."""
    with open(fname) as handle:
        items = [line.strip() for line in handle]
    return items


def make_lof_table(data_table, my_genes, my_samples):
    """Calculate gene-level mutational statistics from a table of mutations.

    Input: nested dict of genes -> samples -> list of mut. type, NMAF, Polyphen
    Output: table stratifying the mutational status of a gene in each sample.

    The output table has a row for each gene and a column for each sample, in
    which there is a number ranging from 0-2 that corresponds to the estimated
    number of alleles lost in the sample. This value is calculated by summing
    the normalized mutant alleles frequencies (NMAF) of all non-synonymous
    mutations striking the gene in this sample, capped at 2.  In addition, the
    final 9 columns of output are the counts of each mutation type (not weighted
    by MAF).

    This output is used as input to Step 2 to calculate the LOF burden.
    """
    # Header
    yield ["Gene"] + my_samples + [
        "Missense:Benign", "Missense:Possibly", "Missense:Probably",
        "MissenseNA", "Indel", "Nonsense", "Frameshift", "Splice-site",
        "Synonymous"]

    gs_lookup = group_data_by_gs(data_table)
    for gene in my_genes:
        synonymous = missense_benign = missense_possibly = missense_probably = \
                missense_na = frameshift = nonsense = splice = indel = 0

        out_row = [gene]
        for sample in my_samples:
            normalized = [0]
            # Count mutations of each type for this gene and sample
            for entry in gs_lookup[gene][sample]:
                if entry['muttype'] == 'Silent':
                    synonymous += 1
                    continue
                if entry['muttype'] == 'Intron':
                    # Shouldn't be here; ignore
                    continue

                if entry['muttype'] == 'Missense_Mutation':
                    if entry['consequence'] == 'benign':
                        missense_benign += 1
                    elif entry['consequence'] == 'possibly':
                        missense_possibly += 1
                    elif entry['consequence'] == 'probably':
                        missense_probably += 1
                    elif entry['consequence'] == 'NA':
                        missense_na += 1
                    else:
                        print("Unhandled missense consequence level:",
                              entry['consequence'], file=sys.stderr)
                elif entry['muttype'] == 'Nonsense_Mutation':
                    nonsense += 1
                elif entry['muttype'] == 'Splice_Site':
                    splice += 1
                elif entry['muttype'] in ('Frame_Shift_Ins', 'Frame_Shift_Del'):
                    frameshift += 1
                elif entry['muttype'] in ('In_Frame_Ins', 'In_Frame_Del'):
                    indel += 1
                else:
                    print("Unhandled mutation type:", entry['muttype'],
                          file=sys.stderr)
                    continue

                normalized.append(entry['normalized'])
            # Add up the normalized mutation counts for this gene and sample
            out_row.append(min(2, sum(normalized)))
        out_row.extend((missense_benign, missense_possibly, missense_probably,
                        missense_na, indel, nonsense, frameshift, splice,
                        synonymous))
        yield out_row


def group_data_by_gs(data_table):
    """Group relevant fields in a data table by gene and sample."""
    gene_data = collections.defaultdict(lambda: collections.defaultdict(list))
    for _idx, row in data_table.iterrows():
        samp = row['sample']
        gene = row['gene']
        gene_data[gene][samp].append({
            'muttype': row['type'].strip(),
            'normalized': row['Normalized'], # NMAF in the manuscript
            'consequence': row['MissenseConsequence'].strip(),
        })
    return gene_data


def rows2dframe(rows):
    """Convert an iterable of table rows to a pandas.DataFrame."""
    header = next(rows)
    return pandas.DataFrame.from_records(rows, columns=header)


# _____________________________________________________________________________

def main(args):
    """."""
    genes = read_list(args.genes)
    samples = read_list(args.samples)
    data_table = pandas.read_table(args.data_table, na_filter=False)

    if args.permute:
        permute_table(data_table)

    lof_table = rows2dframe(make_lof_table(data_table, genes, samples))
    print("Processed", len(lof_table.values), "genes in data table",
          file=sys.stderr)
    lof_table.to_csv(sys.stdout, sep='\t', index=False)


if __name__ == '__main__':
    import argparse
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument('data_table',
                    help="""Mutation data table with NMAF values and Polyphen-2
                    predictions (Data.txt)""")
    AP.add_argument('-g', '--genes', default="Genes.txt",
                    help="List of gene names, one per line (Genes.txt")
    AP.add_argument('-s', '--samples', default="Samples.txt",
                    help="List of sample names, one per line (Samples.txt")
    AP.add_argument('-p', '--permute', action='store_true',
                    help="""Permute the input data table to simulate the
                    background mutation frequencies.""")
    main(AP.parse_args())

#!/usr/bin/env python

"""Step 1: Calculate gene-level mutational statistics from a table of mutations.

Input: table of individual mutations and their Polyphen annotations.
Output: table stratifying the mutational status of a gene in each sample.

- In this table, there is a number ranging from 0-2 that corresponds to the
  estimated number of alleles lost in a given sample. This value is calculated
  by summing the normalized mutant alleles frequencies of all non-synonymous
  mutations striking a gene in a given sample and capping them at 2.
- Additionally, writes the distribution of all mutation types (not weighted
  by MAF) in the final 8 columns of output.

The output of this script is used as input to calculate the LOF burden.

Source: LOF.pl and LOFpermute.pl
"""
from __future__ import print_function

import sys
import collections


def read_list(fname):
    """Parse a "list" file of one string per line."""
    with open(fname) as handle:
        items = [line.strip() for line in handle]
    return items


def read_data(fname):
    """Group relevant fields in (Permuted)Data.txt by gene and sample."""
    gene_data = collections.defaultdict(lambda: collections.defaultdict(list))
    with open(fname, 'rU') as handle:
        lines = iter(handle)
        next(lines)  # Skip the header
        for line in lines:
            (samp, muttype, gene, _contig, norm, _start, _end, _ref, _alt,
             consequence) = line.split('\t', 10)[:10]
            gene_data[gene][samp].append({
                'muttype': muttype.strip(),
                'normalized': float(norm),  # NMAF in the manuscript
                'consequence': consequence.strip(),
            })
    print("Parsed", len(gene_data), "genes from", fname, file=sys.stderr)
    return gene_data


def main(args):
    """."""
    # __________________________________________
    # Step 1

    data_lookup = read_data(args.data)
    my_genes = read_list(args.genes)
    my_samples = read_list(args.samples)

    out_header = [""] + my_samples + [
        "Missense:Benign", "Missense:Possibly", "Missense:Probably",
        "MissenseNA", "Nonsense", "Frameshift", "Splice-site", "Synonymous"]
    print(*out_header, sep='\t')

    for gene in my_genes:
        synonymous = missense_benign = missense_possibly = missense_probably = \
                missense_na = frameshift = nonsense = splice = 0

        out_row = [gene]
        for sample in my_samples:
            normalized = [0]
            # Count mutations of each type for this gene and sample
            for entry in data_lookup[gene][sample]:
                if entry['muttype'] == 'Silent':
                    synonymous += 1
                    continue
                if entry['muttype'] == 'Intron':
                    # XXX
                    continue

                if entry['muttype'] == 'Missense_Mutation':
                    # {'benign': missense_benign,
                    #  'possibly': missense_possibly,
                    #  'probably': missense_probably,
                    #  'NA': missense_na}[entry['consequence'] += 1
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
                    # XXX
                    pass
                else:
                    print("Unhandled mutation type:", entry['muttype'],
                          file=sys.stderr)
                    continue

                normalized.append(entry['normalized'])
            # Add up the normalized mutation counts for this gene and sample
            out_row.append(min(2, sum(normalized)))
        out_row.extend((missense_benign, missense_possibly, missense_probably,
                        missense_na, nonsense, frameshift, splice, synonymous))
        print(*out_row, sep='\t')


if __name__ == '__main__':
    import argparse
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument('data',
                    help="Mutation data table with Polyphen scores (Data.txt)")
    AP.add_argument('-g', '--genes', default="Genes.txt",
                    help="List of gene names, one per line (Genes.txt")
    AP.add_argument('-s', '--samples', default="Samples.txt",
                    help="List of sample names, one per line (Samples.txt")
    main(AP.parse_args())

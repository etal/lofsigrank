#!/usr/bin/env python

"""Identify significantly mutated genes in a set of many WES samples."""
from __future__ import print_function, division

import itertools
import sys

import pandas

from lof import read_list


def shainsig(table, samples):
    mut_probdam = 'Missense:Probably'
    mut_syn = 'Synonymous'
    mut_trunc = ['Nonsense', 'Frameshift', 'Splice-site']
    mut_other = ['Missense:Benign', 'Missense:Possibly', 'MissenseNA']
    mut_all = [mut_probdam, mut_syn] + mut_trunc + mut_other

    # _________________________________________________________________________
    # 1. Calculate the global nonsynonymous:synonymous ratio

    # Within each mutation category, sum counts (across all genes)
    tot_count_probdam = sum(table[mut_probdam])
    tot_count_syn = sum(table[mut_syn])
    tot_count_trunc = sum(itertools.chain(*(list(table[col])
                                            for col in mut_trunc)))
    tot_count_other = sum(itertools.chain(*(list(table[col])
                                            for col in mut_other)))

    # Global mutation count across all categories and genes (= 3504)
    tot_count_all = sum((tot_count_probdam, tot_count_syn, tot_count_trunc,
                         tot_count_other))
    print("Counted", tot_count_all, "mutations across", len(table), "genes",
          "and", len(samples), "samples", file=sys.stderr)

    # Fraction of global mutations in each category of interest
    tot_frac_probdam = tot_count_probdam / tot_count_all
    tot_frac_syn = tot_count_syn / tot_count_all
    tot_frac_trunc = tot_count_trunc / tot_count_all

    # Global nonsynonymous:synonymous ratio = (1-syn)/syn (= 2.13697)
    tot_ns_s_ratio = (1 - tot_frac_syn) / tot_frac_syn

    # _________________________________________________________________________
    # 2. Calculate each gene's mutation score
    for _idx, row in table.iterrows():
        gene_count_all = sum([row[col] for col in mut_all])
        if not gene_count_all:
            # Gene is not mutated at all --> zero score
            yield (row['Gene'], 0.0)
            continue

        # Initial score is the sum the 'Normalized' values across all samples
        raw_score = sum(row[sid] for sid in samples)

        # Adjust for NS:S ratio
        gene_count_syn = row[mut_syn]
        syn_factor = max(1 - tot_ns_s_ratio * gene_count_syn / gene_count_all,
                         0)
        new_score = raw_score * syn_factor

        # Adjust for "probably damaging" missense and truncating mutations
        gene_frac_probdam = row[mut_probdam] / gene_count_all
        probdam_factor = 1 + gene_frac_probdam - tot_frac_probdam
        gene_frac_trunc = sum([row[col] for col in mut_trunc]) / gene_count_all
        trunc_factor = gene_frac_trunc / tot_frac_trunc
        final_score = new_score * probdam_factor * trunc_factor
        yield (row['Gene'], final_score)


def main(args):
    """."""
    table = pandas.read_table(args.lof_table, na_filter=False)
    samples = read_list(args.samples)
    results = shainsig(table, samples)
    for gene, score in sorted(results, key=lambda pair: pair[1]):
        print(gene, score, sep='\t')


if __name__ == '__main__':
    import argparse
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument('lof_table', help="LOF data table (TSV)")
    AP.add_argument('-s', '--samples', default="Samples.txt",
                    help="List of sample names, one per line (Samples.txt")
    # AP.add_argument('-p', '--permutations', type=int, default=7,
    #                 help="""Number of times to permute the input data to
    #                 simulate the background mutation frequencies.""")
    main(AP.parse_args())

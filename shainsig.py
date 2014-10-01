#!/usr/bin/env python

"""Identify significantly mutated genes in a set of many WES samples.

Prints a table of each gene's observed and expected mutation burdens and
estimated FDR for predicted tumor suppressors.
"""
from __future__ import print_function, division

import collections
import itertools
import random
import sys

import pandas
import numpy


# _____________________________________________________________________________
# Step_1: Calculate gene-level mutational statistics

def make_lof_table(data_table, my_genes, my_samples, summary_func):
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
            # Summarize the normalized mutation counts for this gene and sample
            out_row.append(summary_func(normalized))
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
    return pandas.DataFrame.from_records(list(rows), columns=header)


# _____________________________________________________________________________
# Step_2: Identify significantly mutated genes

def shainsig(table, samples, verbose=True):
    """Identify significantly mutated genes in the processed LOF table."""
    mut_probdam = 'Missense:Probably'
    mut_syn = 'Synonymous'
    mut_trunc = ['Nonsense', 'Frameshift', 'Splice-site']
    mut_other = ['Missense:Benign', 'Missense:Possibly', 'MissenseNA', 'Indel']
    mut_all = [mut_probdam, mut_syn] + mut_trunc + mut_other

    # Calculate the global nonsynonymous:synonymous ratio ---------------------
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
    if verbose:
        print("Counted", tot_count_all, "mutations across", len(table), "genes",
            "and", len(samples), "samples", file=sys.stderr)

    # Fraction of global mutations in each category of interest
    tot_frac_probdam = tot_count_probdam / tot_count_all
    tot_frac_syn = tot_count_syn / tot_count_all
    tot_frac_trunc = tot_count_trunc / tot_count_all

    # Global nonsynonymous:synonymous ratio = (1-syn)/syn (= 2.13697)
    tot_ns_s_ratio = (1 - tot_frac_syn) / tot_frac_syn

    # Calculate each gene's mutation score ------------------------------------
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


# _____________________________________________________________________________
# Step_3: False Discovery Rate (FDR) calculation

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

def read_list(fname):
    """Parse a "list" file of one string per line."""
    with open(fname) as handle:
        items = [line.strip() for line in handle]
    return items


def main(args):
    """Run the script."""
    genes = read_list(args.genes)
    samples = read_list(args.samples)
    data_table = pandas.read_table(args.data_table, na_filter=False)
    summary_function = {'sumcap': lambda x: min(2, sum(x)),
                        'mean': numpy.mean,
                        'max': max}[args.function]

    # Step_1
    lof_table = rows2dframe(make_lof_table(data_table, genes, samples,
                                           summary_function))
    print("Processed", len(lof_table.values), "genes in data table",
          file=sys.stderr)

    # Step_2
    results = list(shainsig(lof_table, samples))
    results.sort(key=lambda pair: pair[1])

    # Step_3
    if args.permutations:
        # Calculate gene score percentiles
        orig_pctiles = numpy.arange(1, 0, -1. / len(results))

        # Calculate percentiles for simulated "background" scores
        perm_scores = []
        print("Permuting mutation data", args.permutations, "times:", end=' ',
              file=sys.stderr)
        for idx in range(args.permutations):
            print(idx + 1, end=' ', file=sys.stderr)
            permute_table(data_table)
            ptable = rows2dframe(make_lof_table(data_table, genes, samples,
                                                summary_function))
            perm_scores.extend(s for g, s in shainsig(ptable, samples, False))
        perm_scores = numpy.asfarray(sorted(perm_scores))
        perm_pctiles = numpy.arange(1, 0, -1. / len(perm_scores))
        print("\nMax permutation score:", max(perm_scores), file=sys.stderr)

        # Do the math! Output a table!
        print("Gene", "Obs.Score", "Obs.Pctile", "Sim.Score", "Sim.Pctile", "FDR", sep='\t')
        perm_pctiles_rev = perm_pctiles[::-1]
        for (gene, obs_score), obs_pctile in zip(results, orig_pctiles):
            score_rank = perm_scores.searchsorted(obs_score)
            if score_rank == len(perm_scores):
                exp_pctile = 0
                fdr = 0.0
            else:
                exp_pctile = perm_pctiles[score_rank]
                # FDR: % false positives / % true positives
                fdr = min(1.0, exp_pctile / obs_pctile)
            exp_score = perm_scores[len(perm_scores) - 1 -
                                    perm_pctiles_rev.searchsorted(obs_pctile)]
            print(gene, obs_score, obs_pctile, exp_score,
                  exp_pctile, fdr,
                  sep='\t')

    else:
        print("Gene", "Score", sep='\t')
        for gene, score in results:
            print(gene, score, sep='\t')


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
    AP.add_argument('-p', '--permutations', type=int, default=20,
                    help="""Number of times to permute the input data to
                    simulate the background mutation frequencies.""")
    AP.add_argument('-f', '--function', default='sumcap',
                    choices=['sumcap', 'max', 'mean'],
                    help="Summary function for gene-level NMAF counts.")
    main(AP.parse_args())

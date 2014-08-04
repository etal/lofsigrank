#!/usr/bin/env python

"""Permute a mutation data table's gene, sample and NMAF columns."""

import random
import sys

import pandas


def shuffle_field(dframe, field):
    """Shuffle a column of a pandas DataFrame in-place."""
    column = list(dframe[field])
    random.shuffle(column)
    dframe[field] = column


fname = sys.argv[1]
data = pandas.read_table(fname, na_filter=False)
shuffle_field(data, 'gene')
shuffle_field(data, 'sample')
shuffle_field(data, 'Normalized')

if 'Filler' in data:
    del data['Filler']

data.to_csv(sys.stdout, sep='\t', index=False)

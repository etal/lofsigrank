ShainSig
========

Identify potential tumor suppressors from the relative burden of
loss-of-function mutations. 

Usage 
-----

To run the simulation with 100 iterations:

    python shainsig.py Data.csv -g Genes.txt -s Samples.txt -p 100 > out.tsv

Input files are included in this directory.

Output is a tab-delimited table of:

- Gene name
- Observed score, based on mutation counts
- Observed percentile of this gene, ranked among all genes in Samples.txt
- Expected score based on simulations of random mutations
- Expected percentile based on simulations of random mutations
- False discovery rate (FDR), the ratio of false positives (expected percentile)
  to true positives (observed percentile)


# Somatic variants and deep mutational scanning analysis to map protein-level consequences of different mutational signatures

This repository contains all code (python, R markdown workbooks, bash scripts) used to perform this analysis.

## Content

`Analysis`: R markdown and rendered HTML workbooks containing code and output of analyses presented in the manuscript.

`CountContexts`: Python code used to count number of codons in coding sequences containing each trinucleotide context separately for protein surface, core and interacting interface.

`CountContexts_Species`: Python code used to count number of codons in coding sequences containing each trinucleotide context separately for six different species.

`SigEnrich`: Python code used to calculate the enrichment of given mutational signature contexts across all human protein surface, core and interacting interface.

`VariantProcessing`: Python/R code used to extract relevant columns from TCGA Mutation Annotation Files (MAFs) and annotate trinucleotide mutational contexts and signatures.

Please consult README files in each subfolder for more detailed descriptions of files within each subfolder.

## Citation

Ng JCF & Fraternali F. Meta-analysis of mutations and deep mutational scanning data reveal protein-level consequences of somatic mutational signatures. *Manuscript in preparation*, 2021.


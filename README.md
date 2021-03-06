# Somatic variants and deep mutational scanning analysis to identify & predict damaging protein variants

This repository contains all code (python, R markdown workbooks, bash scripts) used to perform this analysis.

All data files are included in this repository or in the Zenodo repository ([link here](https://doi.org/10.5281/zenodo.4726168)) accompanying this manuscript.

## Content

`Analysis`: R markdown and rendered HTML workbooks containing code and output of analyses presented in the manuscript.

`CountContexts`: Python code used to count number of codons in coding sequences containing each trinucleotide context separately for protein surface, core and interacting interface.

`CountContexts_Species`: Python code used to count number of codons in coding sequences containing each trinucleotide context separately for six different species.

`SigEnrich`: Python code used to calculate the enrichment of given mutational signature contexts across all human protein surface, core and interacting interface.

`VariantProcessing`: Python/R code used to extract relevant columns from TCGA Mutation Annotation Files (MAFs) and annotate trinucleotide mutational contexts and signatures.

`Prediction`: Python code used to assemble and annotate data, and training/evaluating predictors trained on deep mutational scanning data with different feature combinations (AA, DNA Mutational signature, conservation, structure).

Please consult README files in each subfolder for more detailed descriptions of files within each subfolder.

## Citation

Ng JCF & Fraternali F. The 'dark matter' of protein variants carries a distinct DNA signature and improves the prediction of damaging variant effects, 2022.


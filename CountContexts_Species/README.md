# Counting number of codons in coding sequences from different species containing each trinucleotide context

A set of python functions have been written to count the number of codons in coding sequences (obtained from UniProt) from six different species, these codons contain a given trinucleotide context.

## Contents

`cds_noStructMap.py`: The python code which performs this analysis.

`cds_allcontexts_counts.bash`: A bash script which processes all the coding sequence files and perform this counting.

`*.fasta`: Coding sequences from UniProt used in this analysis. The Human file is a symbolic link to the fasta file stored in the `CountContexts` subfolder of this repository.

This python code generates a CSV with the requested counts.

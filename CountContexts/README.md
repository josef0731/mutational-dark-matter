# Counting number of codons in coding sequences containing each trinucleotide context

A set of python functions have been written to count the number of codons in coding sequences (obtained from UniProt) which (1) are localised to a given protein structural region (surface/core/interface) and (2) contain a given trinucleotide context.

## Contents

`cds_StructMap_context.py`: The python code which performs this analysis.

`UP000005640_9606_DNA.fasta`: Coding sequences obtained from UniProt used in this analysis.

`zoomvarCDS_32contexts.bash`: A bash script to invoke the python code to count a user-supplied trinucleotide context.

`zoomvarCDS_32contexts_submit.bash`: A bash script which loops through `zoomvarCDS_32contexts.bash` and considers all 32 trinucleotide contexts.

These code require a tab-separated file with annotation of structural and interaction mapping from ZoomVar on each individual amino acid positions over all human proteins. This is available in the data repository on zenodo accompanying this manuscript.

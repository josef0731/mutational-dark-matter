# Scripts to process Mutation Annotation Files

This set of scripts process, clean and annotate TCGA Mutation Annotation Files downloaded from GDAC Firehose.

## Contents

`variant_processing.*`: R markdown and rendered HTML workbook detailing steps performed to process the variant data. The code require a tab-separated file with annotation of structural and interaction mapping from ZoomVar on each individual amino acid positions over all human proteins. This is available in the data repository on zenodo accompanying this manuscript.

`extractMAF.py`: Python functions to extract relevant columns from the TCGA Mutation Annotation Files.

`extractMAF.sh`: Bash script used to invoke `extractMAF.py`. Columns extracted are indicated in the python call.

`classify_signatures.py`: Python functions to clean variants to parse possible doublet-nucleotide substitutions, and annotate mutational signatures (Aging, APOBEC, POLE, UV, 5-FU, Platinum) based on the DNA sequence contexts.


# Enrichment of given mutational signature contexts across all human protein surface, core and interacting interface.

A set of python functions and code to calculate the enrichment of DNA motifs for a given mutational signature (Aging, APOBEC, POLE, UV, 5-FU, Platinum) over codons mapped to protein surface/core/interacting interfaces of all human proteins.

## Contents

### Annotation of per-protein availability of trinucleotide motifs over protein surface/core/interacting interface

`proteinCDS_struct_5fu_NEW.csv`: Availability of 5-FU-preferred DNA motifs over all human proteins and respective structural regions.
`proteinCDS_struct_aging_NEW.csv`: Availability of aging-preferred DNA motifs over all human proteins and respective structural regions.
`proteinCDS_struct_apobec_NEW.csv`: Availability of APOBEC-preferred DNA motifs over all human proteins and respective structural regions.
`proteinCDS_struct_platinum_NEW.csv`: Availability of Platinum-therapy-preferred DNA motifs over all human proteins and respective structural regions.
`proteinCDS_struct_pole_NEW.csv`: Availability of POLE-preferred DNA motifs over all human proteins and respective structural regions.
`proteinCDS_struct_uv_NEW.csv`: Availability of UV-preferred DNA motifs over all human proteins and respective structural regions.

### Python code

`signature_structEnrich4.py`: Python code to perform the calculation.

### Case lists

`TCGA_5fu.txt`: List of TCGA cases which underwent 5-FU treatment.

`TCGA_apobec.txt`: List of TCGA cases with more-than-expected APOBEC signature enrichment.

`TCGA_hypermutated_POLEonly.txt`: List of POLE-mutated TCGA cases which are of the hypermutated state.

`TCGA_over70.txt`: List of TCGA cases where age of cancer onset is at least 70.

`TCGA_platinum.txt`: List of TCGA cases which underwent platinum-based therapies.

`TCGA_SKCM.txt`: List of TCGA cases which are in the skin cutaneous melanoma (SKCM) cohort.

### Mutations

`TCGA_maf_allmuttype_mapped_withSig_new.txt`: Mutation Annotation File (MAF) with all TCGA cohorts merged together. Each mutation is annotated to be attributable to any of the 6 mutational signatures.

### Bash script

`TCGA_sig_zoomvar2.sh`: Bash script to invoke the python code.

`TCGA_sig_zoomvar_submit2.sh`: Bash script to submit jobs on a slurm cluster to calculate enrichment for each of the 6 mutational signatures.

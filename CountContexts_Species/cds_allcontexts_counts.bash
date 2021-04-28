#!/bin/bash

# count number of contexts (all 32) in CDS of proteomes from selected organisms

folder="."
CDS=(Celegans_UP000001940_6239_DNA.fasta Drerio_UP000000437_7955_DNA.fasta Mouse_UP000000589_10090_DNA.fasta Chimp_UP000002277_9598_DNA.fasta Human_UP000005640_9606_DNA.fasta Yeast_UP000002311_559292_DNA.fasta)

for cds in ${CDS[@]}; do
	organism=${cds%_UP*.fasta}
	echo $organism
	python cds_noStructMap.py $folder/$cds $organism"_CDS_all32contexts_counts.csv" --num_core 1
done

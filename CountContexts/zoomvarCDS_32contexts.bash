#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16000
module load devtools/anaconda/2019.3-python2.7.16

python -u cds_StructMap_context.py protein_position_annotations_zoomvar_category_NEW.tsv \
	UP000005640_9606_DNA.fasta proteinCDS_zoomvar_`echo $context1`.csv \
	--sig $context1 $context2 --num_core 1

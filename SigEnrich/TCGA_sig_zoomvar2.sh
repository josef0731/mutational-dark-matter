#!/bin/bash
#SBATCH --time=4-0:00
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=64000
module load devtools/anaconda/2019.3-python2.7.16 

# Run thing that uses multiple cores 
python -u signature_structEnrich4.py --maf $input --sigCol $sigCol --sigVal 1 \
	--countsTb $countsTb --jobName $jobName \
	--filter all --cohortfile TCGA_sigs_stats_new.txt --caseList $caseList

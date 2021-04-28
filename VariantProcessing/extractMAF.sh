#!/bin/bash

#$ -S /bin/sh

# Work at current directory.
#$ -cwd

# How much memory do you need **per core**?
#$ -l h_rt=04:00:00,h_vmem=16G

# Number of cores you need 
#$ -pe smp 1

#load modules
module load general/python/2.7.10

#extract useful columns from TCGA mafs
for f in `ls ../TCGA.firehose_maf/*_maf.txt`; do
    echo $f
    outloc=${f#../TCGA.firehose_maf/}
    outloc=${outloc%_maf.txt} 
    python extractMAF.py ../TCGA.firehose_maf/$f proteinCDS_zoomvar_apobec_proteinLevel.csv ../TCGA.firehose_maf/`echo $outloc`_maf_part.txt Hugo_Symbol Chromosome Start_Position End_position Strand Variant_Classification Reference_Allele Tumor_Seq_Allele1 Tumor_Seq_Allele2 APOBEC_mutation "CONTEXT(+/-20)" Tumor_Sample_Barcode Matched_Norm_Sample_Barcode Match_Norm_Seq_Allele1 Match_Norm_Seq_Allele2 Annotation_Transcript Protein_Change SwissProt_acc_Id UniProt_AApos i_HGNC_UniProtIDsuppliedbyUniProt i_dbNSFP_Uniprot_acc i_dbNSFP_Uniprot_id i_dbNSFP_Polyphen2_HVAR_score i_dbNSFP_SIFT_score 
done



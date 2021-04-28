#!/bin/bash

sbatch -p brc --export=input=TCGA_maf_allmuttype_mapped_withSig_new.txt,sigCol=PlatSig,countsTb=proteinCDS_struct_platinum_NEW.csv,caseList=TCGA_platinum.txt,jobName=plat --job-name=platinum TCGA_sig_zoomvar2.sh
sbatch -p brc --export=input=TCGA_maf_allmuttype_mapped_withSig_new.txt,sigCol=ApobecSig,countsTb=proteinCDS_struct_apobec_NEW.csv,caseList=TCGA_apobec.txt,jobName=apobec --job-name=apobec TCGA_sig_zoomvar2.sh
sbatch -p brc --export=input=TCGA_maf_allmuttype_mapped_withSig_new.txt,sigCol=AgingSig,countsTb=proteinCDS_struct_aging_NEW.csv,caseList=TCGA_over70.txt,jobName=aging --job-name=aging TCGA_sig_zoomvar2.sh
sbatch -p brc --export=input=TCGA_maf_allmuttype_mapped_withSig_new.txt,sigCol=FiveFUSig,countsTb=proteinCDS_struct_5fu_NEW.csv,caseList=TCGA_5fu.txt,jobName=FiveFU --job-name=FiveFU TCGA_sig_zoomvar2.sh
sbatch -p brc --export=input=TCGA_maf_allmuttype_mapped_withSig_new.txt,sigCol=POLESig,countsTb=proteinCDS_struct_pole_NEW.csv,caseList=TCGA_hypermutated_POLEonly.txt,jobName=POLE --job-name=POLE TCGA_sig_zoomvar2.sh
sbatch -p brc --export=input=TCGA_maf_allmuttype_mapped_withSig_new.txt,sigCol=UVSig,countsTb=proteinCDS_struct_uv_NEW.csv,jobName=UV,caseList=TCGA_SKCM.txt --job-name=UV TCGA_sig_zoomvar2.sh


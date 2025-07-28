#!/bin/bash
#SBATCH --job-name=PRS_inf_UGLI2_script1_CS
#SBATCH --nodes=1 --tasks-per-node=1 --time=12:00:00
#SBATCH -o PRS_inflammation_UGLI2_script1.out

#-------------------------------------------------------

cd ## ** ADD WORKING DIR HERE ** ##

#-------------------------------------------------------
module load PLINK/1.9-beta6-20190617 

## location of genetic data
UGLI2="## ADD LOCATION IMPUTED GENETIC DATA HERE ##"

## extract snps
for i in {1..22}
do
plink --vcf ${UGLI2}/${i}.UGLI2.vcf.gz --extract said_ahluw_bor_rosa_swer_sar_SNPs_duplicatesremoved_lifelines.txt --out bg_chr${i}_adult --make-bed
done

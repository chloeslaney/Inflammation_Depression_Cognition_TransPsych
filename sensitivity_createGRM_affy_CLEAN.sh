#!/bin/bash
#SBATCH --job-name=affy_grm
#SBATCH --nodes=1 --tasks-per-node=1 
#SBATCH --mem=60G
#SBATCH --time=26:00:00
#SBATCH -o affy_grm_final_step4_gcta_relcutoff.out

#--------------------------------------------------------

cd ## ** ADD DIR HERE ** ##

#--------------------------------------------------------
## module load
module load PLINK/1.9-beta6-20190617
module load cluster-utils
module list
ctop

## affy geno location
affy_geno= ## ** ADD DIR FOR GENO DATA HERE ** ##

## Step 1: quality control
## time: mins, mem=1GB (2chr)
for i in {1..22} 
do
plink --bfile ${affy_geno}/chr_${i} --maf 0.01 --geno 0.05 --mind 0.05 --hwe 1e-6 --snps-only --make-bed --out chr${i}_affy
done

## Step2: combine to one dataset
## time: mins, mem=1GB (2chr)
plink --bfile chr1_affy --merge-list chr_list.txt --make-bed --out affy_snps_combined_final 

## Step3: keep independent snps 
plink --bfile affy_snps_combined_final --indep-pairwise 50 10 0.1 --make-bed --out affy_snps_for_grm_final 

# Step4: create genetic relationship matrix
plink --bfile affy_snps_for_grm_final --exclude affy_snps_for_grm_final.prune.out --make-grm-bin --out GRM_AFFY_GCTA
plink --bfile affy_snps_for_grm_final --exclude affy_snps_for_grm_final.prune.out --make-grm-bin --rel-cutoff 0.125 --out GRM_AFFY_GCTA_RELCUTOFF_0125
plink --bfile affy_snps_for_grm_final --exclude affy_snps_for_grm_final.prune.out --make-grm-bin --rel-cutoff 0.25 --out GRM_AFFY_GCTA_RELCUTOFF_025
plink --bfile affy_snps_for_grm_final --exclude affy_snps_for_grm_final.prune.out --make-grm-bin --rel-cutoff 0.5 --out GRM_AFFY_GCTA_RELCUTOFF_05


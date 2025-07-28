#!/bin/bash
#SBATCH --job-name=PRS_inf_UGLI2_script3_CS
#SBATCH --nodes=1 --tasks-per-node=1 --time=12:00:00
#SBATCH -o PRS_inflammation_UGLI2_script3.out

#---------------------------------------------------

cd ## ** ADD WORKING DIR HERE ** ##

#---------------------------------------------------
module load PLINK/1.9-beta6-20190617

## Use --score file.txt no-mean-imputation
## Generate sum scores = file inc var id, EA, beta - ensure all positive
## creating sum score (beta for allele x number of risk allele) - check no missing snps (CNT col) 
plink --bfile all_adult_ugli2_missnpremoved --score-no-mean-imputation --score saidcis_for_plink.txt sum --out saidcis_score_sum_ugli2
plink --bfile all_adult_ugli2_missnpremoved --score-no-mean-imputation --score said_genomewide_for_plink.txt sum --out said_genomewide_score_sum_ugli2
plink --bfile all_adult_ugli2_missnpremoved --score-no-mean-imputation --score ahluwaliacis_for_plink.txt sum --out ahluwaliacis_score_sum_ugli2
plink --bfile all_adult_ugli2_missnpremoved --score-no-mean-imputation --score ahluwalia_genomewide_for_plink.txt sum --out ahluwalia_genomewide_score_sum_ugli2
plink --bfile all_adult_ugli2_missnpremoved --score-no-mean-imputation --score swerdlow_for_plink.txt sum --out swerdlow_score_sum_ugli2
plink --bfile all_adult_ugli2_missnpremoved --score-no-mean-imputation --score sarwar_for_plink.txt sum --out sarwar_score_sum_ugli2
plink --bfile all_adult_ugli2_missnpremoved --score-no-mean-imputation --score rosa_for_plink.txt sum --out rosa_score_sum_ugli2
plink --bfile all_adult_ugli2_missnpremoved --score-no-mean-imputation --score borges_eafchecked_for_plink.txt sum --out borges_score_sum_ugli2

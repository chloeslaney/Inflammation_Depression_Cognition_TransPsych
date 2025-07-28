#!/bin/bash
#SBATCH --job-name=affy_king
#SBATCH --nodes=1 --tasks-per-node=1 
#SBATCH --mem=60G
#SBATCH --time=26:00:00
#SBATCH -o king_affy_phenoavail.out

#--------------------------------------------------------

cd ## ** ADD WORKING DIR HERE ** ##

#--------------------------------------------------------
## module load
module load PLINK/2.0-alpha2-20191006
module load cluster-utils
module list
ctop

## Get list of unrelated indivdiuals (1st degree, > 2nd degree, > 3rd degree)
## Note: used geometric mean of kinship coefficient cutoffs (e.g., 1stdegree=0.25, 2nd degree=0.125, 3rd degree=0.0625)
## e.g., for 1st degree (0.25 x 0.125) = 0.03125. sqrt(0.03125) = 0.177
plink2 --bfile affy_snps_for_grm_final --exclude affy_snps_for_grm_final.prune.out --keep ugli2_FID_IID.txt --double-id --king-cutoff 0.177 --out KING_UGLI2_CUTOFF_177_phenoavail
plink2 --bfile affy_snps_for_grm_final --exclude affy_snps_for_grm_final.prune.out --keep ugli2_FID_IID.txt --double-id --king-cutoff 0.088 --out KING_UGLI2_CUTOFF_088_phenoavail
plink2 --bfile affy_snps_for_grm_final --exclude affy_snps_for_grm_final.prune.out --keep ugli2_FID_IID.txt --double-id --king-cutoff 0.044 --out KING_UGLI2_CUTOFF_044_phenoavail

#plink2 --bfile affy_snps_for_grm_final --exclude affy_snps_for_grm_final.prune.out --king-cutoff 0.177 --out KING_UGLI2_CUTOFF_177
#plink2 --bfile affy_snps_for_grm_final --exclude affy_snps_for_grm_final.prune.out --king-cutoff 0.088 --out KING_UGLI2_CUTOFF_088
#plink2 --bfile affy_snps_for_grm_final --exclude affy_snps_for_grm_final.prune.out --king-cutoff 0.044 --out KING_UGLI2_CUTOFF_044

#!/bin/bash
#SBATCH --job-name=affy_greml
#SBATCH --nodes=1 --tasks-per-node=1 
#SBATCH --mem=60G
#SBATCH --time=29:00:00
#SBATCH -o affy_grammar_part2.out

#--------------------------------------------------------

cd ## ** ADD DIR HERE ** ##

#--------------------------------------------------------
module load PLINK/1.9-beta6-20190617
module load GCTA/1.93.2beta-GCCcore-7.3.0
module load cluster-utils
module list
ctop

## Create sparse GRM from GRM above (set off-diagonals that are < 0.125 to 0)
gcta64 --grm GRM_AFFY_GCTA --make-bK 0.125 --out GRM_AFFY_GCTA_SPARSE

## GREML analysis with --reml-pred-rand option to output residuals of individuals (predicting phenotype from GRM)
# Note: 30G with 12hrs timed out. Finished up to i=9. Next time run=24hrs, 60G.
for i in {1..17}
do
gcta64 --reml --grm GRM_AFFY_GCTA_SPARSE --pheno ugli2_merged_pheno_for_gcta.txt --mpheno ${i} --reml-pred-rand --qcovar ugli2_merged_pcs_for_gcta.txt --out greml_ugli2_pcsadded_${i}
done




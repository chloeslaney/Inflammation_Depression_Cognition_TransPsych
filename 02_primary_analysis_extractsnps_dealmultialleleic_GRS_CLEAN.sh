#!/bin/bash
#SBATCH --job-name=PRS_inf_UGLI2_script2_part2
#SBATCH --nodes=1 --tasks-per-node=1 --time=12:00:00
#SBATCH -o PRS_inflammation_UGLI2_script2_part2.out

#--------------------------------------------------------

cd ## ** ADD WORKING DIR HERE ** ##

#--------------------------------------------------------
module load PLINK/1.9-beta6-20190617

##After creating edited list for chromosomes containing any multi-allelic snps, run code below with edited chr included
plink --bfile bg_chr1_adult --merge-list allfiles_adult_edited.txt --make-bed --out all_adult_ugli2_missnpremoved

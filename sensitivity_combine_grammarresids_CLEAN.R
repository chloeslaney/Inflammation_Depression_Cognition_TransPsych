## Project: GRS for inflammation markers on outcomes in Lifelines (cognition, depression, anxiety phenotypes)
## Extract residuals from greml for GRS analysis
## Outputs: pheno resids for each chip
## Chloe Slaney

###################
## Load packages ##
###################
library(cat)
library(data.table)
library(ieugwasr)
library(tidyverse)
library(readr)
library(dplyr)
library(utils)
library(httr)
library(rlist)
library(haven)
library(devtools)
library(PubHelper)
library(psych)

##########################################################
#################### PART 1: Load data ###################
##########################################################
setwd("## ADD WORKING DIR HERE ##")
rm(list = ls())

################################################
## A. Load greml files containing resids data ##
################################################
cyto_resids_VA_rt <- fread("greml_cyto_pcsadded_1.indi.blp", header=FALSE)
cyto_resids_PS_rt <- fread("greml_cyto_pcsadded_2.indi.blp", header=FALSE)
cyto_resids_EM_acc <- fread("greml_cyto_pcsadded_3.indi.blp", header=FALSE)
cyto_resids_WM_acc <- fread("greml_cyto_pcsadded_4.indi.blp", header=FALSE)
cyto_resids_WM_rt <- fread("greml_cyto_pcsadded_5.indi.blp", header=FALSE)
cyto_resids_RFFT <- fread("greml_cyto_pcsadded_6.indi.blp", header=FALSE)
cyto_resids_panas_neg <- fread("greml_cyto_pcsadded_7.indi.blp", header=FALSE)
cyto_resids_panas_pos <- fread("greml_cyto_pcsadded_8.indi.blp", header=FALSE)
cyto_resids_minidep <- fread("greml_cyto_pcsadded_9.indi.blp", header=FALSE)
cyto_resids_minianx <- fread("greml_cyto_pcsadded_10.indi.blp", header=FALSE)
cyto_resids_minimdd <- fread("greml_cyto_pcsadded_11.indi.blp", header=FALSE)
cyto_resids_minigad <- fread("greml_cyto_pcsadded_12.indi.blp", header=FALSE)
cyto_resids_hscrp <- fread("greml_cyto_pcsadded_13.indi.blp", header=FALSE)
cyto_resids_bmi <- fread("greml_cyto_pcsadded_14.indi.blp", header=FALSE)
cyto_resids_alcohol <- fread("greml_cyto_pcsadded_15.indi.blp", header=FALSE)
cyto_resids_smoke <- fread("greml_cyto_pcsadded_16.indi.blp", header=FALSE)
cyto_resids_education <- fread("greml_cyto_pcsadded_17.indi.blp", header=FALSE)

gsa_resids_VA_rt <- fread("greml_gsa_pcsadded_1.indi.blp", header=FALSE)
gsa_resids_PS_rt <- fread("greml_gsa_pcsadded_2.indi.blp", header=FALSE)
gsa_resids_EM_acc <- fread("greml_gsa_pcsadded_3.indi.blp", header=FALSE)
gsa_resids_WM_acc <- fread("greml_gsa_pcsadded_4.indi.blp", header=FALSE)
gsa_resids_WM_rt <- fread("greml_gsa_pcsadded_5.indi.blp", header=FALSE)
gsa_resids_RFFT <- fread("greml_gsa_pcsadded_6.indi.blp", header=FALSE)
gsa_resids_panas_neg <- fread("greml_gsa_pcsadded_7.indi.blp", header=FALSE)
gsa_resids_panas_pos <- fread("greml_gsa_pcsadded_8.indi.blp", header=FALSE)
gsa_resids_minidep <- fread("greml_gsa_pcsadded_9.indi.blp", header=FALSE)
gsa_resids_minianx <- fread("greml_gsa_pcsadded_10.indi.blp", header=FALSE)
gsa_resids_minimdd <- fread("greml_gsa_pcsadded_11.indi.blp", header=FALSE)
gsa_resids_minigad <- fread("greml_gsa_pcsadded_12.indi.blp", header=FALSE)
gsa_resids_hscrp <- fread("greml_gsa_pcsadded_13.indi.blp", header=FALSE)
gsa_resids_bmi <- fread("greml_gsa_pcsadded_14.indi.blp", header=FALSE)
gsa_resids_alcohol <- fread("greml_gsa_pcsadded_15.indi.blp", header=FALSE)
gsa_resids_smoke <- fread("greml_gsa_pcsadded_16.indi.blp", header=FALSE)
gsa_resids_education <- fread("greml_gsa_pcsadded_17.indi.blp", header=FALSE)

ugli2_resids_VA_rt <- fread("greml_ugli2_pcsadded_1.indi.blp", header=FALSE)
ugli2_resids_PS_rt <- fread("greml_ugli2_pcsadded_2.indi.blp", header=FALSE)
ugli2_resids_EM_acc <- fread("greml_ugli2_pcsadded_3.indi.blp", header=FALSE)
ugli2_resids_WM_acc <- fread("greml_ugli2_pcsadded_4.indi.blp", header=FALSE)
ugli2_resids_WM_rt <- fread("greml_ugli2_pcsadded_5.indi.blp", header=FALSE)
ugli2_resids_RFFT <- fread("greml_ugli2_pcsadded_6.indi.blp", header=FALSE)
ugli2_resids_panas_neg <- fread("greml_ugli2_pcsadded_7.indi.blp", header=FALSE)
ugli2_resids_panas_pos <- fread("greml_ugli2_pcsadded_8.indi.blp", header=FALSE)
ugli2_resids_minidep <- fread("greml_ugli2_pcsadded_9.indi.blp", header=FALSE)
ugli2_resids_minianx <- fread("greml_ugli2_pcsadded_10.indi.blp", header=FALSE)
ugli2_resids_minimdd <- fread("greml_ugli2_pcsadded_11.indi.blp", header=FALSE)
ugli2_resids_minigad <- fread("greml_ugli2_pcsadded_12.indi.blp", header=FALSE)
ugli2_resids_hscrp <- fread("greml_ugli2_pcsadded_13.indi.blp", header=FALSE)
ugli2_resids_bmi <- fread("greml_ugli2_pcsadded_14.indi.blp", header=FALSE)
ugli2_resids_alcohol <- fread("greml_ugli2_pcsadded_15.indi.blp", header=FALSE)
ugli2_resids_smoke <- fread("greml_ugli2_pcsadded_16.indi.blp", header=FALSE)
ugli2_resids_education <- fread("greml_ugli2_pcsadded_17.indi.blp", header=FALSE)

#######################################
## B. Only keep IID and resid column ##
#######################################
cat("checking headers given")
print(head(cyto_resids_bmi),10)

cyto_resids_VA_rt <- cyto_resids_VA_rt[,c('V2','V6')]
cyto_resids_PS_rt <- cyto_resids_PS_rt[,c('V2','V6')] 
cyto_resids_EM_acc <- cyto_resids_EM_acc[,c('V2','V6')] 
cyto_resids_WM_acc <- cyto_resids_WM_acc[,c('V2','V6')] 
cyto_resids_WM_rt <- cyto_resids_WM_rt[,c('V2','V6')] 
cyto_resids_RFFT <- cyto_resids_RFFT[,c('V2','V6')]
cyto_resids_panas_neg <- cyto_resids_panas_neg[,c('V2','V6')] 
cyto_resids_panas_pos <- cyto_resids_panas_pos[,c('V2','V6')] 
cyto_resids_minidep <- cyto_resids_minidep[,c('V2','V6')] 
cyto_resids_minianx <- cyto_resids_minianx[,c('V2','V6')] 
cyto_resids_minimdd <- cyto_resids_minimdd[,c('V2','V6')] 
cyto_resids_minigad <- cyto_resids_minigad[,c('V2','V6')] 
cyto_resids_hscrp <- cyto_resids_hscrp[,c('V2','V6')] 
cyto_resids_bmi <- cyto_resids_bmi[,c('V2','V6')] 
cyto_resids_alcohol <- cyto_resids_alcohol[,c('V2','V6')] 
cyto_resids_smoke <- cyto_resids_smoke[,c('V2','V6')] 
cyto_resids_education <- cyto_resids_education[,c('V2','V6')] 

gsa_resids_VA_rt <- gsa_resids_VA_rt[,c('V2','V6')]
gsa_resids_PS_rt <- gsa_resids_PS_rt[,c('V2','V6')] 
gsa_resids_EM_acc <- gsa_resids_EM_acc[,c('V2','V6')] 
gsa_resids_WM_acc <- gsa_resids_WM_acc[,c('V2','V6')] 
gsa_resids_WM_rt <- gsa_resids_WM_rt[,c('V2','V6')] 
gsa_resids_RFFT <- gsa_resids_RFFT[,c('V2','V6')]
gsa_resids_panas_neg <- gsa_resids_panas_neg[,c('V2','V6')] 
gsa_resids_panas_pos <- gsa_resids_panas_pos[,c('V2','V6')] 
gsa_resids_minidep <- gsa_resids_minidep[,c('V2','V6')] 
gsa_resids_minianx <- gsa_resids_minianx[,c('V2','V6')] 
gsa_resids_minimdd <- gsa_resids_minimdd[,c('V2','V6')] 
gsa_resids_minigad <- gsa_resids_minigad[,c('V2','V6')] 
gsa_resids_hscrp <- gsa_resids_hscrp[,c('V2','V6')] 
gsa_resids_bmi <- gsa_resids_bmi[,c('V2','V6')] 
gsa_resids_alcohol <- gsa_resids_alcohol[,c('V2','V6')] 
gsa_resids_smoke <- gsa_resids_smoke[,c('V2','V6')] 
gsa_resids_education <- gsa_resids_education[,c('V2','V6')] 

ugli2_resids_VA_rt <- ugli2_resids_VA_rt[,c('V2','V6')]
ugli2_resids_PS_rt <- ugli2_resids_PS_rt[,c('V2','V6')] 
ugli2_resids_EM_acc <- ugli2_resids_EM_acc[,c('V2','V6')] 
ugli2_resids_WM_acc <- ugli2_resids_WM_acc[,c('V2','V6')] 
ugli2_resids_WM_rt <- ugli2_resids_WM_rt[,c('V2','V6')] 
ugli2_resids_RFFT <- ugli2_resids_RFFT[,c('V2','V6')]
ugli2_resids_panas_neg <- ugli2_resids_panas_neg[,c('V2','V6')] 
ugli2_resids_panas_pos <- ugli2_resids_panas_pos[,c('V2','V6')] 
ugli2_resids_minidep <- ugli2_resids_minidep[,c('V2','V6')] 
ugli2_resids_minianx <- ugli2_resids_minianx[,c('V2','V6')] 
ugli2_resids_minimdd <- ugli2_resids_minimdd[,c('V2','V6')] 
ugli2_resids_minigad <- ugli2_resids_minigad[,c('V2','V6')] 
ugli2_resids_hscrp <- ugli2_resids_hscrp[,c('V2','V6')] 
ugli2_resids_bmi <- ugli2_resids_bmi[,c('V2','V6')] 
ugli2_resids_alcohol <- ugli2_resids_alcohol[,c('V2','V6')] 
ugli2_resids_smoke <- ugli2_resids_smoke[,c('V2','V6')] 
ugli2_resids_education <- ugli2_resids_education[,c('V2','V6')] 


#################################################
## C. Rename column and merge within each chip ##
#################################################
cyto_resids_VA_rt <- rename(cyto_resids_VA_rt, IID=V2,VA_RES_sd_res=V6)
cyto_resids_PS_rt <-rename(cyto_resids_PS_rt, IID=V2,PS_RES_sd_res=V6)
cyto_resids_EM_acc <- rename(cyto_resids_EM_acc, IID=V2, EM_ACC_sd_res=V6)
cyto_resids_WM_acc <-rename(cyto_resids_WM_acc, IID=V2, WM_ACC_sd_res=V6)
cyto_resids_WM_rt <- rename(cyto_resids_WM_rt, IID=V2, WM_RES_sd_res=V6)
cyto_resids_RFFT <- rename(cyto_resids_RFFT, IID=V2, RFFT_final_sd_res=V6)
cyto_resids_panas_neg <- rename(cyto_resids_panas_neg, IID=V2, panas_neg_total_final_sd_res=V6) 
cyto_resids_panas_pos <- rename(cyto_resids_panas_pos, IID=V2, panas_pos_total_final_sd_res=V6) 
cyto_resids_minidep <- rename(cyto_resids_minidep, IID=V2, minidep_1av1_res=V6)  
cyto_resids_minianx <- rename(cyto_resids_minianx, IID=V2, minianx_1av1_res=V6)  
cyto_resids_minimdd <- rename(cyto_resids_minimdd, IID=V2, minimdd_1av1_res=V6)   
cyto_resids_minigad <- rename(cyto_resids_minigad, IID=V2, minigad_1av1_res=V6) 
cyto_resids_hscrp <- rename(cyto_resids_hscrp, IID=V2, hsCRP_log_res=V6)
cyto_resids_bmi <- rename(cyto_resids_bmi, IID=V2, BMI_1av1_res=V6) 
cyto_resids_alcohol <- rename(cyto_resids_alcohol, IID=V2, alcohol_freq_1a_q_2_res=V6)  
cyto_resids_smoke <- rename(cyto_resids_smoke, IID=V2, Currentsmok_1a_q_2_res=V6)
cyto_resids_education <- rename(cyto_resids_education, IID=V2, education_res=V6)

gsa_resids_VA_rt <- rename(gsa_resids_VA_rt, IID=V2,VA_RES_sd_res=V6)
gsa_resids_PS_rt <-rename(gsa_resids_PS_rt, IID=V2,PS_RES_sd_res=V6)
gsa_resids_EM_acc <- rename(gsa_resids_EM_acc, IID=V2, EM_ACC_sd_res=V6)
gsa_resids_WM_acc <-rename(gsa_resids_WM_acc, IID=V2, WM_ACC_sd_res=V6)
gsa_resids_WM_rt <- rename(gsa_resids_WM_rt, IID=V2, WM_RES_sd_res=V6)
gsa_resids_RFFT <- rename(gsa_resids_RFFT, IID=V2, RFFT_final_sd_res=V6)
gsa_resids_panas_neg <- rename(gsa_resids_panas_neg, IID=V2, panas_neg_total_final_sd_res=V6) 
gsa_resids_panas_pos <- rename(gsa_resids_panas_pos, IID=V2, panas_pos_total_final_sd_res=V6) 
gsa_resids_minidep <- rename(gsa_resids_minidep, IID=V2, minidep_1av1_res=V6)  
gsa_resids_minianx <- rename(gsa_resids_minianx, IID=V2, minianx_1av1_res=V6)  
gsa_resids_minimdd <- rename(gsa_resids_minimdd, IID=V2, minimdd_1av1_res=V6)   
gsa_resids_minigad <- rename(gsa_resids_minigad, IID=V2, minigad_1av1_res=V6) 
gsa_resids_hscrp <- rename(gsa_resids_hscrp, IID=V2, hsCRP_log_res=V6)
gsa_resids_bmi <- rename(gsa_resids_bmi, IID=V2, BMI_1av1_res=V6) 
gsa_resids_alcohol <- rename(gsa_resids_alcohol, IID=V2, alcohol_freq_1a_q_2_res=V6)  
gsa_resids_smoke <- rename(gsa_resids_smoke, IID=V2, Currentsmok_1a_q_2_res=V6)
gsa_resids_education <- rename(gsa_resids_education, IID=V2, education_res=V6)

ugli2_resids_VA_rt <- rename(ugli2_resids_VA_rt, IID=V2,VA_RES_sd_res=V6)
ugli2_resids_PS_rt <-rename(ugli2_resids_PS_rt, IID=V2,PS_RES_sd_res=V6)
ugli2_resids_EM_acc <- rename(ugli2_resids_EM_acc, IID=V2, EM_ACC_sd_res=V6)
ugli2_resids_WM_acc <-rename(ugli2_resids_WM_acc, IID=V2, WM_ACC_sd_res=V6)
ugli2_resids_WM_rt <- rename(ugli2_resids_WM_rt, IID=V2, WM_RES_sd_res=V6)
ugli2_resids_RFFT <- rename(ugli2_resids_RFFT, IID=V2, RFFT_final_sd_res=V6)
ugli2_resids_panas_neg <- rename(ugli2_resids_panas_neg, IID=V2, panas_neg_total_final_sd_res=V6) 
ugli2_resids_panas_pos <- rename(ugli2_resids_panas_pos, IID=V2, panas_pos_total_final_sd_res=V6) 
ugli2_resids_minidep <- rename(ugli2_resids_minidep, IID=V2, minidep_1av1_res=V6)  
ugli2_resids_minianx <- rename(ugli2_resids_minianx, IID=V2, minianx_1av1_res=V6)  
ugli2_resids_minimdd <- rename(ugli2_resids_minimdd, IID=V2, minimdd_1av1_res=V6)   
ugli2_resids_minigad <- rename(ugli2_resids_minigad, IID=V2, minigad_1av1_res=V6) 
ugli2_resids_hscrp <- rename(ugli2_resids_hscrp, IID=V2, hsCRP_log_res=V6)
ugli2_resids_bmi <- rename(ugli2_resids_bmi, IID=V2, BMI_1av1_res=V6) 
ugli2_resids_alcohol <- rename(ugli2_resids_alcohol, IID=V2, alcohol_freq_1a_q_2_res=V6)  
ugli2_resids_smoke <- rename(ugli2_resids_smoke, IID=V2, Currentsmok_1a_q_2_res=V6)
ugli2_resids_education <- rename(ugli2_resids_education, IID=V2, education_res=V6)

################################################
## E. Merge data (residuals) within each chip ##
################################################
cyto_list <- list(cyto_resids_VA_rt,cyto_resids_PS_rt,cyto_resids_EM_acc,cyto_resids_WM_acc,cyto_resids_WM_rt,cyto_resids_RFFT,cyto_resids_panas_neg,cyto_resids_panas_pos,cyto_resids_minidep,cyto_resids_minianx,cyto_resids_minimdd,cyto_resids_minigad,cyto_resids_hscrp,cyto_resids_bmi,cyto_resids_alcohol,cyto_resids_smoke,cyto_resids_education)
cyto_resids_merged <- Reduce(function(x,y) merge(x,y, all=TRUE),cyto_list)

gsa_list <- list(gsa_resids_VA_rt,gsa_resids_PS_rt,gsa_resids_EM_acc,gsa_resids_WM_acc,gsa_resids_WM_rt,gsa_resids_RFFT,gsa_resids_panas_neg,gsa_resids_panas_pos,gsa_resids_minidep,gsa_resids_minianx,gsa_resids_minimdd,gsa_resids_minigad,gsa_resids_hscrp,gsa_resids_bmi,gsa_resids_alcohol,gsa_resids_smoke,gsa_resids_education)
gsa_resids_merged <- Reduce(function(x,y) merge(x,y, all=TRUE),gsa_list)

ugli2_list <- list(ugli2_resids_VA_rt,ugli2_resids_PS_rt,ugli2_resids_EM_acc,ugli2_resids_WM_acc,ugli2_resids_WM_rt,ugli2_resids_RFFT,ugli2_resids_panas_neg,ugli2_resids_panas_pos,ugli2_resids_minidep,ugli2_resids_minianx,ugli2_resids_minimdd,ugli2_resids_minigad,ugli2_resids_hscrp,ugli2_resids_bmi,ugli2_resids_alcohol,ugli2_resids_smoke,ugli2_resids_education)
ugli2_resids_merged <- Reduce(function(x,y) merge(x,y, all=TRUE),ugli2_list)

print(head(cyto_resids_merged),20)
print(head(gsa_resids_merged),20)
print(head(ugli2_resids_merged),20)
cat("cyto resid dim")
print(dim(cyto_resids_merged))
cat("gsa resid dim")
print(dim(gsa_resids_merged))
cat("ugli2 resid dim")
print(dim(ugli2_resids_merged))

##########################################
## E. Save file of resids for each chip ##
##########################################
write.table(cyto_resids_merged, "cyto_resids_merged.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(gsa_resids_merged, "gsa_resids_merged.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(ugli2_resids_merged, "ugli2_resids_merged.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

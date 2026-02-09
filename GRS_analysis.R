# Project: GRS for inflammation markers on outcomes in Lifelines (cognition, depression, anxiety phenotypes)
# Chloe Slaney

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
library(ivreg)

##########################################################
#################### PART 1: Load data ###################
##########################################################
rm(list = ls())

###############################
## 1. Cleaned phenotype data ##
###############################
lifelines_phenotype <- fread("Master_Database_FINAL_CLEAN.csv", header=TRUE)

#Standardize phenotype data
lifelines_phenotype$VA_RES_sd <- scale(lifelines_phenotype$VA_RES)
lifelines_phenotype$PS_RES_sd <- scale(lifelines_phenotype$PS_RES)
lifelines_phenotype$EM_ACC_sd <- scale(lifelines_phenotype$EM_ACC)
lifelines_phenotype$WM_RES_sd <- scale(lifelines_phenotype$WM_RES)
lifelines_phenotype$WM_ACC_sd <- scale(lifelines_phenotype$WM_ACC)
lifelines_phenotype$RFFT_final_sd <- scale(lifelines_phenotype$RFFT_final)
lifelines_phenotype$panas_neg_total_final_sd <- scale(lifelines_phenotype$panas_neg_total_final)
lifelines_phenotype$panas_pos_total_final_sd <- scale(lifelines_phenotype$panas_pos_total_final)
lifelines_phenotype$hsCRP_log_sd <- scale(lifelines_phenotype$hsCRP_log)
#cat("check standardised vars")
#describe(lifelines_phenotype$VA_RES_sd)
#summary(lifelines_phenotype$VA_RES_sd)
#sum(! is.na(lifelines_phenotype$VA_RES_sd))

cat("health summed score details")
describe(lifelines_phenotype$health_excluded_summedscore)
sum(! is.na(lifelines_phenotype$health_excluded_summedscore))
summary(lifelines_phenotype$health_excluded_summedscore)

#Create factor
lifelines_phenotype$education <- as.factor(lifelines_phenotype$education_attain_q_1.x)
lifelines_phenotype$GENDER_1a_v_1 <- as.factor(lifelines_phenotype$GENDER_1a_v_1)

############################
## 2. Genetic risk scores ##
############################
### CytoSNP ###
ahluwalia_cis_cyto <- fread("ahluwaliacis_score_sum_cyto.profile", header=TRUE)
ahluwalia_gw_cyto <- fread("ahluwalia_genomewide_score_sum_cyto.profile", header=TRUE)
borges_proxyadd_cyto <- fread("borges_proxiesadded_score_sum_cyto.profile", header=TRUE)
rosa_cyto <- fread("rosa_score_sum_cyto.profile", header=TRUE)
said_cis_proxyadd_cyto <- fread("said_cis_proxiesadded_score_sum_cyto.profile", header=TRUE)
said_gw_proxyadd_cyto <- fread("said_genomewide_proxiesadded_score_sum_cyto.profile", header=TRUE)
sarwar_cyto <- fread("sarwar_score_sum_cyto.profile", header=TRUE)

##Renamed
ahluwalia_cis_cyto_prs <- rename(ahluwalia_cis_cyto, ahluwalia_cis_cyto=SCORESUM)
AHLUWALIA_CIS_CYTO_PRS <- ahluwalia_cis_cyto_prs[, c('IID', 'ahluwalia_cis_cyto')]
ahluwalia_gw_cyto_prs <- rename(ahluwalia_gw_cyto, ahluwalia_gw_cyto=SCORESUM)
AHLUWALIA_GW_CYTO_PRS <- ahluwalia_gw_cyto_prs[, c('IID', 'ahluwalia_gw_cyto')]
borges_proxyadd_cyto_prs <- rename(borges_proxyadd_cyto, borges_proxyadd_cyto=SCORESUM)
BORGES_PROXYADD_CYTO_PRS <- borges_proxyadd_cyto_prs[, c('IID', 'borges_proxyadd_cyto')]
rosa_cyto_prs <- rename(rosa_cyto, rosa_cyto=SCORESUM)
ROSA_CYTO_PRS <- rosa_cyto_prs[, c('IID', 'rosa_cyto')]
said_cis_proxyadd_cyto_prs <- rename(said_cis_proxyadd_cyto, said_cis_proxyadd_cyto=SCORESUM)
SAID_CIS_PROXYADD_CYTO_PRS <- said_cis_proxyadd_cyto_prs[, c('IID', 'said_cis_proxyadd_cyto')]
said_gw_proxyadd_cyto_prs <- rename(said_gw_proxyadd_cyto, said_gw_proxyadd_cyto=SCORESUM)
SAID_GW_PROXYADD_CYTO_PRS <- said_gw_proxyadd_cyto_prs[, c('IID', 'said_gw_proxyadd_cyto')]
sarwar_cyto_prs <- rename(sarwar_cyto, sarwar_cyto=SCORESUM)
SARWAR_CYTO_PRS <- sarwar_cyto_prs[, c('IID', 'sarwar_cyto')]

### GSA ###
ahluwalia_cis_gsa <- fread("ahluwaliacis_score_sum_gsa.profile", header=TRUE)
ahluwalia_gw_gsa <- fread("ahluwalia_genomewide_score_sum_gsa.profile", header=TRUE)
borges_proxyadd_gsa <- fread("borges_proxiesadded_score_sum_gsa.profile", header=TRUE)
rosa_gsa <- fread("rosa_score_sum_gsa.profile", header=TRUE)
said_cis_proxyadd_gsa <- fread("said_cis_proxiesadded_score_sum_gsa.profile", header=TRUE)
said_gw_proxyadd_gsa <- fread("said_genomewide_proxiesadded_score_sum_gsa.profile", header=TRUE)
sarwar_gsa <- fread("sarwar_score_sum_gsa.profile", header=TRUE)

##Renamed
ahluwalia_cis_gsa_prs <- rename(ahluwalia_cis_gsa, ahluwalia_cis_gsa=SCORESUM)
AHLUWALIA_CIS_GSA_PRS <- ahluwalia_cis_gsa_prs[, c('IID', 'ahluwalia_cis_gsa')]
ahluwalia_gw_gsa_prs <- rename(ahluwalia_gw_gsa, ahluwalia_gw_gsa=SCORESUM)
AHLUWALIA_GW_GSA_PRS <- ahluwalia_gw_gsa_prs[, c('IID', 'ahluwalia_gw_gsa')]
borges_proxyadd_gsa_prs <- rename(borges_proxyadd_gsa, borges_proxyadd_gsa=SCORESUM)
BORGES_PROXYADD_GSA_PRS <- borges_proxyadd_gsa_prs[, c('IID', 'borges_proxyadd_gsa')]
rosa_gsa_prs <- rename(rosa_gsa, rosa_gsa=SCORESUM)
ROSA_GSA_PRS <- rosa_gsa_prs[, c('IID', 'rosa_gsa')]
said_cis_proxyadd_gsa_prs <- rename(said_cis_proxyadd_gsa, said_cis_proxyadd_gsa=SCORESUM)
SAID_CIS_PROXYADD_GSA_PRS <- said_cis_proxyadd_gsa_prs[, c('IID', 'said_cis_proxyadd_gsa')]
said_gw_proxyadd_gsa_prs <- rename(said_gw_proxyadd_gsa, said_gw_proxyadd_gsa=SCORESUM)
SAID_GW_PROXYADD_GSA_PRS <- said_gw_proxyadd_gsa_prs[, c('IID', 'said_gw_proxyadd_gsa')]
sarwar_gsa_prs <- rename(sarwar_gsa, sarwar_gsa=SCORESUM)
SARWAR_GSA_PRS <- sarwar_gsa_prs[, c('IID', 'sarwar_gsa')]

### Affymetrix ###
ahluwalia_cis_ugli2 <- fread("ahluwaliacis_score_sum_ugli2.profile", header=TRUE)
ahluwalia_gw_ugli2 <- fread("ahluwalia_genomewide_score_sum_ugli2.profile", header=TRUE)
borges_proxyadd_ugli2 <- fread("borges_proxiesadded_score_sum_ugli2.profile", header=TRUE)
rosa_ugli2 <- fread("rosa_score_sum_ugli2.profile", header=TRUE)
said_cis_proxyadd_ugli2 <- fread("said_cis_proxiesadded_score_sum_ugli2.profile", header=TRUE)
said_gw_proxyadd_ugli2 <- fread("said_genomewide_proxiesadded_score_sum_ugli2.profile", header=TRUE)
sarwar_ugli2 <- fread("sarwar_score_sum_ugli2.profile", header=TRUE)

##Renamed
ahluwalia_cis_ugli2_prs <- rename(ahluwalia_cis_ugli2, ahluwalia_cis_ugli2=SCORESUM)
AHLUWALIA_CIS_UGLI2_PRS <- ahluwalia_cis_ugli2_prs[, c('IID', 'ahluwalia_cis_ugli2')]
ahluwalia_gw_ugli2_prs <- rename(ahluwalia_gw_ugli2, ahluwalia_gw_ugli2=SCORESUM)
AHLUWALIA_GW_UGLI2_PRS <- ahluwalia_gw_ugli2_prs[, c('IID', 'ahluwalia_gw_ugli2')]
borges_proxyadd_ugli2_prs <- rename(borges_proxyadd_ugli2, borges_proxyadd_ugli2=SCORESUM)
BORGES_PROXYADD_UGLI2_PRS <- borges_proxyadd_ugli2_prs[, c('IID', 'borges_proxyadd_ugli2')]
rosa_ugli2_prs <- rename(rosa_ugli2, rosa_ugli2=SCORESUM)
ROSA_UGLI2_PRS <- rosa_ugli2_prs[, c('IID', 'rosa_ugli2')]
said_cis_proxyadd_ugli2_prs <- rename(said_cis_proxyadd_ugli2, said_cis_proxyadd_ugli2=SCORESUM)
SAID_CIS_PROXYADD_UGLI2_PRS <- said_cis_proxyadd_ugli2_prs[, c('IID', 'said_cis_proxyadd_ugli2')]
said_gw_proxyadd_ugli2_prs <- rename(said_gw_proxyadd_ugli2, said_gw_proxyadd_ugli2=SCORESUM)
SAID_GW_PROXYADD_UGLI2_PRS <- said_gw_proxyadd_ugli2_prs[, c('IID', 'said_gw_proxyadd_ugli2')]
sarwar_ugli2_prs <- rename(sarwar_ugli2, sarwar_ugli2=SCORESUM)
SARWAR_UGLI2_PRS <- sarwar_ugli2_prs[, c('IID', 'sarwar_ugli2')]

################################################################
## 3. Linkage files - link ID in pheno data with genetic data ##
################################################################
cytosnp_linkages <- fread("cytosnp_linkage_file_project_pseudo_id.txt", header=TRUE)
gsa_linkages <- fread("gsa_linkage_file_project_pseudo_id.txt", header=TRUE)
ugli2_linkages <-fread("UGLI2_samples_afterQC.txt", header=TRUE)
cytosnp_linkages <- rename(cytosnp_linkages, IID=cytosnp_ID)
gsa_linkages <- rename(gsa_linkages, IID=UGLI_ID)
ugli2_linkages <-rename(ugli2_linkages,IID=Barcode)

#############################
## 4. Principal components ##
#############################
## Used for grammar method.
PC_cyto <- fread("LL_CytoSNP_PCs.txt", header=TRUE)
PC_gsa <- fread("PCA_eur.UGLI.eigenvec")
PC_ugli2 <- fread("UGLI2_PCs.txt", header=TRUE)
#cat("UGLI PC file V1 and V2 columns identical = ", identical(PC_gsa$V1, PC_gsa$V2))

PC_cyto <- PC_cyto[,c('IID','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')]
PC_gsa <- rename(PC_gsa, IID=V2,PC1=V3, PC2=V4, PC3=V5, PC4=V6, PC5=V7, PC6=V8, PC7=V9, PC8=V10, PC9=V11, PC10=V12)
PC_gsa <- PC_gsa[,c('IID','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')]
PC_ugli2 <- PC_ugli2[,c('IID','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')]

#################################################################
#################### PART TWO: CytoSNP Prep  ####################
#################################################################
#################
## Merge files ##
#################
#linkage ID, genetic data merged with phenotype data (n=14,942)
lifelines_cyto_merged <- merge(lifelines_phenotype, cytosnp_linkages, by = "PROJECT_PSEUDO_ID") 
lifelines_cyto_merged <- merge(lifelines_cyto_merged, SAID_CIS_PROXYADD_CYTO_PRS, by = "IID") 
lifelines_cyto_merged <- merge(lifelines_cyto_merged, SAID_GW_PROXYADD_CYTO_PRS, by = "IID") 
lifelines_cyto_merged <- merge(lifelines_cyto_merged, BORGES_PROXYADD_CYTO_PRS, by = "IID") 
lifelines_cyto_merged <- merge(lifelines_cyto_merged, AHLUWALIA_CIS_CYTO_PRS, by = "IID") 
lifelines_cyto_merged <- merge(lifelines_cyto_merged, AHLUWALIA_GW_CYTO_PRS, by = "IID") 
lifelines_cyto_merged <- merge(lifelines_cyto_merged, ROSA_CYTO_PRS, by = "IID") 
lifelines_cyto_merged <- merge(lifelines_cyto_merged, SARWAR_CYTO_PRS, by = "IID") 
lifelines_cyto_merged <- merge(lifelines_cyto_merged, PC_cyto, by = "IID")
lifelines_cyto_merged$chip <- 1
#print(head(lifelines_cyto_merged),5)
#cat("number of rows in cyto + pheno merged = ")
#print(nrow(lifelines_cyto_merged))

#####################
## Standardise GRS ##
#####################
lifelines_cyto_merged$ahluwalia_cis_cyt_sd <- scale(lifelines_cyto_merged$ahluwalia_cis_cyto)
lifelines_cyto_merged$ahluwalia_gw_cyt_sd <- scale(lifelines_cyto_merged$ahluwalia_gw_cyto)
lifelines_cyto_merged$borges_proxyadd_cyt_sd <- scale(lifelines_cyto_merged$borges_proxyadd_cyto)
lifelines_cyto_merged$said_cis_proxyadd_cyt_sd <- scale(lifelines_cyto_merged$said_cis_proxyadd_cyto)
lifelines_cyto_merged$said_gw_proxyadd_cyt_sd <- scale(lifelines_cyto_merged$said_gw_proxyadd_cyto)
lifelines_cyto_merged$rosa_cyt_sd <- scale(lifelines_cyto_merged$rosa_cyto)
lifelines_cyto_merged$sarwar_cyt_sd <- scale(lifelines_cyto_merged$sarwar_cyto)

## NOTE: log CRP has negative values - as some CRP values are near 0
hist(lifelines_cyto_merged$said_cis_proxyadd_cyt_sd)
hist(lifelines_cyto_merged$said_gw_proxyadd_cyt_sd) 
hist(lifelines_cyto_merged$borges_proxyadd_cyt_sd) 
hist(lifelines_cyto_merged$ahluwalia_cis_cyt_sd) 
hist(lifelines_cyto_merged$ahluwalia_gw_cyt_sd) 
hist(lifelines_cyto_merged$rosa_cyt_sd) 
hist(lifelines_cyto_merged$sarwar_cyt_sd) 

#cat("check standardised prs")
#describe(lifelines_cyto_merged$ahluwalia_cis_cyt_sd)
#describe(lifelines_cyto_merged$ahluwalia_gw_cyt_sd)

###################################################################################
## Remove duplicates and 1st-degree relatives in CytoSNP who are in  GSA + UGLI2 ##
###################################################################################
## Due to poorer imputation quality of CytoSNP, removed individuals from this chip where possible.
## NOTE: old list of duplics and 1st-degree relatives in CytoSNP for UGLI had 29 more people excluded - unclear why.
#cat("number cytosnp + pheno = ") #14,942
#sum(! is.na(lifelines_cyto_merged$IID))

# List of duplicates and 1st-degree relatives in UGLI and UGLI2
cytosnp_duplicates_and_1stdegree_UGLI <- fread("CytoSNP_duplicates+1stdgr_in_UGLI.txt", header=FALSE)
cytosnp_duplicates_and_1stdegree_UGLI2 <- fread("CytoSNP_duplicates+1stdgr_in_UGLI2.txt", header=FALSE)
cytosnp_duplicates_and_1stdegree_UGLI <- rename(cytosnp_duplicates_and_1stdegree_UGLI, IID=V1)
cytosnp_duplicates_and_1stdegree_UGLI2 <- rename(cytosnp_duplicates_and_1stdegree_UGLI2, IID=V1)
cytosnp_duplicates_and_1stdegree_UGLI <- cytosnp_duplicates_and_1stdegree_UGLI[, c('IID')]
cytosnp_duplicates_and_1stdegree_UGLI2 <- cytosnp_duplicates_and_1stdegree_UGLI2[, c('IID')]

# Sanity checks - UGLI file all distinct, but UGLI2 file has some individuals not distinct in list
#cat("distinct in cyto-ugli list")
#cyto_distinct_ugli <- distinct(cytosnp_duplicates_and_1stdegree_UGLI)
#sum(! is.na(cyto_distinct_ugli))
#cat("distinct in cyto-ugli2 list")
#cyto_distinct_ugli2 <- distinct(cytosnp_duplicates_and_1stdegree_UGLI2)
#sum(! is.na(cyto_distinct_ugli2))

# Remove duplicates and 1st-degree relatives who are in UGLI and UGLI2
cat("Sample size after UGLI list removed") 
lifelines_cyto_merged <- lifelines_cyto_merged[! lifelines_cyto_merged$IID %in% cytosnp_duplicates_and_1stdegree_UGLI$IID, ]
sum(! is.na (lifelines_cyto_merged$IID))
#cat("How many unique UGLI2 need to be removed") 
#check <- merge(lifelines_cyto_merged, cyto_distinct_ugli2, by="IID")
#sum(! is.na(check$IID))
cat("Sample size after UGLI2 additionally removed")
lifelines_cyto_merged <- lifelines_cyto_merged[! lifelines_cyto_merged$IID %in% cytosnp_duplicates_and_1stdegree_UGLI2$IID, ]
sum(! is.na(lifelines_cyto_merged$IID))

#Check - no duplicates or 1st-degree relatives remain = yes.
#cyto_ugli_removed <- merge(lifelines_cyto_merged, cytosnp_duplicates_and_1stdegree_UGLI, by="IID")
#cat("number cyto with ugli removed")
#sum(! is.na(cyto_ugli_removed$IID)) #0
#cyto_ugli2_removed <-merge(lifelines_cyto_merged, cytosnp_duplicates_and_1stdegree_UGLI2, by="IID")
#cat("number cyto with ugli2 removed") #0
#sum(! is.na(cyto_ugli2_removed$IID))

################################
## Prepare for GRAMMAR Method ##
################################
### Step 1: Extract files for REML (gcta to extract residuals)
#Phenotype file (contains FID, IID, outcomes, covariates, crp)
#Note composite cog not found and not included.
#Get FID (differs in cyto) as needed for running reml
cyto_grm_ids <- fread("GRM_CYTO_GCTA_SPARSE.grm.id", header=FALSE)
cyto_grm_ids <- rename(cyto_grm_ids, FID=V1, IID=V2)
lifelines_cyto_merged <- merge(lifelines_cyto_merged, cyto_grm_ids, by="IID")

#pheno file
#cyto_merged_phenotypes_gcta <- lifelines_cyto_merged[,c("FID","IID","VA_RES_sd","PS_RES_sd","EM_ACC_sd","WM_ACC_sd","WM_RES_sd", "RFFT_final_sd","panas_neg_total_final_sd","panas_pos_total_final_sd","minidep_1av1","minianx_1av1","minimdd_1av1","minigad_1av1","hsCRP_log","BMI_1av1","alcohol_freq_1a_q_2","Currentsmok_1a_q_2","education")]
#print(head(cyto_merged_phenotypes_gcta),50)
#write.table(cyto_merged_phenotypes_gcta, "cyto_merged_pheno_for_gcta_FIDchanged.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

#covariate pc file
#cyto_merged_pcs_gcta <- lifelines_cyto_merged[,c("FID","IID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")]
#write.table(cyto_merged_pcs_gcta, "cyto_merged_pcs_for_gcta_FIDchanged.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

### Step 2: Load in residuals to be used as outcomes in adjusted analysis and merge with data
cyto_resids <- fread("cyto_resids_merged.txt", header=TRUE)
lifelines_cyto_merged <- merge(lifelines_cyto_merged, cyto_resids, by="IID") 
#cat("lifelines cyto and resids merged")
#print(head(lifelines_cyto_merged),20)
#cat("lifelines cyto + resids dims")
#print(dim(lifelines_cyto_merged))

################################################
## Create file for REGENIE - another analysis ##
################################################
cyto_regenie_pheno_covar <- lifelines_cyto_merged[,c("FID", "IID", "VA_RES_sd","PS_RES_sd", "EM_ACC_sd", "WM_ACC_sd", "RFFT_final_sd", "panas_neg_total_final_sd", "panas_pos_total_final_sd", "AGE_1a_q_1", "GENDER_1a_v_1", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")]
#write.table(cyto_regenie_pheno_covar, "cyto_regenie_pheno_covar.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
#summary(cyto_regenie_pheno_covar)
#print(head(cyto_regenie_pheno_covar))

#########################################################
#################### PART 3: GSA prep ###################
#########################################################
##Upload list of non-European individuals
non_euro_gsa <- fread("ugli_qc_release1_samples.csv",header=TRUE)
non_euro_gsa <- rename(non_euro_gsa, genotyping_name=Sample_ID)

#################
## Merge files ##
#################
#Linkage ID, genetic data merged with phenotype data (n=31,810)
lifelines_gsa_merged <- merge(lifelines_phenotype, gsa_linkages, by = "PROJECT_PSEUDO_ID") 
lifelines_gsa_merged <- merge(lifelines_gsa_merged, SAID_CIS_PROXYADD_GSA_PRS, by = "IID") 
lifelines_gsa_merged <- merge(lifelines_gsa_merged, SAID_GW_PROXYADD_GSA_PRS, by = "IID") 
lifelines_gsa_merged <- merge(lifelines_gsa_merged, BORGES_PROXYADD_GSA_PRS, by = "IID") 
lifelines_gsa_merged <- merge(lifelines_gsa_merged, AHLUWALIA_CIS_GSA_PRS, by = "IID") 
lifelines_gsa_merged <- merge(lifelines_gsa_merged, AHLUWALIA_GW_GSA_PRS, by = "IID") 
lifelines_gsa_merged <- merge(lifelines_gsa_merged, ROSA_GSA_PRS, by = "IID") 
lifelines_gsa_merged <- merge(lifelines_gsa_merged, SARWAR_GSA_PRS, by = "IID") 
lifelines_gsa_merged <- merge(lifelines_gsa_merged, PC_gsa, by = "IID")
lifelines_gsa_merged <- merge(lifelines_gsa_merged, non_euro_gsa, by = "genotyping_name")
lifelines_gsa_merged$chip <- 2

##############################################################
## Non-european individuals are not included - checked here ##
##############################################################
cat("number of rows in gsa + pheno merged = ") #31,810
print(nrow(lifelines_gsa_merged))

cat("number non-euro in gsa merged file")
non_euro_gsa_merged <- subset(lifelines_gsa_merged, lifelines_gsa_merged$PCA_european=="FALSE")
print(nrow(non_euro_gsa_merged))

#cat("number euro in gsa merged file")
#euro_gsa_merged <- subset(lifelines_gsa_merged, lifelines_gsa_merged$PCA_european=="TRUE")
#print(nrow(euro_gsa_merged))

#check n non-euro in genetics data (not merged file) #34
#lifelines_gsa_genetic <- merge(gsa_linkages, SAID_CIS_PROXYADD_GSA_PRS, by="IID")
#lifelines_gsa_genetic <- merge(lifelines_gsa_genetic, non_euro_gsa, by = "genotyping_name")
#cat("number non-euro in genetic gsa file")
#genetic_gsa_noneuro <- subset(lifelines_gsa_genetic, lifelines_gsa_genetic$PCA_european=="FALSE")
#print(nrow(genetic_gsa_noneuro))

#################################################################################
## Remove duplicates and 1st-degree relatives in GSA who are present in  UGLI2 ##
#################################################################################
cat("number GSA + pheno = ") #n=26,334
sum(! is.na(lifelines_gsa_merged$IID))

# List of duplicates and 1st-degree relatives in UGLI2
GSA_duplicates_and_1stdegree_UGLI2 <- fread("UGLI_duplicates+1stdgr_in_UGLI2.txt", header=FALSE)
GSA_duplicates_and_1stdegree_UGLI2 <- rename(GSA_duplicates_and_1stdegree_UGLI2, IID=V1)
GSA_duplicates_and_1stdegree_UGLI2 <- GSA_duplicates_and_1stdegree_UGLI2[, c('IID')]

# List of distinct individuals which must be removed 
cat("distinct in ugli-ugli2 list")
ugli_distinct_ugli2 <- distinct(GSA_duplicates_and_1stdegree_UGLI2)
sum(! is.na(ugli_distinct_ugli2))

# Remove duplicates and 1st-degree relatives who are in UGLI and UGLI2
lifelines_gsa_merged <- lifelines_gsa_merged[! lifelines_gsa_merged$IID %in% GSA_duplicates_and_1stdegree_UGLI2$IID, ]
sum(! is.na(lifelines_gsa_merged$IID))

#Check - no duplicates or 1st-degree relatives remain.
gsa_ugli2_removed <- merge(lifelines_gsa_merged, GSA_duplicates_and_1stdegree_UGLI2, by="IID")
cat("number gsa with ugli2 included")
sum(! is.na(gsa_ugli2_removed$IID))

#####################
## Standardise GRS ##
#####################
lifelines_gsa_merged$ahluwalia_cis_gsa_sd <- scale(lifelines_gsa_merged$ahluwalia_cis_gsa)
lifelines_gsa_merged$ahluwalia_gw_gsa_sd <- scale(lifelines_gsa_merged$ahluwalia_gw_gsa)
lifelines_gsa_merged$borges_proxyadd_gsa_sd <- scale(lifelines_gsa_merged$borges_proxyadd_gsa)
lifelines_gsa_merged$said_cis_proxyadd_gsa_sd <- scale(lifelines_gsa_merged$said_cis_proxyadd_gsa)
lifelines_gsa_merged$said_gw_proxyadd_gsa_sd <- scale(lifelines_gsa_merged$said_gw_proxyadd_gsa)
lifelines_gsa_merged$rosa_gsa_sd <- scale(lifelines_gsa_merged$rosa_gsa)
lifelines_gsa_merged$sarwar_gsa_sd <- scale(lifelines_gsa_merged$sarwar_gsa)

################################
## Prepare for GRAMMAR Method ##
################################
lifelines_gsa_merged$FID <- lifelines_gsa_merged$IID

### Step 1: Extract files for REML (gcta to extract residuals)
#pheno file
#gsa_merged_phenotypes_gcta <- lifelines_gsa_merged[,c("FID","IID","VA_RES_sd","PS_RES_sd","EM_ACC_sd","WM_ACC_sd","WM_RES_sd", "RFFT_final_sd","panas_neg_total_final_sd","panas_pos_total_final_sd","minidep_1av1","minianx_1av1","minimdd_1av1","minigad_1av1","hsCRP_log","BMI_1av1","alcohol_freq_1a_q_2","Currentsmok_1a_q_2","education")]
#print(head(gsa_merged_phenotypes_gcta),50)
#write.table(gsa_merged_phenotypes_gcta, "gsa_merged_pheno_for_gcta.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

#covariate pc file
#gsa_merged_pcs_gcta <- lifelines_gsa_merged[,c("FID","IID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")]
#write.table(gsa_merged_pcs_gcta, "gsa_merged_pcs_for_gcta.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

### Step 2: Load in residuals to be used as outcomes in adjusted analysis and merge with data
gsa_resids <- fread("gsa_resids_merged.txt", header=TRUE)
lifelines_gsa_merged <- merge(lifelines_gsa_merged, gsa_resids, by="IID") 
#cat("lifelines gsa and resids merged")
#print(head(lifelines_gsa_merged),20)
#cat("lifelines gda + resids dims")
#print(dim(lifelines_gsa_merged))

################################################
## Create file for REGENIE - another analysis ##
################################################
gsa_regenie_pheno_covar <- lifelines_gsa_merged[,c("FID", "IID", "VA_RES_sd","PS_RES_sd", "EM_ACC_sd", "WM_ACC_sd", "RFFT_final_sd", "panas_neg_total_final_sd", "panas_pos_total_final_sd", "AGE_1a_q_1", "GENDER_1a_v_1", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")]
#write.table(gsa_regenie_pheno_covar, "gsa_regenie_pheno_covar.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
#summary(gsa_regenie_pheno_covar)
#print(head(gsa_regenie_pheno_covar))

###############################################################
################### PART 4: Affymetrix Prep ##################
###############################################################
#################
## Merge files ##
#################
#Linkage ID, Genetic data merged with phenotype data (n=26,334)
lifelines_ugli2_merged <- merge(lifelines_phenotype, ugli2_linkages, by = "PROJECT_PSEUDO_ID") 
lifelines_ugli2_merged <- merge(lifelines_ugli2_merged, SAID_CIS_PROXYADD_UGLI2_PRS, by = "IID") 
lifelines_ugli2_merged <- merge(lifelines_ugli2_merged, SAID_GW_PROXYADD_UGLI2_PRS, by = "IID") 
lifelines_ugli2_merged <- merge(lifelines_ugli2_merged, BORGES_PROXYADD_UGLI2_PRS, by = "IID") 
lifelines_ugli2_merged <- merge(lifelines_ugli2_merged, AHLUWALIA_CIS_UGLI2_PRS, by = "IID") 
lifelines_ugli2_merged <- merge(lifelines_ugli2_merged, AHLUWALIA_GW_UGLI2_PRS, by = "IID") 
lifelines_ugli2_merged <- merge(lifelines_ugli2_merged, ROSA_UGLI2_PRS, by = "IID") 
lifelines_ugli2_merged <- merge(lifelines_ugli2_merged, SARWAR_UGLI2_PRS, by = "IID") 
lifelines_ugli2_merged <- merge(lifelines_ugli2_merged, PC_ugli2, by = "IID")
lifelines_ugli2_merged$chip <- 3
#print(head(lifelines_ugli2_merged),5)
#cat("number of rows in ugli2 + pheno merged = ")
#print(nrow(lifelines_ugli2_merged))

#####################
## Standardise GRS ##
#####################
lifelines_ugli2_merged$ahluwalia_cis_aff_sd <- scale(lifelines_ugli2_merged$ahluwalia_cis_ugli2)
lifelines_ugli2_merged$ahluwalia_gw_aff_sd <- scale(lifelines_ugli2_merged$ahluwalia_gw_ugli2)
lifelines_ugli2_merged$borges_proxyadd_aff_sd <- scale(lifelines_ugli2_merged$borges_proxyadd_ugli2)
lifelines_ugli2_merged$said_cis_proxyadd_aff_sd <- scale(lifelines_ugli2_merged$said_cis_proxyadd_ugli2)
lifelines_ugli2_merged$said_gw_proxyadd_aff_sd <- scale(lifelines_ugli2_merged$said_gw_proxyadd_ugli2)
lifelines_ugli2_merged$rosa_aff_sd <- scale(lifelines_ugli2_merged$rosa_ugli2)
lifelines_ugli2_merged$sarwar_aff_sd <- scale(lifelines_ugli2_merged$sarwar_ugli2)

###########################################################
## Remove Exclusions - non-European and genetic outliers ##
###########################################################
ugli2_noneuro <- fread("nonEuropeans.flagged.samples", header=FALSE)
ugli2_noneuro <- rename(ugli2_noneuro, FID=V1, IID=V2) 
ugli2_geneticoutliers <- fread("UGLI2_genetic_outliers.flagged.samples", header=TRUE)
#cat("ugli2 non euro list")
#sum(! is.na(ugli2_noneuro$IID))
#cat("ugli2 genetic outliers")
#sum(! is.na(ugli2_geneticoutliers$IID))

#Remove non-European
cat("ugli2 n prior removals")
sum(! is.na(lifelines_ugli2_merged$IID))
lifelines_ugli2_merged <- lifelines_ugli2_merged[! lifelines_ugli2_merged$IID %in% ugli2_noneuro$IID, ]
cat("ugli2 n once non-euro removed")
sum(! is.na(lifelines_ugli2_merged$IID))

#Remove genetic outliers (run included then removed).
lifelines_ugli2_merged <- lifelines_ugli2_merged[! lifelines_ugli2_merged$IID %in% ugli2_geneticoutliers$IID, ]
cat("ugli2 n once genetic outliers removed")
sum(! is.na(lifelines_ugli2_merged$IID))

################################
## Prepare for GRAMMAR Method ##
################################
#lifelines_ugli2_merged$FID <- 1
#Get FID as needed for running reml
ugli2_grm_ids <- fread("affy_snps_for_grm_final.fam", header=FALSE)
ugli2_grm_ids <- rename(ugli2_grm_ids, FID=V1, IID=V2)
lifelines_ugli2_merged <- merge(lifelines_ugli2_merged, ugli2_grm_ids, by="IID")

### Step 1: Extract files for REML (gcta to extract residuals)
#pheno file
#ugli2_merged_phenotypes_gcta <- lifelines_ugli2_merged[,c("FID","IID","VA_RES_sd","PS_RES_sd","EM_ACC_sd","WM_ACC_sd","WM_RES_sd", "RFFT_final_sd","panas_neg_total_final_sd","panas_pos_total_final_sd","minidep_1av1","minianx_1av1","minimdd_1av1","minigad_1av1","hsCRP_log","BMI_1av1","alcohol_freq_1a_q_2","Currentsmok_1a_q_2","education")]
#print(head(ugli2_merged_phenotypes_gcta),50)
#write.table(ugli2_merged_phenotypes_gcta, "ugli2_merged_pheno_for_gcta.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

#covariate pc file
#ugli2_merged_pcs_gcta <- lifelines_ugli2_merged[,c("FID","IID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")]
#write.table(ugli2_merged_pcs_gcta, "ugli2_merged_pcs_for_gcta.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

### Step 2: Load in residuals to be used as outcomes in adjusted analysis and merge with data
ugli2_resids <- fread("ugli2_resids_merged.txt", header=TRUE)
lifelines_ugli2_merged <- merge(lifelines_ugli2_merged, ugli2_resids, by="IID") 
#cat("lifelines ugli2 and resids merged")
#print(head(lifelines_ugli2_merged),20)
#cat("lifelines ugli2 + resids dims")
#print(dim(lifelines_ugli2_merged))

################################################
## Create file for REGENIE - another analysis ##
################################################
ugli2_regenie_pheno_covar <- lifelines_ugli2_merged[,c("FID", "IID", "VA_RES_sd","PS_RES_sd", "EM_ACC_sd", "WM_ACC_sd", "RFFT_final_sd", "panas_neg_total_final_sd", "panas_pos_total_final_sd", "AGE_1a_q_1", "GENDER_1a_v_1", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")]
#write.table(ugli2_regenie_pheno_covar, "ugli2_regenie_pheno_covar.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
#summary(ugli2_regenie_pheno_covar)
#print(head(ugli2_regenie_pheno_covar))

###################################################################################################
################## PART 5: Get lists of exclusions for analysis on unrelated ppl ##################
###################################################################################################
## List of IIDs for each chip - needed for KING to remove relatives for secondary analysis (only included those who can be included to maximise sample size). Note: to get list must have FID first, then IID. For CytoSNP, ran fine without FID included.
#cyto_IID <- lifelines_cyto_merged[,c("IID")]
#gsa_IID <- lifelines_gsa_merged[,c("FID", "IID")]
#ugli2_IID <- lifelines_ugli2_merged[,c("FID", "IID")]
#write.table(cyto_IID, file="cyto_IID.txt", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
#write.table(gsa_IID, file="gsa_IID_FID.txt", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
#write.table(ugli2_IID, file="ugli2_FID_IID.txt", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

## When extracting lists for exclusion based on relatedness (KING kinship matrix), only including people who have phenotype data (and thus could be included in analysis) maximises sample size, so this approach is taken. Comparison below.
## Load list of IIDs to be excluded based on KING (includes all genotyped people - removes more ppl)
#cyto_kingexclusions_177 <- fread("KING_CYTO_CUTOFF_177.king.cutoff.out.id", header=TRUE)
#cyto_kingexclusions_088 <- fread("KING_CYTO_CUTOFF_088.king.cutoff.out.id", header=TRUE)
#cyto_kingexclusions_044 <- fread("KING_CYTO_CUTOFF_044.king.cutoff.out.id", header=TRUE)
## List of IIDs to be excluded (only includes IIDs who could be included - i.e., dups/1st deg between removed)
#cyto_kingexclusions_177_phenoavail <- fread("KING_CYTO_CUTOFF_177_phenoavail.king.cutoff.out.id", header=TRUE)
#cyto_kingexclusions_088_phenoavail <- fread("KING_CYTO_CUTOFF_088_phenoavail.king.cutoff.out.id", header=TRUE)
#cyto_kingexclusions_044_phenoavail <- fread("KING_CYTO_CUTOFF_044_phenoavail.king.cutoff.out.id", header=TRUE)

#print(head(cyto_kingexclusions_177))

## Get sample size once removed related (using all geno vs IIDs who have pheno) - using those who have pheno = > sample size.
#cat("N cyto - remove related")
#lifelines_cyto_king_177 <- lifelines_cyto_merged[! lifelines_cyto_merged$IID %in% cyto_kingexclusions_177$IID, ]
#lifelines_cyto_king_088 <- lifelines_cyto_merged[! lifelines_cyto_merged$IID %in% cyto_kingexclusions_088$IID, ]
#lifelines_cyto_king_044 <- lifelines_cyto_merged[! lifelines_cyto_merged$IID %in% cyto_kingexclusions_044$IID, ]
#sum(! is.na(lifelines_cyto_king_177$IID))
#sum(! is.na(lifelines_cyto_king_088$IID))
#sum(! is.na(lifelines_cyto_king_044$IID))
#cat("N cyto - remove related when king done on IIDs who have pheno data available")
#lifelines_cyto_king_177_max <- lifelines_cyto_merged[! lifelines_cyto_merged$IID %in% cyto_kingexclusions_177_phenoavail$IID, ]
#lifelines_cyto_king_088_max <- lifelines_cyto_merged[! lifelines_cyto_merged$IID %in% cyto_kingexclusions_088_phenoavail$IID, ]
#lifelines_cyto_king_044_max <- lifelines_cyto_merged[! lifelines_cyto_merged$IID %in% cyto_kingexclusions_044_phenoavail$IID, ]
#sum(! is.na(lifelines_cyto_king_177_max$IID))
#sum(! is.na(lifelines_cyto_king_088_max$IID))
#sum(! is.na(lifelines_cyto_king_044_max$IID))

## Load list of people to be removed from analysis
cyto_kingexclusions_177_phenoavail <- fread("KING_CYTO_CUTOFF_177_phenoavail.king.cutoff.out.id", header=TRUE)
cyto_kingexclusions_088_phenoavail <- fread("KING_CYTO_CUTOFF_088_phenoavail.king.cutoff.out.id", header=TRUE)
cyto_kingexclusions_044_phenoavail <- fread("KING_CYTO_CUTOFF_044_phenoavail.king.cutoff.out.id", header=TRUE)
gsa_kingexclusions_177_phenoavail <- fread("KING_GSA_CUTOFF_177_phenoavail.king.cutoff.out.id", header=TRUE)
gsa_kingexclusions_088_phenoavail <- fread("KING_GSA_CUTOFF_088_phenoavail.king.cutoff.out.id", header=TRUE)
gsa_kingexclusions_044_phenoavail <- fread("KING_GSA_CUTOFF_044_phenoavail.king.cutoff.out.id", header=TRUE)
ugli2_kingexclusions_177_phenoavail <- fread("KING_UGLI2_CUTOFF_177_phenoavail.king.cutoff.out.id", header=TRUE)
ugli2_kingexclusions_088_phenoavail <- fread("KING_UGLI2_CUTOFF_088_phenoavail.king.cutoff.out.id", header=TRUE)
ugli2_kingexclusions_044_phenoavail <- fread("KING_UGLI2_CUTOFF_044_phenoavail.king.cutoff.out.id", header=TRUE)

#############################################################################
#################### PART 6: Analysis Prep - Merge data  ####################
#############################################################################
###############################
## Combine data across chips ##
###############################
#rename genetic risk scores to enable combining data
lifelines_cyto_merged <- rename(lifelines_cyto_merged, 
	said_cis_proxyadd_sd=said_cis_proxyadd_cyt_sd, 
	borges_proxyadd_sd=borges_proxyadd_cyt_sd,
	ahluwalia_cis_sd=ahluwalia_cis_cyt_sd,
	rosa_sd=rosa_cyt_sd,
	sarwar_sd=sarwar_cyt_sd,
	said_gw_proxyadd_sd=said_gw_proxyadd_cyt_sd,
	ahluwalia_gw_sd=ahluwalia_gw_cyt_sd)

lifelines_gsa_merged <- rename(lifelines_gsa_merged,
	said_cis_proxyadd_sd=said_cis_proxyadd_gsa_sd,
	borges_proxyadd_sd=borges_proxyadd_gsa_sd,
	ahluwalia_cis_sd=ahluwalia_cis_gsa_sd,
	rosa_sd=rosa_gsa_sd,
	sarwar_sd=sarwar_gsa_sd,
	said_gw_proxyadd_sd=said_gw_proxyadd_gsa_sd,
	ahluwalia_gw_sd=ahluwalia_gw_gsa_sd)
	
lifelines_ugli2_merged <- rename(lifelines_ugli2_merged,
	said_cis_proxyadd_sd=said_cis_proxyadd_aff_sd,
	borges_proxyadd_sd=borges_proxyadd_aff_sd,
	ahluwalia_cis_sd=ahluwalia_cis_aff_sd,
	rosa_sd=rosa_aff_sd,
	sarwar_sd=sarwar_aff_sd,
	said_gw_proxyadd_sd=said_gw_proxyadd_aff_sd,
	ahluwalia_gw_sd=ahluwalia_gw_aff_sd)

lifelines_cyto_merged_clean <- lifelines_cyto_merged[,c("IID","VA_RES_sd","PS_RES_sd","EM_ACC_sd","WM_ACC_sd","WM_RES_sd", "RFFT_final_sd","panas_neg_total_final_sd","panas_pos_total_final_sd","minidep_1av1","minianx_1av1","minimdd_1av1","minigad_1av1","hsCRP_log","BMI_1av1","alcohol_freq_1a_q_2","Currentsmok_1a_q_2","education","said_cis_proxyadd_sd","borges_proxyadd_sd","ahluwalia_cis_sd","rosa_sd","sarwar_sd","said_gw_proxyadd_sd","ahluwalia_gw_sd","VA_RES_sd_res","PS_RES_sd_res","EM_ACC_sd_res","WM_ACC_sd_res","WM_RES_sd_res", "RFFT_final_sd_res","panas_neg_total_final_sd_res","panas_pos_total_final_sd_res","BMI_1av1_res","alcohol_freq_1a_q_2_res","hsCRP_log_res","Currentsmok_1a_q_2_res","education_res","chip","AGE_1a_q_1","GENDER_1a_v_1","minidep_1av1_res","minianx_1av1_res","minimdd_1av1_res","minigad_1av1_res")]

lifelines_gsa_merged_clean <- lifelines_gsa_merged[,c("IID","VA_RES_sd","PS_RES_sd","EM_ACC_sd","WM_ACC_sd","WM_RES_sd", "RFFT_final_sd","panas_neg_total_final_sd","panas_pos_total_final_sd","minidep_1av1","minianx_1av1","minimdd_1av1","minigad_1av1","hsCRP_log","BMI_1av1","alcohol_freq_1a_q_2","Currentsmok_1a_q_2","education","said_cis_proxyadd_sd","borges_proxyadd_sd","ahluwalia_cis_sd","rosa_sd","sarwar_sd","said_gw_proxyadd_sd","ahluwalia_gw_sd","VA_RES_sd_res","PS_RES_sd_res","EM_ACC_sd_res","WM_ACC_sd_res","WM_RES_sd_res", "RFFT_final_sd_res","panas_neg_total_final_sd_res","panas_pos_total_final_sd_res","BMI_1av1_res","alcohol_freq_1a_q_2_res","hsCRP_log_res","Currentsmok_1a_q_2_res","education_res","chip","AGE_1a_q_1","GENDER_1a_v_1","minidep_1av1_res","minianx_1av1_res","minimdd_1av1_res","minigad_1av1_res")]

lifelines_ugli2_merged_clean <- lifelines_ugli2_merged[,c("IID","VA_RES_sd","PS_RES_sd","EM_ACC_sd","WM_ACC_sd","WM_RES_sd", "RFFT_final_sd","panas_neg_total_final_sd","panas_pos_total_final_sd","minidep_1av1","minianx_1av1","minimdd_1av1","minigad_1av1","hsCRP_log","BMI_1av1","alcohol_freq_1a_q_2","Currentsmok_1a_q_2","education","said_cis_proxyadd_sd","borges_proxyadd_sd","ahluwalia_cis_sd","rosa_sd","sarwar_sd","said_gw_proxyadd_sd","ahluwalia_gw_sd","VA_RES_sd_res","PS_RES_sd_res","EM_ACC_sd_res","WM_ACC_sd_res","WM_RES_sd_res", "RFFT_final_sd_res","panas_neg_total_final_sd_res","panas_pos_total_final_sd_res","BMI_1av1_res","alcohol_freq_1a_q_2_res","hsCRP_log_res","Currentsmok_1a_q_2_res","education_res","chip","AGE_1a_q_1","GENDER_1a_v_1","minidep_1av1_res","minianx_1av1_res","minimdd_1av1_res","minigad_1av1_res")]

lifelines_combined <- rbind(lifelines_cyto_merged_clean,lifelines_gsa_merged_clean,lifelines_ugli2_merged_clean)
#print(head(lifelines_combined))

#Ensure gender and chip are factors
lifelines_combined$chip <- as.factor(lifelines_combined$chip)
lifelines_combined$GENDER_1a_v_1 <- as.factor(lifelines_combined$GENDER_1a_v_1)
lifelines_combined$education <- as.factor(lifelines_combined$education)
#summary(lifelines_combined)

#Check sample size as expect
#cat("number panas_neg_total_final_sd_CYTO")
#sum(! is.na(lifelines_cyto_merged_clean$panas_neg_total_final_sd_res))
#cat("number panas_neg_total_final_sd_GSA")
#sum(! is.na(lifelines_gsa_merged_clean$panas_neg_total_final_sd_res))
#cat("number panas_neg_total_final_sd_UGLI2")
#sum(! is.na(lifelines_ugli2_merged_clean$panas_neg_total_final_sd_res))
#cat("number panas_neg_total_final_sd_combined")
#sum(! is.na(lifelines_combined$panas_neg_total_final_sd_res))
cat("N for combined data")
sum(! is.na(lifelines_combined$IID))

#Load merged PCs for non-grammar analyses
PCs_on_merged_data <- fread("PCs_merged_cytosnp+gsa+ugli2.txt", header=TRUE)
PCs_on_merged_data <- PCs_on_merged_data[,c('IID','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')]
lifelines_combined <- merge(lifelines_combined, PCs_on_merged_data, by="IID")
#cat("N for combined once PCs added")
#sum(! is.na(lifelines_combined$IID))
#print(head(lifelines_combined), 50)

##################################################################################################################
#################### Create combined dataset with related IIDs removed for sensitivity analysis ##################
##################################################################################################################
##First-degree relatives removed
combined_king_1stdeg <- fread("combined_king_1stdegree.txt",header=TRUE)
lifelines_combined_1stremoved <- lifelines_combined[! lifelines_combined$IID %in% combined_king_1stdeg$IID,]
cat("dim 1st deg king")
dim(combined_king_1stdeg)
#lifelines_combined_1stremoved <- lifelines_combined[! lifelines_combined$IID %in% cyto_kingexclusions_177_phenoavail$IID,]
#lifelines_combined_1stremoved <- lifelines_combined_1stremoved[! lifelines_combined_1stremoved$IID %in% gsa_kingexclusions_177_phenoavail$IID,]
#lifelines_combined_1stremoved <- lifelines_combined_1stremoved[! lifelines_combined_1stremoved$IID %in% ugli2_kingexclusions_177_phenoavail$IID,]
cat("sample size once 1st degree relatives removed")
sum(! is.na(lifelines_combined_1stremoved$IID))

##Second-degree relatives removed
combined_king_2nddeg <- fread("combined_king_2nddegree.txt",header=TRUE)
lifelines_combined_2ndremoved <- lifelines_combined[! lifelines_combined$IID %in% combined_king_2nddeg$IID,]
cat("dim 2nd deg king")
dim(combined_king_2nddeg)
#lifelines_combined_2ndremoved <- lifelines_combined[! lifelines_combined$IID %in% cyto_kingexclusions_088_phenoavail$IID,]
#lifelines_combined_2ndremoved <- lifelines_combined_2ndremoved[! lifelines_combined_2ndremoved$IID %in% gsa_kingexclusions_088_phenoavail$IID,]
#lifelines_combined_2ndremoved <- lifelines_combined_2ndremoved[! lifelines_combined_2ndremoved$IID %in% ugli2_kingexclusions_088_phenoavail$IID,]
cat("sample size once 2nd degree relatives removed")
sum(! is.na(lifelines_combined_2ndremoved$IID))

##Third-degree relatives removed
combined_king_3rddeg <- fread("combined_king_3rddegree.txt",header=TRUE)
lifelines_combined_3rdremoved <- lifelines_combined[! lifelines_combined$IID %in% combined_king_3rddeg$IID,]
cat("dim 3rd deg king")
dim(combined_king_3rddeg)
#lifelines_combined_3rdremoved <- lifelines_combined[! lifelines_combined$IID %in% cyto_kingexclusions_044_phenoavail$IID,]
#lifelines_combined_3rdremoved <- lifelines_combined_3rdremoved[! lifelines_combined_3rdremoved$IID %in% gsa_kingexclusions_044_phenoavail$IID,]
#lifelines_combined_3rdremoved <- lifelines_combined_3rdremoved[! lifelines_combined_3rdremoved$IID %in% gsa_kingexclusions_044_phenoavail$IID,]
cat("sample size once 3rd degree relatives removed")
sum(! is.na(lifelines_combined_3rdremoved$IID))

##############################################################
#################### PART 7: Run Analysis ####################
##############################################################
#Analyses will be run using 
#1. Adjusted (Grammar method) - residuals used as outcomes (adjusted for relatedness up to 0.125 within chip)
#2. Adjusted (Relatives removed) - used KING to identify threshold for relatedness 
#3. Unadjusted - unadjusted for relatedness within chip

######################################
## a. Check GRS associated with CRP ##
######################################
#cat("number crp in dataset")
#sum(! is.na(lifelines_combined$hsCRP_log))

####### Adjusted (grammar) #######
##Residuals from REML used as outcome (sparse GRM included in model to predict phenotype, as pcs included in grammar method, it is not included here.)
#cat("CRP GRS association with CRP - grammar adjusted, predictors = chip, grs")
#mapGLMTables(data=lifelines_combined, x=c("said_cis_proxyadd_sd","said_gw_proxyadd_sd"), y=c("hsCRP_log_res"),z=c("chip"), model.type="lm")
#cat("CRP GRS association with CRP - grammar adjusted, predictors =  chip, age, sex, grs")
#mapGLMTables(data=lifelines_combined, x=c("said_cis_proxyadd_sd","said_gw_proxyadd_sd"), y=c("hsCRP_log_res"),z=c("AGE_1a_q_1","GENDER_1a_v_1","chip"), model.type="lm")
#cat("CRP GRS association with CRP - grammar adjusted, no covariates")
#said_crp_model_cis <- lm(hsCRP_log_res~said_cis_proxyadd_sd, data=lifelines_combined)
#summary(said_crp_model_cis)
#said_crp_model_gw <- lm(hsCRP_log_res~said_gw_proxyadd_sd, data=lifelines_combined)
#summary(said_crp_model_gw)
#write.table(GRS_crp_grammar, file="LL_GRS_crp_grammar_adjustedchiponly.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
#write.table(GRS_crp_grammar_adjusted, file="LL_GRS_crp_grammar_covaradjusted.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

####### Adjusted (exclusions) #######
#Excluded one of each pair of individuals who are 1st-degree rel, up to 2nd-deg, up to 3rd
#GRS_crp_1stdeg_rem <- mapGLMTables(data=lifelines_combined_1stremoved, x=c("said_cis_proxyadd_sd","said_gw_proxyadd_sd"), y=c("hsCRP_log"),z=c("chip"), model.type="lm")
#GRS_crp_2nddeg_rem <- mapGLMTables(data=lifelines_combined_2ndremoved, x=c("said_cis_proxyadd_sd","said_gw_proxyadd_sd"), y=c("hsCRP_log"),z=c("chip"), model.type="lm")
#GRS_crp_3rddeg_rem <- mapGLMTables(data=lifelines_combined_3rdremoved, x=c("said_cis_proxyadd_sd","said_gw_proxyadd_sd"), y=c("hsCRP_log"), z=c("chip"), model.type="lm")
#write.table(GRS_crp_1stdeg_rem, file="LL_GRS_crp_1stdeg_rem_final_adjustchiponly.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
#write.table(GRS_crp_2nddeg_rem, file="LL_GRS_crp_2nddeg_rem_final_adjustchiponly.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
#write.table(GRS_crp_3rddeg_rem, file="LL_GRS_crp_3rddeg_rem_final_adjustchiponly.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
#GRS_crp_1stdeg_rem_adjusted <- mapGLMTables(data=lifelines_combined_1stremoved, x=c("said_cis_proxyadd_sd","said_gw_proxyadd_sd"), y=c("hsCRP_log"),z=c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","chip","AGE_1a_q_1","GENDER_1a_v_1"), model.type="lm")
#GRS_crp_2nddeg_rem_adjusted <- mapGLMTables(data=lifelines_combined_2ndremoved, x=c("said_cis_proxyadd_sd","said_gw_proxyadd_sd"), y=c("hsCRP_log"),z=c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","chip","AGE_1a_q_1","GENDER_1a_v_1"), model.type="lm")
#GRS_crp_3rddeg_rem_adjusted <- mapGLMTables(data=lifelines_combined_3rdremoved, x=c("said_cis_proxyadd_sd","said_gw_proxyadd_sd"), y=c("hsCRP_log"),z=c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","chip","AGE_1a_q_1","GENDER_1a_v_1"), model.type="lm")
#write.table(GRS_crp_1stdeg_rem_adjusted, file="LL_GRS_crp_1stdeg_rem_final_covaradjusted.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
#write.table(GRS_crp_2nddeg_rem_adjusted, file="LL_GRS_crp_2nddeg_rem_final_covaradjusted.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
#write.table(GRS_crp_3rddeg_rem_adjusted, file="LL_GRS_crp_3rddeg_rem_final_covaradjusted.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

########################################################
## b. Check GRS associated with potential confounders ##
########################################################
######### Adjusted (grammar) #######
#GRS_confounders_grammar <- mapGLMTables(data=lifelines_combined, x=c("said_cis_proxyadd_sd","borges_proxyadd_sd","ahluwalia_cis_sd","rosa_sd","sarwar_sd","said_gw_proxyadd_sd","ahluwalia_gw_sd"), y=c("BMI_1av1_res","Currentsmok_1a_q_2_res","education_res"),z=c("chip","AGE_1a_q_1","GENDER_1a_v_1"), model.type="lm")
#write.table(GRS_confounders_grammar, file="LL_GRS_confounders_grammar_geneticoutliersexcluded.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

####### Adjusted (exclusions) #######
##1st-degree relatives removed
#GRS_confounders_continuous_1stdeg_rem <- mapGLMTables(data=lifelines_combined_1stremoved, x=c("said_cis_proxyadd_sd","borges_proxyadd_sd","ahluwalia_cis_sd","rosa_sd","sarwar_sd","said_gw_proxyadd_sd","ahluwalia_gw_sd"), y=c("BMI_1av1"),z=c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","chip","AGE_1a_q_1","GENDER_1a_v_1"), model.type="lm")
#GRS_confounders_ordinal_1stdeg_rem <- mapGLMTables(data=lifelines_combined_1stremoved, x=c("said_cis_proxyadd_sd","borges_proxyadd_sd","ahluwalia_cis_sd","rosa_sd","sarwar_sd","said_gw_proxyadd_sd","ahluwalia_gw_sd"), y=c("education"),z=c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","chip","AGE_1a_q_1","GENDER_1a_v_1"), model.type="polr")
#GRS_confounders_binary_1stdeg_rem <- mapGLMTables(data=lifelines_combined_1stremoved, x=c("said_cis_proxyadd_sd","borges_proxyadd_sd","ahluwalia_cis_sd","rosa_sd","sarwar_sd","said_gw_proxyadd_sd","ahluwalia_gw_sd"), y=c("Currentsmok_1a_q_2"),z=c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","chip","AGE_1a_q_1","GENDER_1a_v_1"), model.type="glm")
#write.table(GRS_confounders_continuous_1stdeg_rem, file="LL_GRS_confounders_continuous_1stdeg_rem.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
#write.table(GRS_confounders_binary_1stdeg_rem, file="LL_GRS_confounders_binary_1stdeg_rem.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
#write.table(GRS_confounders_ordinal_1stdeg_rem, file="LL_GRS_confounders_ordinal_1stdeg_rem.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

##2nd-degree relatives removed
#GRS_confounders_continuous_2nddeg_rem <- mapGLMTables(data=lifelines_combined_2ndremoved, x=c("said_cis_proxyadd_sd","borges_proxyadd_sd","ahluwalia_cis_sd","rosa_sd","sarwar_sd","said_gw_proxyadd_sd","ahluwalia_gw_sd"), y=c("BMI_1av1"),z=c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","chip","AGE_1a_q_1","GENDER_1a_v_1"), model.type="lm")
#GRS_confounders_ordinal_2nddeg_rem <- mapGLMTables(data=lifelines_combined_2ndremoved, x=c("said_cis_proxyadd_sd","borges_proxyadd_sd","ahluwalia_cis_sd","rosa_sd","sarwar_sd","said_gw_proxyadd_sd","ahluwalia_gw_sd"), y=c("education"),z=c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","chip","AGE_1a_q_1","GENDER_1a_v_1"), model.type="polr")
#GRS_confounders_binary_2nddeg_rem <- mapGLMTables(data=lifelines_combined_2ndremoved, x=c("said_cis_proxyadd_sd","borges_proxyadd_sd","ahluwalia_cis_sd","rosa_sd","sarwar_sd","said_gw_proxyadd_sd","ahluwalia_gw_sd"), y=c("Currentsmok_1a_q_2"),z=c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","chip","AGE_1a_q_1","GENDER_1a_v_1"), model.type="glm")
#write.table(GRS_confounders_continuous_2nddeg_rem, file="LL_GRS_confounders_continuous_2nddeg_rem.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
#write.table(GRS_confounders_binary_2nddeg_rem, file="LL_GRS_confounders_binary_2nddeg_rem.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
#write.table(GRS_confounders_ordinal_2nddeg_rem, file="LL_GRS_confounders_ordinal_2nddeg_rem.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

##3rd-degree relatives removed
#GRS_confounders_continuous_3rddeg_rem <- mapGLMTables(data=lifelines_combined_3rdremoved, x=c("said_cis_proxyadd_sd","borges_proxyadd_sd","ahluwalia_cis_sd","rosa_sd","sarwar_sd","said_gw_proxyadd_sd","ahluwalia_gw_sd"), y=c("BMI_1av1"),z=c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","chip","AGE_1a_q_1","GENDER_1a_v_1"), model.type="lm")
#GRS_confounders_ordinal_3rddeg_rem <- mapGLMTables(data=lifelines_combined_3rdremoved, x=c("said_cis_proxyadd_sd","borges_proxyadd_sd","ahluwalia_cis_sd","rosa_sd","sarwar_sd","said_gw_proxyadd_sd","ahluwalia_gw_sd"), y=c("education"),z=c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","chip","AGE_1a_q_1","GENDER_1a_v_1"), model.type="polr")
#GRS_confounders_binary_3rddeg_rem <- mapGLMTables(data=lifelines_combined_3rdremoved, x=c("said_cis_proxyadd_sd","borges_proxyadd_sd","ahluwalia_cis_sd","rosa_sd","sarwar_sd","said_gw_proxyadd_sd","ahluwalia_gw_sd"), y=c("Currentsmok_1a_q_2"),z=c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","chip","AGE_1a_q_1","GENDER_1a_v_1"), model.type="glm")
#write.table(GRS_confounders_continuous_3rddeg_rem, file="LL_GRS_confounders_continuous_3rddeg_rem.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
#write.table(GRS_confounders_binary_3rddeg_rem, file="LL_GRS_confounders_binary_3rddeg_rem.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
#write.table(GRS_confounders_ordinal_3rddeg_rem, file="LL_GRS_confounders_ordinal_3rddeg_rem.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

###########################################
## c. Check GRS associated with outcomes ##
###########################################
# Note: adjusted analysis (sex and age) slightly lower N (n=18 without gender data)
####################################
######### Adjusted (grammar) #######
####################################
#GRS_outcomes_continuous_grammar <- mapGLMTables(data=lifelines_combined, y=c("VA_RES_sd_res","PS_RES_sd_res","EM_ACC_sd_res","WM_ACC_sd_res","RFFT_final_sd_res","panas_neg_total_final_sd_res","panas_pos_total_final_sd_res","minidep_1av1_res","minianx_1av1_res","minimdd_1av1_res","minigad_1av1_res"), x=c("said_cis_proxyadd_sd","borges_proxyadd_sd","ahluwalia_cis_sd","rosa_sd","sarwar_sd","said_gw_proxyadd_sd","ahluwalia_gw_sd"),z=c("AGE_1a_q_1","GENDER_1a_v_1","chip"))
#write.table(GRS_outcomes_continuous_grammar, file="LL_GRS_outcomes_continuous_grammar_geneticoutliersexcluded.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

##Sanity checks - why lower when include age and gender - n=18 without gender data
#cat("n without age data")
#sum(is.na(lifelines_combined$AGE_1a_q_1))
#cat("n without gender data")
#sum(is.na(lifelines_combined$GENDER_1a_v_1))

#####################################
######### Adjusted (exclusions) #####
#####################################
##1st-degree relatives (based on KING thresholds) removed
#GRS_outcomes_continuous_1stdeg_rem <- mapGLMTables(data=lifelines_combined_1stremoved, y=c("VA_RES_sd","PS_RES_sd","EM_ACC_sd","WM_ACC_sd","RFFT_final_sd","panas_neg_total_final_sd","panas_pos_total_final_sd"), x=c("said_cis_proxyadd_sd","borges_proxyadd_sd","ahluwalia_cis_sd","rosa_sd","sarwar_sd","said_gw_proxyadd_sd","ahluwalia_gw_sd"),z=c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","AGE_1a_q_1","GENDER_1a_v_1","chip"))
#GRS_outcomes_binary_1stdeg_rem <- mapGLMTables(data=lifelines_combined_1stremoved, y=c("minidep_1av1","minianx_1av1","minimdd_1av1","minigad_1av1"), x=c("said_cis_proxyadd_sd","borges_proxyadd_sd","ahluwalia_cis_sd","rosa_sd","sarwar_sd","said_gw_proxyadd_sd","ahluwalia_gw_sd"),z=c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","AGE_1a_q_1","GENDER_1a_v_1","chip"), model.type="glm")
#write.table(GRS_outcomes_continuous_1stdeg_rem, file="sens_correct_unrel/LL_GRS_outcomes_continuous_1stdeg_removed_31072025.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
#write.table(GRS_outcomes_binary_1stdeg_rem, file="sens_correct_unrel/LL_GRS_outcomes_binary_1stdeg_removed_31072025.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

##2nd-degree relatives and above removed
#GRS_outcomes_continuous_2nddeg_rem <- mapGLMTables(data=lifelines_combined_2ndremoved, y=c("VA_RES_sd","PS_RES_sd","EM_ACC_sd","WM_ACC_sd","RFFT_final_sd","panas_neg_total_final_sd","panas_pos_total_final_sd"), x=c("said_cis_proxyadd_sd","borges_proxyadd_sd","ahluwalia_cis_sd","rosa_sd","sarwar_sd","said_gw_proxyadd_sd","ahluwalia_gw_sd"),z=c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","AGE_1a_q_1","GENDER_1a_v_1","chip"))
#GRS_outcomes_binary_2nddeg_rem <- mapGLMTables(data=lifelines_combined_2ndremoved, y=c("minidep_1av1","minianx_1av1","minimdd_1av1","minigad_1av1"), x=c("said_cis_proxyadd_sd","borges_proxyadd_sd","ahluwalia_cis_sd","rosa_sd","sarwar_sd","said_gw_proxyadd_sd","ahluwalia_gw_sd"),z=c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","AGE_1a_q_1","GENDER_1a_v_1","chip"), model.type="glm")
#write.table(GRS_outcomes_continuous_2nddeg_rem, file="sens_correct_unrel/LL_GRS_outcomes_continuous_2nddeg_removed_31072025.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
#write.table(GRS_outcomes_binary_2nddeg_rem, file="sens_correct_unrel/LL_GRS_outcomes_binary_2nddeg_removed_31072025.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

##3rd-degree relatives and above removed
#GRS_outcomes_continuous_3rddeg_rem <- mapGLMTables(data=lifelines_combined_3rdremoved, y=c("VA_RES_sd","PS_RES_sd","EM_ACC_sd","WM_ACC_sd","RFFT_final_sd","panas_neg_total_final_sd","panas_pos_total_final_sd"), x=c("said_cis_proxyadd_sd","borges_proxyadd_sd","ahluwalia_cis_sd","rosa_sd","sarwar_sd","said_gw_proxyadd_sd","ahluwalia_gw_sd"),z=c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","AGE_1a_q_1","GENDER_1a_v_1","chip"))
#GRS_outcomes_binary_3rddeg_rem <- mapGLMTables(data=lifelines_combined_3rdremoved, y=c("minidep_1av1","minianx_1av1","minimdd_1av1","minigad_1av1"), x=c("said_cis_proxyadd_sd","borges_proxyadd_sd","ahluwalia_cis_sd","rosa_sd","sarwar_sd","said_gw_proxyadd_sd","ahluwalia_gw_sd"),z=c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","AGE_1a_q_1","GENDER_1a_v_1","chip"), model.type="glm")
#write.table(GRS_outcomes_continuous_3rddeg_rem, file="sens_correct_unrel/LL_GRS_outcomes_continuous_3rddeg_removed_31072025.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
#write.table(GRS_outcomes_binary_3rddeg_rem, file="sens_correct_unrel/LL_GRS_outcomes_binary_3rddeg_removed_31072025.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

##########################
####### Unadjusted #######
##########################
#GRS_outcomes_continuous_unadjusted <- mapGLMTables(data=lifelines_combined, y=c("VA_RES_sd","PS_RES_sd","EM_ACC_sd","WM_ACC_sd", "RFFT_final_sd","panas_neg_total_final_sd","panas_pos_total_final_sd"), x=c("said_cis_proxyadd_sd","borges_proxyadd_sd","ahluwalia_cis_sd","rosa_sd","sarwar_sd","said_gw_proxyadd_sd","ahluwalia_gw_sd"),z=c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","AGE_1a_q_1","GENDER_1a_v_1","chip"))
#GRS_outcomes_binary_unadjusted <- mapGLMTables(data=lifelines_combined, y=c("minidep_1av1","minianx_1av1","minimdd_1av1","minigad_1av1"), x=c("said_cis_proxyadd_sd","borges_proxyadd_sd","ahluwalia_cis_sd","rosa_sd","sarwar_sd","said_gw_proxyadd_sd","ahluwalia_gw_sd"),z=c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","AGE_1a_q_1","GENDER_1a_v_1","chip"), model.type="glm")
#write.table(GRS_outcomes_continuous_unadjusted, file="LL_GRS_outcomes_continuous_unadjusted_geneticoutliersexcluded.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
#write.table(GRS_outcomes_binary_unadjusted, file="LL_GRS_outcomes_binary_unadjusted_geneticoutliersexcluded.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

######################
### d. 2sls for crp ##
######################
####### adjusted (grammar) #######
#tsls_panasneg_saidcis_grammar_agesex <- summary(ivreg(panas_neg_total_final_sd_res ~ hsCRP_log_res + chip + AGE_1a_q_1 + GENDER_1a_v_1 | said_cis_proxyadd_sd + chip + AGE_1a_q_1 + GENDER_1a_v_1, data=lifelines_combined), diagnostics=T)
#capture.output(tsls_panasneg_saidcis_grammar_agesex, file="LL_tsls_panasneg_saidcis_grammar_agesex.txt")
#tsls_minianx_saidcis_grammar_agesex <- summary(ivreg(minianx_1av1_res ~ hsCRP_log_res + chip + AGE_1a_q_1 + GENDER_1a_v_1 | said_cis_proxyadd_sd + chip + AGE_1a_q_1 + GENDER_1a_v_1, data=lifelines_combined), diagnostics=T)
#capture.output(tsls_minianx_saidcis_grammar_agesex, file="LL_tsls_minianx_saidcis_grammar_agesex.txt")

#added PCs - have been adjusted for in outcome and CRP so should not be necessary. But check here.
#tsls_panasneg_saidcis_grammar_PCadded <- summary(ivreg(panas_neg_total_final_sd_res ~ hsCRP_log_res + chip + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + AGE_1a_q_1 + GENDER_1a_v_1 | said_cis_proxyadd_sd + chip + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + AGE_1a_q_1 + GENDER_1a_v_1, data=lifelines_combined), diagnostics=T)
#capture.output(tsls_panasneg_saidcis_grammar_PCadded, file="LL_tsls_panasneg_saidcis_grammar_PCadded.txt")
#tsls_minianx_saidcis_grammar_PCadded <- summary(ivreg(minianx_1av1_res ~ hsCRP_log_res + chip + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + AGE_1a_q_1 + GENDER_1a_v_1 | said_cis_proxyadd_sd + chip + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +  AGE_1a_q_1 + GENDER_1a_v_1, data=lifelines_combined), diagnostics=T)
#capture.output(tsls_minianx_saidcis_grammar_PCadded, file="LL_tsls_minianx_saidcis_grammar_PCadded.txt")

# 1st removed
#tsls_panasneg_saidcis_1stremoved <- summary(ivreg(panas_neg_total_final_sd ~ hsCRP_log + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + chip + AGE_1a_q_1 + GENDER_1a_v_1 | said_cis_proxyadd_sd + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + chip + AGE_1a_q_1 + GENDER_1a_v_1, data=lifelines_combined_1stremoved), diagnostics=T)
#capture.output(tsls_panasneg_saidcis_1stremoved, file="sens_correct_unrel/LL_tsls_panasneg_saidcis_1stremoved_31072025.txt")
#tsls_minianx_saidcis_1stremoved <- summary(ivreg(minianx_1av1 ~ hsCRP_log + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + chip + AGE_1a_q_1 + GENDER_1a_v_1 | said_cis_proxyadd_sd + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + chip + AGE_1a_q_1 + GENDER_1a_v_1, data=lifelines_combined_1stremoved), diagnostics=T)
#capture.output(tsls_minianx_saidcis_1stremoved, file="sens_correct_unrel/LL_tsls_minianx_saidcis_1stremoved_31072025.txt")

#2nd removed
#tsls_panasneg_saidcis_2ndremoved <- summary(ivreg(panas_neg_total_final_sd ~ hsCRP_log + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + chip + AGE_1a_q_1 + GENDER_1a_v_1 | said_cis_proxyadd_sd + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + chip + AGE_1a_q_1 + GENDER_1a_v_1, data=lifelines_combined_2ndremoved), diagnostics=T)
#capture.output(tsls_panasneg_saidcis_2ndremoved, file="sens_correct_unrel/LL_tsls_panasneg_saidcis_2ndremoved_31072025.txt")
#tsls_minianx_saidcis_2ndremoved <- summary(ivreg(minianx_1av1 ~ hsCRP_log + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + chip + AGE_1a_q_1 + GENDER_1a_v_1 | said_cis_proxyadd_sd + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + chip + AGE_1a_q_1 + GENDER_1a_v_1, data=lifelines_combined_2ndremoved), diagnostics=T)
#capture.output(tsls_minianx_saidcis_2ndremoved, file="sens_correct_unrel/LL_tsls_minianx_saidcis_2ndremoved_31072025.txt")

# 3rd removed
#tsls_panasneg_saidcis_3rdremoved <- summary(ivreg(panas_neg_total_final_sd ~ hsCRP_log + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + chip + AGE_1a_q_1 + GENDER_1a_v_1 | said_cis_proxyadd_sd + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + chip + AGE_1a_q_1 + GENDER_1a_v_1, data=lifelines_combined_3rdremoved), diagnostics=T)
#capture.output(tsls_panasneg_saidcis_3rdremoved, file="sens_correct_unrel/LL_tsls_panasneg_saidcis_3rdremoved_31072025.txt")
#tsls_minianx_saidcis_3rdremoved <- summary(ivreg(minianx_1av1 ~ hsCRP_log + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + chip + AGE_1a_q_1 + GENDER_1a_v_1 | said_cis_proxyadd_sd + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + chip + AGE_1a_q_1 + GENDER_1a_v_1, data=lifelines_combined_3rdremoved), diagnostics=T)
#capture.output(tsls_minianx_saidcis_3rdremoved, file="sens_correct_unrel/LL_tsls_minianx_saidcis_3rdremoved_31072025.txt")

#############################################
## Check Phenotype data value on each chip ##
#############################################
#check sample size
#cat("n cyto")
#sum(! is.na(lifelines_cyto_merged$IID))
#cat("n gsa")
#sum(! is.na(lifelines_gsa_merged$IID))
#cat("n ugli2")
#sum(! is.na(lifelines_ugli2_merged$IID))

cat("n below <18")
sum((lifelines_combined$AGE_1a_q_1 == 18) & ! is.na(lifelines_combined$RFFT_final_sd))
cat("n 18-60")
sum((lifelines_combined$AGE_1a_q_1 >18) & (lifelines_combined$AGE_1a_q_1 <= 60) & ! is.na(lifelines_combined$RFFT_final_sd))
cat("n 60+")
sum((lifelines_combined$AGE_1a_q_1 >60) & ! is.na(lifelines_combined$RFFT_final_sd))


cat("n combined")
sum(! is.na(lifelines_combined$IID))


cat("n CRP and any anxiety")
sum(! is.na(lifelines_combined$hsCRP_log_res) & ! is.na(lifelines_combined$minianx_1av1) & ! is.na(lifelines_combined$AGE_1a_q_1) & ! is.na(lifelines_combined$GENDER_1a_v_1))
cat("n CRP and panas neg")
sum(! is.na(lifelines_combined$hsCRP_log_res) & ! is.na(lifelines_combined$panas_neg_total_final_sd) & ! is.na(lifelines_combined$AGE_1a_q_1) & ! is.na(lifelines_combined$GENDER_1a_v_1))

#check duplicates = 0
#cat("number cyto duplicates")
#print(sum(duplicated(lifelines_cyto_merged$IID)))
#cat("number gsa duplicates")
#print(sum(duplicated(lifelines_gsa_merged$IID)))
#cat("number ugli2 duplicates")
#print(sum(duplicated(lifelines_ugli2_merged$IID)))

#check dimensions = as expect
cat("cyto pheno dim")
print(dim(lifelines_cyto_merged))
cat("gsa pheno dim")
print(dim(lifelines_gsa_merged))
cat("ugli2 pheno dim")
print(dim(lifelines_ugli2_merged))
cat("lifelines combined dim")
print(dim(lifelines_combined))


### Get descriptives on raw data
lifelines_cyto_merged_raw <- lifelines_cyto_merged[,c("IID","VA_RES","PS_RES","EM_ACC","WM_ACC","WM_RES", "RFFT_final","panas_neg_total_final","panas_pos_total_final","minidep_1av1","minianx_1av1","minimdd_1av1","minigad_1av1","hsCRP","BMI_1av1","alcohol_freq_1a_q_2","Currentsmok_1a_q_2","education","AGE_1a_q_1","GENDER_1a_v_1")]
lifelines_gsa_merged_raw <- lifelines_gsa_merged[,c("IID","VA_RES","PS_RES","EM_ACC","WM_ACC","WM_RES", "RFFT_final","panas_neg_total_final","panas_pos_total_final","minidep_1av1","minianx_1av1","minimdd_1av1","minigad_1av1","hsCRP","BMI_1av1","alcohol_freq_1a_q_2","Currentsmok_1a_q_2","education","AGE_1a_q_1","GENDER_1a_v_1")]
lifelines_ugli2_merged_raw <- lifelines_ugli2_merged[,c("IID","VA_RES","PS_RES","EM_ACC","WM_ACC","WM_RES", "RFFT_final","panas_neg_total_final","panas_pos_total_final","minidep_1av1","minianx_1av1","minimdd_1av1","minigad_1av1","hsCRP","BMI_1av1","alcohol_freq_1a_q_2","Currentsmok_1a_q_2","education","AGE_1a_q_1","GENDER_1a_v_1")]
lifelines_combined_raw <- rbind(lifelines_cyto_merged_raw,lifelines_gsa_merged_raw,lifelines_ugli2_merged_raw)

##AGE
cat("age")
describe(lifelines_combined_raw$AGE_1a_q_1, na.rm=TRUE)
mean(lifelines_combined_raw$AGE_1a_q_1, na.rm=TRUE)
sd(lifelines_combined_raw$AGE_1a_q_1, na.rm=TRUE)
sum(! is.na(lifelines_combined_raw$AGE_1a_q_1))

##SEX
cat("sex")
cat("sex outputs = ", table(lifelines_combined_raw$GENDER_1a_v_1))
cat("sex is 0 = ", sum(lifelines_combined_raw$GENDER_1a_v_1 == "0", na.rm=TRUE))
cat("sex is 1 =", sum(lifelines_combined_raw$GENDER_1a_v_1 == "1", na.rm=TRUE))
sum(! is.na(lifelines_combined_raw$GENDER_1a_v_1))

##BMI
cat("BMI")
mean(lifelines_combined_raw$BMI_1av1, na.rm=TRUE)
sd(lifelines_combined_raw$BMI_1av1, na.rm=TRUE)
sum(! is.na(lifelines_combined_raw$BMI_1av1))

##Education
cat("education")
#describe(lifelines_combined_raw$education, na.rm=TRUE)
cat("education outputs = ", table(lifelines_combined_raw$education))
cat("education is 0 = ", sum(lifelines_combined_raw$education == "0", na.rm=TRUE))
cat("education is 1 =", sum(lifelines_combined_raw$education == "1", na.rm=TRUE))
cat("education is 2 =", sum(lifelines_combined_raw$education == "2", na.rm=TRUE))
cat("education is 3 =", sum(lifelines_combined_raw$education == "3", na.rm=TRUE))
sum(! is.na(lifelines_combined_raw$education))

##CRP
cat("CRP data")
#describe(lifelines_combined_raw$hsCRP, na.rm=TRUE)
sum(! is.na(lifelines_combined_raw$hsCRP))
mean(lifelines_combined_raw$hsCRP, na.rm=TRUE)
sd(lifelines_combined_raw$hsCRP, na.rm=TRUE)
min(lifelines_combined_raw$hsCRP, na.rm=TRUE)
max(lifelines_combined_raw$hsCRP, na.rm=TRUE)
median(lifelines_combined_raw$hsCRP, na.rm=TRUE)
IQR(lifelines_combined_raw$hsCRP, na.rm=TRUE)

#typeof(lifelines_combined_raw$hsCRP)

##Anxiety
cat("Anxiety outputs = ", table(lifelines_combined_raw$minianx_1av1))
cat("Anxiety diagnosis = ", sum(lifelines_combined_raw$minianx_1av1 == "1", na.rm=TRUE))
cat("Anxiety non diagnosis =", sum(lifelines_combined_raw$minianx_1av1 == "0", na.rm = TRUE))
#typeof(lifelines_combined_raw$minianx_1av1)

##GAD
cat("GAD outputs = ", table(lifelines_combined_raw$minigad_1av1))
cat("GAD diagnosis = ", sum(lifelines_combined_raw$minigad_1av1 == "1", na.rm=TRUE))
cat("GAD non diagnosis =", sum(lifelines_combined_raw$minigad_1av1 == "0", na.rm = TRUE))
#typeof(lifelines_combined_raw$minigad_1av1)

##MDD
cat("MDD outputs = ", table(lifelines_combined_raw$minimdd_1av1))
cat("MDD diagnosis = ", sum(lifelines_combined_raw$minimdd_1av1 == "1", na.rm=TRUE))
cat("MDD non diagnosis =", sum(lifelines_combined_raw$minimdd_1av1 == "0", na.rm = TRUE))
#typeof(lifelines_combined_raw$minimdd_1av1)
#
##Depression
cat("Dep outputs = ", table(lifelines_combined_raw$minidep_1av1))
cat("Dep diagnosis = ", sum(lifelines_combined_raw$minidep_1av1 == "1", na.rm=TRUE))
cat("Dep non diagnosis =", sum(lifelines_combined_raw$minidep_1av1 == "0", na.rm = TRUE))
#typeof(lifelines_combined_raw$minidep_1av1)

##Cognition
##Executive function (RFFT)
cat("RFFT ")
#describe(lifelines_combined_raw$RFFT_final, na.rm=TRUE)
sum(! is.na(lifelines_combined_raw$RFFT_final))
mean(lifelines_combined_raw$RFFT_final, na.rm=TRUE)
sd(lifelines_combined_raw$RFFT_final, na.rm=TRUE)
min(lifelines_combined_raw$RFFT_final, na.rm=TRUE)
max(lifelines_combined_raw$RFFT_final, na.rm=TRUE) 
#typeof(lifelines_combined_raw$RFFT_final)

##Identification Task (01 code) - visual attention
##primary outcome = reaction time (ms) normalised using log10 transform 
cat("identification task")
#describe(lifelines_combined_raw$VA_RES, na.rm=TRUE)
sum(! is.na(lifelines_combined_raw$VA_RES))
mean(lifelines_combined_raw$VA_RES, na.rm=TRUE)
sd(lifelines_combined_raw$VA_RES, na.rm=TRUE)
min(lifelines_combined_raw$VA_RES, na.rm=TRUE)
max(lifelines_combined_raw$VA_RES, na.rm=TRUE) 
#typeof(lifelines_combined_raw$VA_RES)

##Detection Task (02 code) - processing/psychomotor speed
##primary outcome = reaction time (ms) normalised using log10 transform
cat("detection task")
#describe(lifelines_combined_raw$PS_RES, na.rm=TRUE)
sum(! is.na(lifelines_combined_raw$PS_RES))
mean(lifelines_combined_raw$PS_RES, na.rm=TRUE)
sd(lifelines_combined_raw$PS_RES, na.rm=TRUE)
min(lifelines_combined_raw$PS_RES, na.rm=TRUE)
max(lifelines_combined_raw$PS_RES, na.rm=TRUE) 
#typeof(lifelines_combined_raw$PS_RES)

##One Card Learning Task (03 code) - visual learning
##primary outcome = accuracy - arcsine sq root hit rate (correct/total)
cat("One card learning task")
describe(lifelines_combined_raw$EM_ACC, na.rm=TRUE)
sum(! is.na(lifelines_combined_raw$EM_ACC))
mean(lifelines_combined_raw$EM_ACC, na.rm=TRUE)
sd(lifelines_combined_raw$EM_ACC, na.rm=TRUE)
min(lifelines_combined_raw$EM_ACC, na.rm=TRUE)
max(lifelines_combined_raw$EM_ACC, na.rm=TRUE) 
#typeof(lifelines_combined_raw$EM_ACC)

##One Back Task (04 code) - working memory (RT)
##primary outcome = reaction time transformed
cat("one back task - RT")
describe(lifelines_combined_raw$WM_RES, na.rm=TRUE)
sum(! is.na(lifelines_combined_raw$WM_RES))
mean(lifelines_combined_raw$WM_RES, na.rm=TRUE)
sd(lifelines_combined_raw$WM_RES, na.rm=TRUE)
min(lifelines_combined_raw$WM_RES, na.rm=TRUE)
max(lifelines_combined_raw$WM_RES, na.rm=TRUE) 
#typeof(lifelines_combined_raw$WM_RES)

##One Back Task (04 code) - working memory (ACCURACY)
##primary outcome = accuracy - arcsine sq root hit rate (correct/total)
cat("one back task - accuracy transformed")
describe(lifelines_combined_raw$WM_ACC, na.rm=TRUE)
sum(! is.na(lifelines_combined_raw$WM_ACC))
mean(lifelines_combined_raw$WM_ACC, na.rm=TRUE)
sd(lifelines_combined_raw$WM_ACC, na.rm=TRUE)
min(lifelines_combined_raw$WM_ACC, na.rm=TRUE)
max(lifelines_combined_raw$WM_ACC, na.rm=TRUE) 
#typeof(lifelines_combined_raw$WM_ACC)

##PANAS - positive
cat("panas - positive")
#describe(lifelines_combined_raw$panas_pos_total_final, na.rm=TRUE)
sum(! is.na(lifelines_combined_raw$panas_pos_total_final))
mean(lifelines_combined_raw$panas_pos_total_final, na.rm=TRUE)
sd(lifelines_combined_raw$panas_pos_total_final, na.rm=TRUE)
min(lifelines_combined_raw$panas_pos_total_final, na.rm=TRUE)
max(lifelines_combined_raw$panas_pos_total_final, na.rm=TRUE)
#typeof(lifelines_combined_raw$panas_pos_total_final)

##PANAS - negative
cat("panas - negative")
#describe(lifelines_combined_raw$panas_neg_total_final, na.rm=TRUE)
sum(! is.na(lifelines_combined_raw$panas_neg_total_final))
mean(lifelines_combined_raw$panas_neg_total_final, na.rm=TRUE)
sd(lifelines_combined_raw$panas_neg_total_final, na.rm=TRUE)
min(lifelines_combined_raw$panas_neg_total_final, na.rm=TRUE)
max(lifelines_combined_raw$panas_neg_total_final, na.rm=TRUE)
#typeof(lifelines_combined_raw$panas_neg_total_final)

#cat("is gender a factor")
#is.factor(lifelines_combined$GENDER_1a_v_1)
#cat("is chip factor")
#is.factor(lifelines_combined$chip)
#cat("is age factor")
#is.factor(lifelines_combined$AGE_1a_q_1)
#
#
#cat("cyto FID") # distinct values
#describe(lifelines_cyto_merged$FID)
#cat("gsa FID") # IID and FID identical
#sum(! lifelines_gsa_merged$IID == lifelines_gsa_merged$FID)
#cat("ugli2 FID") #all same
#sum(! duplicated(lifelines_ugli2_merged$FID))==1
#
#cat("count unique gsa values")
#length(unique(lifelines_gsa_merged$FID))
#

########################################################################
## Old list of duplicate and 1st-degree rel file between cyto and gsa ##
########################################################################
#cytosnp_duplicates <- fread("LL_Cyto_with_duplicates_in_UGLI.txt", header=TRUE)
#cytosnp_1stdegree_relatives_UGLI <- fread("LL_Cyto_with_1st-degree-relatives_in_UGLI.txt")
#cytosnp_duplicates_and_1stdegree_UGLI <- fread("CytoSNP_duplicates+1stdgr_in_UGLI.txt", header=FALSE)

#cytosnp_duplicates <- rename(cytosnp_duplicates, FID=FID1, IID=IID1)
#cytosnp_1stdegree_relatives_UGLI <- rename(cytosnp_1stdegree_relatives_UGLI, IID=V2)
#cytosnp_duplicates_and_1stdegree_UGLI <- rename(cytosnp_duplicates_and_1stdegree_UGLI, IID=V1)
#cytosnp_duplicates <- cytosnp_duplicates[, c('IID')]
#cytosnp_1stdegree_relatives_UGLI <- cytosnp_1stdegree_relatives_UGLI[, c('IID')]


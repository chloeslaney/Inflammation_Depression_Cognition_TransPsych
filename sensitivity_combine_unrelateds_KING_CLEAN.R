## Load KING lists to be removed for GRS secondary analysis - combine across chips to form single lists
## Chloe Slaney

#####################
# PRELIMINARY STEPS #
#####################
library(cat)
library(data.table)
library(tidyverse)
library(readr)
library(dplyr)
library(utils)
library(httr)
library(rlist)
library(haven)
library(devtools)

############################
## Set working directory  ##
############################
setwd(" ## ADD WORKING DIR HERE ##")
rm(list = ls())

####################################################
## 1. Load KING list of individuals to be removed ##
###################################################
cyto_1st_degree <- fread("KING_CYTO_CUTOFF_177_phenoavail.king.cutoff.out.id", header=TRUE)
cyto_2nd_degree <- fread("KING_CYTO_CUTOFF_088_phenoavail.king.cutoff.out.id", header=TRUE)
cyto_3rd_degree <- fread("KING_CYTO_CUTOFF_044_phenoavail.king.cutoff.out.id", header=TRUE)

gsa_1st_degree <- fread("KING_GSA_CUTOFF_177_phenoavail.king.cutoff.out.id", header=TRUE)
gsa_2nd_degree <- fread("KING_GSA_CUTOFF_088_phenoavail.king.cutoff.out.id", header=TRUE)
gsa_3rd_degree <- fread("KING_GSA_CUTOFF_044_phenoavail.king.cutoff.out.id", header=TRUE)
 
affy_1st_degree <- fread("KING_AFFY_CUTOFF_177_phenoavail.king.cutoff.out.id", header=TRUE)
affy_2nd_degree <- fread("KING_AFFY_CUTOFF_088_phenoavail.king.cutoff.out.id", header=TRUE)
affy_3rd_degree <- fread("KING_AFFY_CUTOFF_044_phenoavail.king.cutoff.out.id", header=TRUE)

###################################
## 2. Combine lists across chips ##
###################################
combined_1st_degree <- rbind(cyto_1st_degree,gsa_1st_degree,affy_1st_degree) 
combined_2nd_degree <- rbind(cyto_2nd_degree,gsa_2nd_degree,affy_2nd_degree)
combined_3rd_degree <- rbind(cyto_3rd_degree,gsa_3rd_degree,affy_3rd_degree)

combined_1st_degree <-combined_1st_degree[,c("IID")]
combined_2nd_degree <- combined_2nd_degree[,c("IID")]
combined_3rd_degree <- combined_3rd_degree[,c("IID")]

############################
## 3. Check combined data ##
############################
cat("combined 1st degree list")
print(head(combined_1st_degree))
cat("combined 2nd degree list")
print(head(combined_2nd_degree))
cat("combined 3rd degree list")
print(head(combined_3rd_degree))

############################
## 4. Save combined lists ##
############################
#write.table(combined_1st_degree, file="combined_king_1stdegree_out.txt", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
#write.table(combined_2nd_degree, file="combined_king_2nddegree_out.txt", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
#write.table(combined_3rd_degree, file="combined_king_3rddegree_out.txt", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

#####################################################################
# ~~~~~~~~~~~~~~~~~~
# Tumor subtype and cell type independent DNA methylation alterations 
# associated with stage progression in invasive breast carcinoma 
# ~~~~~~~~~~~~~~~~~~
# Way, G., Johnson, K., Christensen, B. 2015

# The script will apply RefFreeEWAS to the given model subset
# Using the packages 'RefFreeEWAS' and 'isva' as distributed by 
# Houseman et al 2014 and Teschendorff et al 2011
#####################################################################

################################
# Read in command line arguments
################################
args <- commandArgs(trailingOnly = T)
bootstraps <- args[1]
stage <- args[2]
subtype <- args[3]

set.seed(123)
################################
# Load Libraries
################################
library(RefFreeEWAS)
# library(RnBeads)
library(readr)
library(isva)
library(ggplot2)
library(reshape)
library(gridExtra)
library(qvalue)
source("II.RefFreeEWAS/Scripts/Functions/doRefFree_functions2.R")

################################
# Load and Subset Data
################################
f.location <- "I.Data_Processing/Data/"
files <- list.files(f.location)
beta.file <- paste(f.location, files[grep(subtype, files)], sep = "")
cat(beta.file, "\n")

# Load the file of PAM50 subtype specific betas
beta2 <- read.table(beta.file, row.names = 1, header = T, sep = "\t", stringsAsFactors = F, nrows = 400000, comment.char = "")

# Load covariate file
covariates <- read.table("I.Data_Processing/Files/BRCAtarget_covariates.csv", row.names = 1, header = T, sep = ",", stringsAsFactors = F)

# The colnames for the beta file have an "X" appended to the beginning of each basename, remove it
colnames(beta2) <- substr(colnames(beta2), 2, nchar(colnames(beta2)))
rownames(covariates) <- covariates$Basename

# Subset the covariate data to only the samples in the beta file
covariates <- covariates[intersect(rownames(covariates), colnames(beta2)),]

# Subset covariate data based on stage assignment
if(stage == "low") {
  stageSub <- c("Stage I", "Stage IA", "Stage IB", "Stage II", "Stage IIA", "Stage IIB")
  } else if(stage == "high") {
    stageSub <- c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC", "Stage IV")
  }

# This function is taken from "doRefFree_functions.R"
stageCov <- subsetStage(covariates, stageSub)

# How many samples
n <- nrow(stageCov[stageCov$sample.type == "Primary Tumor", ])

################################
# Run customRefFree
################################
cat("stage: ", stage, "\n")
cat("bootstraps: ", bootstraps, "\n")
file = paste(subtype, stage, "_TCGA-BRCA", sep = "")
cat(file, "\n")

# subset betas
newBeta <- beta2[ ,rownames(stageCov)]
newBeta <- data.matrix(newBeta)

# Choose only the top 10,000 rows (i.e. those with the largest var)
DMGR_Var <- apply(newBeta, 1, var)
rankvar<-rank(-DMGR_Var)
Y_shortened<-newBeta[rankvar<=10000,]

# Step 1 - 2: Alternate fixing Mu and Omega by iterating from 2 to 10 (Kmax) cell types
DMGR_RefFree_Array <- RefFreeCellMixArray(Y_shortened, Klist=2:10, iters=25)
Y_shortened2<-Y_shortened
sapply(DMGR_RefFree_Array,deviance,Y=Y_shortened)

#Using the full betas plus the shortened most variable probes to infer the cell types
DMGR_RefFree_Array2 <- RefFreeCellMixArray(newBeta, Klist=2:10, iters=25, Yfinal=Y_shortened2)
sapply(DMGR_RefFree_Array2,deviance,Y=Y_shortened2)

celprop<-DMGR_RefFree_Array2[[2]]$Omega


# Step 3: Bootstrap method for determining the optimal number of Classes K
# May need to increase the number of bootstraps though?
# R is the number of bootstrapped vectors by default is five, boots is the number of iterations, by default is five #1500 to obtain 0.001 p-value
RefFree_DMGR_Boots = RefFreeCellMixArrayDevianceBoots(DMGR_RefFree_Array2, Y_shortened2, R=1500, bootstrapIterations=10)

RefFree_DMGR_Boots
RefFree_DMGR_Boots2<-RefFree_DMGR_Boots
apply(RefFree_DMGR_Boots[-1,],2,mean,trim=0.25)
which.min(apply(RefFree_DMGR_Boots[-1,],2,mean,trim=0.25))

class(RefFree_DMGR_Boots)

# Save Results
save(DMGR_RefFree_Array, DMGR_RefFree_Array2, RefFree_DMGR_Boots,
     file=paste("/global/scratch/atitus/data/", subtype, "_", stage, "_RefFree2.0List_05Oct2016.RData", sep = ""), 
     compress=TRUE)







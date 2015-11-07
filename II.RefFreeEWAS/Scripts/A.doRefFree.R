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
library(isva)
library(ggplot2)
library(reshape)
library(gridExtra)
library(qvalue)
source("II.RefFreeEWAS/Scripts/Functions/doRefFree_functions.R")

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
returnList <- customRefFree(covariates = stageCov, betas = beta2, bootstraps = bootstraps)

################################
# Extract Results
################################
results <- data.frame(returnList[[2]])

adj.q <- qvalue(results[ ,3])
unadj.q <- qvalue(results[ ,4])

# Prepare dataframe for export
adjusted <- cbind(adj.q$pvalues, adj.q$qvalues, results[ ,1])
unadjusted <- cbind(unadj.q$pvalues, unadj.q$qvalues, results[ ,2])

# Name the rows and columns for each model
rownames(adjusted) <- rownames(unadjusted) <- rownames(beta2)
colnames(adjusted) <- colnames(unadjusted) <- c("pvalues", "qvalues", "beta")

# print out summary of the adjusted and unadjusted results
summary(adj.q)
summary(unadj.q)

# Extract deltas
deltas <- results[ ,5:6]
rownames(deltas) <- rownames(beta2)
colnames(deltas) <- c("Delta", "pvDelta")

# Get the custom filename
f.name <- paste("II.RefFreeEWAS/Data/", subtype, "_", stage, sep = "")

# write the adjusted, unadjusted, and the deltas values to file 
write.table(adjusted, file = paste(f.name, "_qvalues_adjusted.csv", sep = ""), sep = ",")
write.table(unadjusted, file = paste(f.name, "_qvalues_unadjusted.csv", sep = ""), sep = ",")
write.table(deltas, file = paste(f.name, "_delta.csv", sep = ""), sep = ",")

################################
# Extract Dimensionality Estimates
################################
string <- c(subtype, stage, returnList[[3]], n, bootstraps)
names(string) <- c("Subtype", "Stage", "Dimensionality", "n", "bootstraps")

write.table(t(string), file = paste("II.RefFreeEWAS/Data/dimension/", subtype, "_", stage, "_dim.txt", sep = ""), 
            row.names = F, sep = ",")
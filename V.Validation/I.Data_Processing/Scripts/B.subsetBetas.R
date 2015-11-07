#####################################################################
# ~~~~~~~~~~~~~~~~~~
# Tumor subtype and cell type independent DNA methylation alterations 
# associated with stage progression in invasive breast carcinoma 
# ~~~~~~~~~~~~~~~~~~
# Way, G., Johnson, K., Christensen, B. 2015

# The script will subset the full beta file into PAM50 subtype specific 
# subsets. Each subset will include all normal samples.
#####################################################################

################################
# Read in command line arguments
args <- commandArgs(trailingOnly = T)
args <- c("PAM50.RNAseq", "Basal", "Her2", "LumA", "LumB", "Normal")

# Large subset is the covariate information (column name) of the covariate file
Lsubset <- args[1]

# Small subset is the specific information in the covariate information 
# In this case: specific PAM50 subtypes
Ssubset <- args[2:length(args)]

################################
# Load Libraries
################################
library(minfi)

################################
# Load Data
################################
# Read in full beta values
beta <- read.table("I.Data_Processing/Data/TCGA_BRCA_Betas.tsv", row.names = 1, header = T, 
                   sep = "\t", stringsAsFactors = F, nrows = 450000, comment.char = "", 
                   colClasses = c("character", rep("numeric", 870)))

# Load covariate file
covariates <- read.table("I.Data_Processing/Files/BRCAtarget_covariates.csv", row.names = 1, 
                         header = T, sep = ",", stringsAsFactors = F)

# The colnames for the beta file have an "X" appended to the beginning of each basename, remove it
colnames(beta) <- substr(colnames(beta), 2, nchar(colnames(beta)))

# Get Normal samples
normals <- covariates[covariates$sample.type == "Solid Tissue Normal", ]

################################
# Subset Data and Write to File
################################
for (i in 1:length(Ssubset)) {
  # Separate models based on PAM50 subtype
  CovSubset <- covariates[covariates[ ,Lsubset] == Ssubset[i], ]
  # Attach all the normal samples to each model
  ToSubset <- unique(c(CovSubset$Basename, normals$Basename))
  # Subset the beta file
  betaSubset <- beta[ ,ToSubset]
  
  # Generate a file name and write to file
  file <- paste("I.Data_Processing/Data/BRCAmethSubset_", Lsubset, "_", Ssubset[i], ".tsv", sep = "")
  write.table(betaSubset, file = file, row.names = T, col.names = NA, sep = "\t")
}
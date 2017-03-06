#####################################################################
# ~~~~~~~~~~~~~~~~~~
# Tumor subtype and cell type independent DNA methylation alterations 
# associated with stage progression in invasive breast carcinoma 
# ~~~~~~~~~~~~~~~~~~
# Way, G., Johnson, K., Christensen, B. 2015
#
# Estimate dimensionality in entire dataset
#####################################################################
################################
# Load Library
################################
#install.packages('readr')
#install.packages('isva')
library(readr)
library(isva)

################################
# Load and Subset Data
################################
# load covariates
covariates <- read.table("I.Data_Processing/Files/BRCAtarget_covariates.csv", row.names = 1, 
                         header = T, sep = ",", stringsAsFactors = F)

# load Full betas
betas <- read_tsv("I.Data_Processing/Data/TCGA_BRCA_Betas.tsv")
rownames(betas) <- betas[["sample_id"]]
betas[["sample_id"]] <- NULL

# subset covariate file
rownames(covariates) <- covariates$Basename
#covariates <- covariates[-1:-19, ]

################################
# Run the dimensionality estimate
################################
# Full Model
mod <- model.matrix(~covariates$sample.type + covariates$age.Dx)

# Subset beta file
covariates <- covariates[-1:-19, ]
newBeta <- betas[ ,rownames(covariates)]
tempL <- as.matrix(newBeta[ , ])
tempR <- solve(t(mod)%*%mod)

# Obtain Bstar
tmpBstar <- tempL %*% mod %*% tempR

# Get Residual Matrix
tmpResid <- tempL-tmpBstar%*%t(mod)

# Perform EstDimRMT from isva package
RMT <- EstDimRMT(as.matrix(tmpResid), plot = F)
k <- RMT$dim

# Get results ready to write to file
string <- c("FullSet", "NA", k, ncol(newBeta), "NA")
names(string) <- c("Subtype", "Stage", "Dimensionality", "n", "bootstraps")
write.table(t(string), "II.RefFreeEWAS/Data/dimension/FullModel.txt", sep = ",")

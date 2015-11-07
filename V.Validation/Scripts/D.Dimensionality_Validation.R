#####################################################################
# ~~~~~~~~~~~~~~~~~~
# Tumor subtype and cell type independent DNA methylation alterations 
# associated with stage progression in invasive breast carcinoma 
# ~~~~~~~~~~~~~~~~~~
# Way, G., Johnson, K., Christensen, B. 2015

# This script will summarize dimensionality of the Validation set
#####################################################################

################################
# Load Libraries
################################
library(isva)
library(readr)
library(plyr)

################################
# Load and subset Data
################################
# Load smallest beta files (all of them contain all normal samples)
beta.file <- "V.Validation/Data/Validation_Betas.tsv"
beta2 <- read_tsv(beta.file)
rownames(beta2) <- beta2[ ,1]
beta2 <- beta2[ ,-1]
colnames(beta2) <- laply(colnames(beta2), function(x){paste(unlist(strsplit(x, "_"))[2], unlist(strsplit(x, "_"))[3], sep = "_")})

# load total covariate file
covariates <- read.table("V.Validation/Data/GSE60185_manifest.csv", row.names = 1, header = T, sep = ",", stringsAsFactors = F)
covariates$Basename <- paste(covariates$Sentrix_ID, covariates$Sentrix_Position, sep = "_")

# Subset the covariate file to the beta samples
covariates <- covariates[covariates$Basename %in% colnames(beta2), ]

# subset to only normals
normals <- covariates[covariates$Sample_Group == "control", ]
IDC <- covariates[covariates$Sample_Group == "IDC", ]

# Get a list of sample types
SampleList <- list("normals" = normals, "IDC" = IDC)

################################
# Get Dimensionality Estimates
################################
Dimensions <- matrix("", nrow = 2, ncol = 3)
for (i in 1:length(SampleList)) {
  # Subset Covariate File
  samps <- SampleList[[i]]
  
  # sample size
  n <- nrow(samps)
  
  # Construct model matrix 
  mod <- model.matrix(~1, data = samps)
  
  # Subset the original beta file 
  newBeta <- beta2[ ,samps$Basename]
  
  # Prepare for Random Matrix Theory Dimensionality
  tmpBstar <- as.matrix(newBeta[ , ]) %*% mod %*% solve(t(mod)%*%mod)
  # Get Residual Matrix
  tmpResid <- as.matrix(newBeta[ , ])-tmpBstar%*%t(mod)
  
  # Run EstDimRMT from package "isva"
  RMT <- EstDimRMT(as.matrix(tmpResid), plot = F)
  k <- RMT$dim
  
  Dimensions[i, ] <- c(n, k, names(SampleList)[i])
}
colnames(Dimensions) <- c("n = ", "Dimensions", "Validation")

# Write to file
write.table(Dimensions, file = "V.Validation/Tables/Dimension_Validation.csv", sep = ",", row.names = F)

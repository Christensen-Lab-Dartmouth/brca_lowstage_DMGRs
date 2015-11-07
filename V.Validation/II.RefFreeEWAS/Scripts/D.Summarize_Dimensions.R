#####################################################################
# ~~~~~~~~~~~~~~~~~~
# Tumor subtype and cell type independent DNA methylation alterations 
# associated with stage progression in invasive breast carcinoma 
# ~~~~~~~~~~~~~~~~~~
# Way, G., Johnson, K., Christensen, B. 2015
#
# Compile a summary table of dimenstionality for the RefFree Analysis
# It will also do an additional analysis on only Normal Samples
#####################################################################

################################
# Load Libraries
################################
library(isva)

################################
# Extract dimension and sample size information
################################
# This piece will load all of the .txt files that hold sample sizes and dimensionality for each restricted subset
loc <- "II.RefFreeEWAS/Data/dimension/"
filesFull <- list.files(loc)
files <- filesFull[!grepl("Full", filesFull)]
dimension <- c()
for(i in 1:length(files)) {
  tmp <- read.table(paste(loc, files[i], sep = ""), header = T, sep = ",", stringsAsFactors = F)
  dimension <- rbind(dimension, tmp)
}

################################
# Run analysis on normal tissue only
################################
# Load smallest beta files (all of them contain all normal samples)
beta.file <- "I.Data_Processing/Data/BRCAmethSubset_PAM50.RNAseq_Normal.tsv"
beta2 <- read.table(beta.file, row.names = 1, header = T, sep = "\t", stringsAsFactors = F, 
                    nrows = 400000, comment.char = "")

# load total covariate file
covariates <- read.table("I.Data_Processing/Files/BRCAtarget_covariates.csv", row.names = 1, header = T, 
                         sep = ",", stringsAsFactors = F)

# The colnames for the beta file have an "X" appended to the beginning of each basename, remove it
colnames(beta2) <- substr(colnames(beta2), 2, nchar(colnames(beta2)))
rownames(covariates) <- covariates$Basename

# Subset the covariate data to only the samples in the beta file
covariates <- covariates[intersect(rownames(covariates), colnames(beta2)), ]

# subset to only normals
normals <- covariates[covariates$sample.type == "Solid Tissue Normal", ]

# sample size
n <- nrow(normals)

# Construct model matrix only correcting for age
mod <- model.matrix(~normals$age.Dx)

# Subset the original beta file to only normal samples
newBeta <- beta2[ ,rownames(normals)]

# Prepare for Random Matrix Theory Dimensionality
tmpBstar <- as.matrix(newBeta[ , ]) %*% mod %*% solve(t(mod)%*%mod)
# Get Residual Matrix
tmpResid <- as.matrix(newBeta[ , ])-tmpBstar%*%t(mod)

# Run EstDimRMT from package "isva"
RMT <- EstDimRMT(as.matrix(tmpResid), plot = F)
k <- RMT$dim

# Compile Info and append to RefFree Summary table
string <- c("Normal Tissue", "Normal", k, n, "No RefFree Run")
names(string) <- c("Subtype", "Stage", "Dimensionality", "n", "bootstraps")

################################
# Combine Data Together to Write
################################
# Load full dimensionality data
all <- read.table(paste(loc, setdiff(filesFull, files), sep = ""), header = T, sep = ",", stringsAsFactors = F)
dimension <- rbind(dimension, string, all)

write.table(dimension, file = "II.RefFreeEWAS/Tables/RefFreeSummary.csv", sep = ",", row.names = T, col.names = NA)

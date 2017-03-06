#########################
# Much of this specific file was manually processed
#########################

# Set up the workspace
# install.packages('corrgram')
# source("https://bioconductor.org/biocLite.R")
# biocLite("qvalue")
require(data.table)
require(limma)
require(qvalue)

rm(list = ls())

# Open and load the data
dir <- '/Users/alexandertitus/Documents/brca_lowstage_DMGRs'
setwd(dir)

subtype <- "Basal"
#subtype <- "Her2"
#subtype <- "LumA"
#subtype <- "LumB"
#subtype <- "Normal"

stage <- "low"
#stage <- "high"

file_name <- paste('II.RefFreeEWAS/Data/RefFreeEWAS/', subtype, '_', stage, '_RefFree2.0List_05Oct2016.RData', sep = '')
load(file_name)

head(RefFree_DMGR_Boots)
which.min(apply(RefFree_DMGR_Boots[-1, ], 2, median)) 
which.min(apply(RefFree_DMGR_Boots[-1, ], 2, mean)) 
which.min(apply(RefFree_DMGR_Boots[-1, ], 2, mean, trim=0.05)) 
which.min(apply(RefFree_DMGR_Boots[-1, ], 2, mean, trim=0.15)) 
which.min(apply(RefFree_DMGR_Boots[-1, ], 2, mean, trim=0.25)) 
which.min(apply(RefFree_DMGR_Boots[-1, ], 2, mean, trim=0.35)) 
which.min(apply(RefFree_DMGR_Boots[-1, ], 2, mean, trim=0.45)) 
k <- 5 # Basal low
#k <- 3 # Basal high
#k <- 5 # Her2 low
#k <- 4 # Her2 high
#k <- 9 # LumA low
#k <- 5 # LumA high
#k <- 7 # LumB low
#k <- 4 # LumB high
#k <- 8 # Normal low
#k <- 7 # Normal high

# Inspect the results
head(DMGR_RefFree_Array2[[k-1]]$Mu)
head(DMGR_RefFree_Array2[[k-1]]$Omega)

png(paste('/Users/alexandertitus/Documents/brca_lowstage_DMGRs/II.RefFreeEWAS/Figures/k_distribution/', 
          subtype, '_', stage, '_', 'k_estimate.png', sep = ''))
boxplot(DMGR_RefFree_Array2[[k-1]]$Omega)
dev.off()

# Update the covariate file specific to this subset
DMGR_RefFree_Array2[[k-1]]$Omega[DMGR_RefFree_Array2[[k-1]]$Omega < 0] = 0
covariates <- read.table("I.Data_Processing/Files/BRCAtarget_covariates.csv", row.names = 1, header = T, sep = ",", stringsAsFactors = F)
inter <- intersect(rownames(DMGR_RefFree_Array2[[k-1]]$Omega), covariates$Basename)
covariates <- covariates[covariates$Basename %in% inter,]
DMGR_cell_proportions <- cbind(covariates, DMGR_RefFree_Array2[[k-1]]$Omega)

DMGR_om <- DMGR_RefFree_Array2[[k-1]]$Omega

# Drop one cell-type to avoid multi-collinearity
DMGR_om <- DMGR_om[, 1:(k-1)]

XX <- model.matrix(~sample.type + age.Dx, data = covariates)
rownames(XX) <- covariates$Basename

# Load betas
beta.file <- paste("I.Data_Processing/Data/BRCAmethSubset_PAM50.RNAseq_", subtype, ".tsv", sep = '')
betas <- data.frame(fread(beta.file), row.names=1)
colnames(betas) <- substring(colnames(betas), 2, length(colnames(betas)))
betas = betas[ ,order(colnames(betas), decreasing=T)]
betas_dmgr <- data.matrix(betas)

# Convert to M-values
betas_dmgrM <- ifelse(betas_dmgr>=1, 1-1E-6, ifelse(betas_dmgr<=0, 1E-6, betas_dmgr))
betas_dmgrM <- log(betas_dmgrM)-log(1-betas_dmgrM)
betas_dmgrM2 <- betas_dmgrM[, colnames(betas_dmgrM) %in% rownames(XX)]
XX = XX[order(rownames(XX), decreasing=T), ]
all(colnames(betas_dmgrM2) == rownames(XX))

# Limma models
lf_Null <-  eBayes(lmFit(betas_dmgrM2, XX))
lf_Omega <- eBayes(lmFit(betas_dmgrM2, cbind(XX, DMGR_om)))

# analysis of q-values
adj.q <- qvalue(lf_Omega$p.value[,2])
unadj.q <- qvalue(lf_Null$p.value[,2])

# Prepare dataframe for export
adjusted <- cbind(adj.q$pvalues, adj.q$qvalues, lf_Omega$coefficients[, 2])
unadjusted <- cbind(unadj.q$pvalues, unadj.q$qvalues, lf_Null$coefficients[, 2])

# Name the rows and columns for each model
rownames(adjusted) <- rownames(unadjusted) <- rownames(betas)
colnames(adjusted) <- colnames(unadjusted) <- c("pvalues", "qvalues", "beta")

# print out summary of the adjusted and unadjusted results
summary(adj.q)
summary(unadj.q)

# Extract deltas 
deltas <- lf_Omega$coefficients[, 2] - lf_Null$coefficients[, 2]
deltas <- cbind(deltas, rep(1, length(deltas)))

rownames(deltas) <- rownames(betas)
colnames(deltas) <- c("Delta", "temp")

# Get the custom filename
f.name <- paste("II.RefFreeEWAS/Data/", subtype, "_", stage, sep = "")

# write the adjusted, unadjusted, and the deltas values to file 
write.table(adjusted, file = paste(f.name, "_qvalues_adjusted.csv", sep = ""), sep = ",")
write.table(unadjusted, file = paste(f.name, "_qvalues_unadjusted.csv", sep = ""), sep = ",")
write.table(deltas, file = paste(f.name, "_delta.csv", sep = ""), sep = ",")


################################
# Extract Dimensionality Estimates
################################
n <- nrow(covariates[covariates$sample.type == "Primary Tumor", ])
bootstraps <- 10
string <- c(subtype, stage, k, n, bootstraps)
names(string) <- c("Subtype", "Stage", "Dimensionality", "n", "bootstraps")

write.table(t(string), file = paste("II.RefFreeEWAS/Data/dimension/", subtype, "_", stage, "_dim.txt", sep = ""), 
            row.names = F, sep = ",")
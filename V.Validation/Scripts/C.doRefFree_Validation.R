#####################################################################
# ~~~~~~~~~~~~~~~~~~
# Tumor subtype and cell type independent DNA methylation alterations 
# associated with stage progression in invasive breast carcinoma 
# ~~~~~~~~~~~~~~~~~~
# Way, G., Johnson, K., Christensen, B. 2015
#
# The script will apply RefFreeEWAS to the Full Validation Set
#####################################################################

################################
# Read in command line arguments
################################
bootstraps <- 10 # commandArgs(trailingOnly = T)[1]
set.seed(123)

################################
# Load Libraries
################################
#install.packages('RefFreeEWAS')
#source("https://bioconductor.org/biocLite.R")
#biocLite("qvalue")
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
# Load the Validation Betas
beta2 <- read.table("V.Validation/Data/Validation_Betas.tsv", row.names = 1, header = T, 
                    sep = "\t", stringsAsFactors = F, nrows = 400000, comment.char = "")

# Load covariate file that holds information on Tumor vs. Normal
covariates <- read.table("V.Validation/Tables/ValidationSampleID.txt", row.names = 1, 
                         header = T, sep = ",", stringsAsFactors = F)

# We only need tumor vs. control info
cov <- data.frame(covariates$Sample_Group)
rownames(cov) <- rownames(covariates)
colnames(cov) <- "sample.type"

# How many samples
n <- length(cov)

################################
# Run customRefFree
################################
cat("bootstraps: ", bootstraps, "\n")
file = "Validation_Set"
cat(file, "\n")


# subset betas
newBeta <- beta2[ ,rownames(cov)]
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
     file=paste("V.Validation/Data/Validation_RefFree2.0List_Dec2016.RData", sep = ""), 
     compress=TRUE)


file_name <- 'V.Validation/Data/Validation_RefFree2.0List_Dec2016.RData'
load(file_name)

head(RefFree_DMGR_Boots)
which.min(apply(RefFree_DMGR_Boots[-1, ], 2, median)) 
which.min(apply(RefFree_DMGR_Boots[-1, ], 2, mean)) 
which.min(apply(RefFree_DMGR_Boots[-1, ], 2, mean, trim=0.05)) 
which.min(apply(RefFree_DMGR_Boots[-1, ], 2, mean, trim=0.15)) 
which.min(apply(RefFree_DMGR_Boots[-1, ], 2, mean, trim=0.25)) 
which.min(apply(RefFree_DMGR_Boots[-1, ], 2, mean, trim=0.35)) 
which.min(apply(RefFree_DMGR_Boots[-1, ], 2, mean, trim=0.45)) 
k <- 8 # Validation  estimates

# Inspect the results
head(DMGR_RefFree_Array2[[k-1]]$Mu)
head(DMGR_RefFree_Array2[[k-1]]$Omega)

png(paste('V.Validation/Figures/k_estimate.png', sep = ''))
  boxplot(DMGR_RefFree_Array2[[k-1]]$Omega, main='Distribution of putative cell-type estimates')
dev.off()

# Update the covariate file specific to this subset
DMGR_RefFree_Array2[[k-1]]$Omega[DMGR_RefFree_Array2[[k-1]]$Omega < 0] = 0
inter <- intersect(rownames(DMGR_RefFree_Array2[[k-1]]$Omega), rownames(cov))
# cov <- cov[rownames(cov) %in% inter,]
DMGR_cell_proportions <- cbind(cov, DMGR_RefFree_Array2[[k-1]]$Omega)

DMGR_om <- DMGR_RefFree_Array2[[k-1]]$Omega

# Drop one cell-type to avoid multi-collinearity
DMGR_om <- DMGR_om[, 1:(k-1)]

XX <- model.matrix(~sample.type, data = cov)

# Load betas
require(data.table)
beta.file <- paste("V.Validation/Data/Validation_Betas.tsv", sep = '')
betas <- data.frame(fread(beta.file), row.names=1)
betas = betas[ ,order(colnames(betas), decreasing=T)]
betas_dmgr <- data.matrix(betas)

# Convert to M-values
betas_dmgrM <- ifelse(betas_dmgr>=1, 1-1E-6, ifelse(betas_dmgr<=0, 1E-6, betas_dmgr))
betas_dmgrM <- log(betas_dmgrM)-log(1-betas_dmgrM)
betas_dmgrM2 <- betas_dmgrM[, colnames(betas_dmgrM) %in% rownames(XX)]
XX = XX[order(rownames(XX), decreasing=T), ]
all(colnames(betas_dmgrM2) == rownames(XX))

# Limma models
require(limma)
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

# Extract deltas - MORE WORK TO DO HERE
deltas <- lf_Omega$coefficients[, 2] - lf_Null$coefficients[, 2]
deltas <- cbind(deltas, rep(1, length(deltas)))

#seDelta <- apply(RefFree_DMGR_Boots[, 4]-rfb[ , , "B", ], 1:2, sd)
#pvDelta <- 2 * pt(-abs(Delta)/seDelta[ ,2], denDegFree)
rownames(deltas) <- rownames(betas)
colnames(deltas) <- c("Delta", "temp")

# Get the custom filename
f.name <- paste("V.Validation/Data/Validation", sep = "")

# write the adjusted, unadjusted, and the deltas values to file 
write.table(adjusted, file = paste(f.name, "_qvalues_adjusted.csv", sep = ""), sep = ",")
write.table(unadjusted, file = paste(f.name, "_qvalues_unadjusted.csv", sep = ""), sep = ",")
write.table(deltas, file = paste(f.name, "_delta.csv", sep = ""), sep = ",")

################################
# Visualize RefFreeEWAS results
################################
# Draw Volcano Plots
# Read in data
adjusted <- read.table("V.Validation/Data/Validation_qvalues_adjusted.csv", sep = ",")
unadjusted <- read.table("V.Validation/Data/Validation_qvalues_unadjusted.csv", sep = ",")
deltas <- read.table("V.Validation/Data/Validation_delta.csv", sep = ",")

# Append adjusted Q values to the results file
# analysis of q-values
qtmp <- qvalue(adjusted$pvalues, fdr.level = 0.01)
qtmp_high <- qvalue(adjusted$pvalues, fdr.level = 0.05)

# Get Q value cutoffs to draw lines in plots
tmp1 <- adjusted[qtmp$significant == T, ]
tmp2 <- adjusted[qtmp_high$significant == T, ]

# Of these cutoffs, what is the lowest p value?
qcut1 <- min(-log10(as.numeric(paste(tmp1[,3])))) 
qcut2 <- min(-log10(as.numeric(paste(tmp2[,3]))))

# Find max and min for plotting margins
# Find max and min for plotting margins
maxX <- max(c(max(as.numeric(paste(adjusted[ ,3]))),as.numeric(paste(max(unadjusted[ ,3])))))
maxY <- max(c(-log10(as.numeric(paste(unadjusted[ ,1])))[!is.infinite(-log10(as.numeric(paste(unadjusted[ ,1]))))], 
              -log10(as.numeric(paste(adjusted[ ,1])))[!is.infinite(-log10(as.numeric(paste(adjusted[ ,1]))))]))

# Plotting margin for adjusted volcano plot
maxY1 <- max(-log10(as.numeric(paste(adjusted[ ,1])))[!is.infinite(-log10(as.numeric(paste(adjusted[ ,1]))))])

# Store the plots to then arrange them on a grid to save as png

# Unadjusted Plot
un <- ggplot(unadjusted, aes(as.numeric(paste(unadjusted[ ,3])), -log10(as.numeric(paste(unadjusted[ ,1]))))) + 
  geom_point(aes(colour = deltas[ ,1]), size = 9) + 
  scale_color_gradient2(low = "blue", mid="grey", high = "red") + 
  labs(list(x = "Beta Coefficient", y = "-log10 p Value", title = paste("Reference Free Undjusted\nValidation"), color = "Delta"), size = 18) + 
  geom_hline(yintercept = qcut1, color = "red", lintetype = "longdash") + 
  geom_hline(yintercept = qcut2, acolor = "black", lintetype = "longdash") +
  xlim((-1 * maxX), maxX) + ylim(0, maxY) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(), 
         axis.text = element_text(size = rel(1.4), color = "black"),
        axis.title = element_text(size = rel(1.6)), plot.title = element_text(size = rel(2)),
        legend.text = element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.3)),
        axis.ticks = element_line(size = rel(1.4), color = "black"),
        axis.line = element_line(size = rel(1.4), color = "black"))

# Adjusted plot
ad <- ggplot(adjusted, aes(as.numeric(paste(adjusted[ ,3])), 
                           -log10(as.numeric(paste(adjusted[ ,1]))))) + 
  geom_point(aes(colour = deltas[ ,1]), size = 9) + 
  scale_color_gradient2(low = "blue", mid = "grey", high = "red") + 
  labs(list(x = "Beta Coefficient", y = "-log10 p Value", title = paste("Reference Free Adjusted\nValidation"), color = "Delta")) + 
  geom_hline(yintercept = qcut1, color = "red", linetype = "dashed", size = 1.8) + 
  geom_hline(yintercept = qcut2, color = "black", linetype = "dashed", size = 1.8) +
  xlim((-1 * maxX), maxX) + ylim(0, maxY) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(), 
        axis.text = element_text(size = rel(1.4), color = "black"),
        axis.title = element_text(size = rel(1.6)), plot.title = element_text(size = rel(2)),
        legend.text = element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.3)),
        axis.ticks = element_line(size = rel(1.4), color = "black"),
        axis.line = element_line(size = rel(1.4), color = "black"))

# Adjusted plot with different max y value
ad1 <- ggplot(adjusted, aes(as.numeric(paste(adjusted[,3])), 
                            -log10(as.numeric(paste(adjusted[,1]))))) + 
  geom_point(aes(colour = deltas[ ,1]), size = 9) + 
  scale_color_gradient2(low = "blue", mid="grey", high = "red") + 
  labs(list(x = "Beta Coefficient", y = "-log10 p Value", title = paste("Reference Free Adjusted (Resized)\nValidation"), color = "Delta")) + 
  geom_hline(yintercept = qcut1, color = "red", linetype = "dashed", size = 1.8) + 
  geom_hline(yintercept = qcut2, color = "black", linetype = "dashed", size = 1.8) +
  xlim((-1 * maxX), maxX) + ylim(0, maxY1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(), 
        axis.text = element_text(size = rel(1.4), color = "black"),
        axis.title = element_text(size = rel(1.6)), plot.title = element_text(size = rel(2)),
        legend.text = element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.3)),
        axis.ticks = element_line(size = rel(1.4), color = "black"),
        axis.line = element_line(size = rel(1.4), color = "black"))

png(paste("V.Validation/Figures/", file, "_", bootstraps, "_volcano.png", sep = ""), width = 900, height = 600)
grid.arrange(un, ad, ad1, ncol = 3, nrow = 1)
dev.off()


################################
# Extract Dimensionality Estimates
################################
string <- c("Validation", "", k, n, bootstraps)
names(string) <- c("Subtype", "Stage", "Dimensionality", "n", "bootstraps")

write.table(t(string), file = paste("II.RefFreeEWAS/Data/dimension/Validation_set_dim.txt", sep = ""), row.names = F, sep = ",")

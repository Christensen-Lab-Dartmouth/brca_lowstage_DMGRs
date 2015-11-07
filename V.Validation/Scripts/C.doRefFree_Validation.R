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
bootstraps <- commandArgs(trailingOnly = T)[1]
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
returnList <- customRefFree(covariates = cov, betas = beta2, age = F, bootstraps = bootstraps)

################################
# Visualize RefFreeEWAS results
################################
# Draw Volcano Plots
results <- data.frame(returnList[[2]])

# Append adjusted Q values to the results file
qtmp <- qvalue(results$pvAdj, fdr.level = 0.01)
qtmp_high <- qvalue(results$pvAdj, fdr.level = 0.05)

# Get Q value cutoffs to draw lines in plots
tmp1 <- results[qtmp$significant == T, ]
tmp2 <- results[qtmp_high$significant == T, ]

# Of these cutoffs, what is the lowest p value?
qcut1 <- min(-log10(as.numeric(paste(tmp1[,3])))) 
qcut2 <- min(-log10(as.numeric(paste(tmp2[,3]))))

# Find max and min for plotting margins
maxX <- max(c(max(as.numeric(paste(results[ ,1]))),as.numeric(paste(max(results[ ,2])))))
maxY <- max(c(-log10(as.numeric(paste(results[ ,3])))[!is.infinite(-log10(as.numeric(paste(results[ ,3]))))], 
              -log10(as.numeric(paste(results[ ,4])))[!is.infinite(-log10(as.numeric(paste(results[ ,4]))))]))

# Plotting margin for adjusted volcano plot
maxY1 <- max(-log10(as.numeric(paste(results[ ,3])))[!is.infinite(-log10(as.numeric(paste(results[ ,3]))))])

# Store the plots to then arrange them on a grid to save as png
# Unadjusted Plot
un <- ggplot(results, aes(as.numeric(paste(results[ ,2])), -log10(as.numeric(paste(results[ ,4]))))) + 
  geom_point(aes(colour = results[ ,colnames(results)[grepl("Delta", colnames(results))][1]])) + 
  scale_color_gradient2(low = "blue", mid="grey", high = "red") + 
  labs(list(x = "Beta Coefficient", y = "-log10 p Value", title = paste("Reference Free Undjusted\nValidation"), color = "Delta"), size = 18) + 
  geom_hline(yintercept = qcut1, color = "red", lintetype = "longdash") + 
  geom_hline(yintercept = qcut2, acolor = "black", lintetype = "longdash") +
  xlim((-1 * maxX), maxX) + ylim(0, maxY) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(), 
        axis.line = element_line(), axis.text = element_text(size = rel(1.4), color = "black"),
        axis.title = element_text(size = rel(1.6)), plot.title = element_text(size = rel(2)),
        legend.text = element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.3)),
        axis.ticks = element_line(size = rel(1.4), color = "black"),
        axis.line = element_line(size = rel(1.4), color = "black"))

# Adjusted plot
ad <- ggplot(results, aes(as.numeric(paste(results[ ,1])), -log10(as.numeric(paste(results[ ,3]))))) + 
  geom_point(aes(colour = results[,colnames(results)[grepl("Delta", colnames(results))][1]])) + 
  scale_color_gradient2(low = "blue", mid="grey", high = "red") + 
  labs(list(x = "Beta Coefficient", y = "-log10 p Value", title = paste("Reference Free Adjusted\nValidation"), color = "Delta"), size = 18) + 
  geom_hline(yintercept = qcut1, color = "red", lintetype = "longdash") +
  geom_hline(yintercept = qcut2, acolor = "black", lintetype = "longdash") +
  xlim((-1 * maxX), maxX) + ylim(0, maxY) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(), 
        axis.line = element_line(), axis.text = element_text(size = rel(1.4), color = "black"),
        axis.title = element_text(size = rel(1.6)), plot.title = element_text(size = rel(2)),
        legend.text = element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.3)),
        axis.ticks = element_line(size = rel(1.4), color = "black"),
        axis.line = element_line(size = rel(1.4), color = "black"))
  
# Adjusted plot with different max y value
ad1 <- ggplot(results, aes(as.numeric(paste(results[ ,1])), -log10(as.numeric(paste(results[ ,3]))))) + 
  geom_point(aes(colour = results[,colnames(results)[grepl("Delta", colnames(results))][1]])) + 
  scale_color_gradient2(low = "blue", mid="grey", high = "red") + 
  labs(list(x = "Beta Coefficient", y = "-log10 p Value", title = paste("Reference Free Adjusted (Resized)\nValidation"), color = "Delta"), size = 18) + 
  geom_hline(yintercept = qcut1, color = "red", lintetype = "longdash") + 
  geom_hline(yintercept = qcut2, acolor = "black", lintetype = "longdash") +
  xlim((-1 * maxX), maxX) + ylim(0, maxY1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(), 
        axis.line = element_line(), axis.text = element_text(size = rel(1.4), color = "black"),
        axis.title = element_text(size = rel(1.6)), plot.title = element_text(size = rel(2)),
        legend.text = element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.3)),
        axis.ticks = element_line(size = rel(1.4), color = "black"),
        axis.line = element_line(size = rel(1.4), color = "black"))

png(paste("V.Validation/Figures/", file, "_", bootstraps, "_volcano.png", sep = ""), width = 900, height = 600)
grid.arrange(un, ad, ad1, ncol = 3, nrow = 1)
dev.off()

################################
# Extract q value summary
################################
adj.q <- qvalue(results[ ,3])
unadj.q <- qvalue(results[ ,4])

# Prepare dataframe for export
adjusted <- cbind(adj.q$pvalues, adj.q$qvalues, results[ ,1])
unadjusted <- cbind(unadj.q$pvalues, unadj.q$qvalues, results[ ,2])

rownames(adjusted) <- rownames(beta2)
colnames(adjusted) <- c("pvalues", "qvalues", "beta")
rownames(unadjusted) <- rownames(beta2)
colnames(unadjusted) <- c("pvalues", "qvalues", "beta")

# print out summary of the adjusted and unadjusted results
summary(adj.q)
summary(unadj.q)

# Extract deltas
deltas <- results[ ,5:6]
rownames(deltas) <- rownames(beta2)
colnames(deltas) <- c("Delta", "pvDelta")

f.name <- "V.Validation/Data/Validation_set"

write.table(adjusted, file = paste(f.name, "_qvalues_adjusted.csv", sep = ""), sep = ",")
write.table(unadjusted, file = paste(f.name, "_qvalues_unadjusted.csv", sep = ""), sep = ",")
write.table(deltas, file = paste(f.name, "_delta.csv", sep = ""), sep = ",")

################################
# Extract Dimensionality Estimates
################################
string <- c("Validation", "", returnList[[3]], n, bootstraps)
names(string) <- c("Subtype", "Stage", "Dimensionality", "n", "bootstraps")

write.table(t(string), file = paste("II.RefFreeEWAS/Data/dimension/Validation_set_dim.txt", sep = ""), row.names = F, sep = ",")

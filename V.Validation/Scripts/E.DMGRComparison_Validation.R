#####################################################################
# ~~~~~~~~~~~~~~~~~~
# Tumor subtype and cell type independent DNA methylation alterations 
# associated with stage progression in invasive breast carcinoma 
# ~~~~~~~~~~~~~~~~~~
# Way, G., Johnson, K., Christensen, B. 2015

# Get q value for all DMGRs in validation set; plot heatmap and validation figure
#####################################################################

################################
# Load Libraries
################################
library(readr)
library(Hmisc)
source("III.DMGR_analysis/Scripts/Functions/make_heatmaps.R") #this will load custom CGHeatmap and heatmap.3 function

source("III.DMGR_analysis/Scripts/Functions/DMcgs_functions.R")

################################
# Load Function
################################
# Extract information for the given DMGR to use to in summary plot
find_min_DMGR <- function (data, DMGR) {
  subset <- data[grepl(DMGR, data$Gene_Region), ]
  
  # Find info about the minimum cpg qvalue
  info <- subset[nrow(subset), ]
  
  # What is the CpG of interest
  cpg <- info$TargetID
  
  # What index of the cpg
  Index <- as.numeric(paste(rownames(info)))
  
  # What is the q value of the cpg
  qval <- info$qvalues
  
  # Return the pertinent data
  return(c(cpg, Index, qval))
}

################################
# Load Data
################################
# Load adjusted data and deltas
adjustedQ <- read.table("V.Validation/Data/Validation_set_qvalues_adjusted.csv", sep = ",", stringsAsFactors = F)
Val.Deltas <- read.table("V.Validation/Data/Validation_set_delta.csv", sep = ",", stringsAsFactors = F)

# Combine them into one variable
adjustedQ <- cbind(adjustedQ, Val.Deltas)

# Load extended annotation file. This file does not ignore multiple genes for a single cpg
annotation <- read.table("I.Data_Processing/Files/Expanded_annotationfile.csv", stringsAsFactors = F, row.names = 1, 
                         header = T, sep = ",", nrows = 1000000, comment.char = "")

# Append a column of only the first gene and first region to the annotation file
annoFirst <- apply(annotation, 1, function (x) {splitGene(x[12])[1]})
annoReg <- apply(annotation, 1, function (x) {splitGene(x[14])[1]})
annotation <- cbind(annotation, paste(annoFirst, annoReg))
colnames(annotation)[ncol(annotation)] <- "Gene_Region"

# subset annotation file
annotation <- annotation[annotation$TargetID %in% rownames(adjustedQ), ]

# Load Test Set low stage overlapping DMGRs
test_set <- read_csv("III.DMGR_analysis/Tables/commonLowStageOverlaps_FullAnnotation_extended.csv")
rownames(test_set) <- test_set[ ,1]
test_set <- test_set[ ,-1]

rm(Val.Deltas, annoFirst, annoReg)
################################
# Combine Data
################################
Combine <- cbind(annotation, adjustedQ[match(annotation$TargetID, rownames(adjustedQ)), ])

# Order Combined Data by increasing q value
CombineOrdered <- Combine[order(Combine$qvalues, decreasing = T), ]

# Remove cpgs that don't map to genes
CombineOrdered <- CombineOrdered[CombineOrdered$UCSC_RefGene_Name != "", ]
# Change the rownames to indicate the index of the cpg
rownames(CombineOrdered) <- 1:nrow(CombineOrdered)

# Extract information from test set
Indices <- CpGs <- Q <- c()
for (i in 1:nrow(test_set)) {
  result <- find_min_DMGR(CombineOrdered, rownames(test_set)[i])
  cat(paste(rownames(test_set)[i], ":", sep = ""), result, "\n")
  Indices <- c(Indices, as.numeric(paste(result[2])))
  CpGs <- c(CpGs, result[1])
  Q <- c(Q, as.numeric(paste(result[3])))
}

# The CombineOrdered file is plotReady
plotReady <- CombineOrdered

################################
# Load and Process Heatmap Data
################################
# Load Betas
Betas <- read_tsv("V.Validation/Data/Validation_Betas.tsv")
rownames(Betas) <- Betas[ ,1]
Betas <- Betas[ ,-1]

# Load Annotation File with Gene Regions
annotationGR <- read.table("I.Data_Processing/Files/Expanded_annotationfile.csv", stringsAsFactors = F, row.names = 1, 
                           header = T, sep = ",", nrows = 1000000, comment.char = "")

# Load Annotation File with Map Info
annotationMap <- read_csv("I.Data_Processing/Files/HumanMethylation450_15017482_v.1.1.csv")

# Load covariate file
covariates <- read.table("V.Validation/Tables/ValidationSampleID.txt", sep = ",", row.names = 1, header = T, 
                         stringsAsFactors = F)
covariates <- covariates[colnames(Betas), ]

################################
# Load Common Overlaps
################################
CommonOverlaps <- read.csv("III.DMGR_analysis/Tables/commonLowStageOverlaps_FullAnnotation_extended.csv", sep = ",", 
                           row.names = 1, header = T, stringsAsFactors = F)

################################
# Run Custom HeatMap Function
################################
# Extract Genes and Gene Regions
Genes <- c()
Regions <- c()
for (i in 1:nrow(CommonOverlaps)) {
  Genes <- c(Genes, unlist(strsplit(rownames(CommonOverlaps)[i], " "))[1])
  Regions <- c(Regions, unlist(strsplit(rownames(CommonOverlaps)[i], " "))[2])
}

################################
# Total CpG Heatmap
################################
# Create Large Heatmap with all CpGs
lowcgs <- c()
regions <- c()
for (i in 1:length(Genes)) {
  cg <- ExtractCommonCGs(Genes[i], CommonOverlaps = CommonOverlaps)
  lowcgs <- c(lowcgs, unique(cg))
  tmpregion <- c(rep(Regions[i], length(unique(cg))))
  regions <- c(regions, tmpregion)
}

################################
# Plot Figures (These will be two panels of the same figure)
################################
# Output heatmap visualization
png("V.Validation/Figures/Heatmap_Validation.png", height = 600, width = 800)
CGHeatmap("Validation Set", lowcgs, regions, annotationGR, annotationMap, Covariates = covariates, Betas, Validation = T)
dev.off()

# Output Line Graph visualization
png("V.Validation/Figures/LineGraph_Validation.png", height = 400, width = 400)
plot(x = 1:nrow(plotReady), y = -log10(as.numeric(paste(plotReady$qvalues))), 
     xlab = "Q value rank", ylab = "-log10 Q value", col = "black", cex = 0.7, 
     pch = 16, bty = "n")

# Extend the length of the x axis
axis(1, at=c(0, 1e5, 2e5, 3e5, 4e5, 5e5, 6e5), col = "black")

# Add dashed lines indicating cutoff
abline(h = -log10(0.01), col = "grey", lwd = 2, lty = "dashed")
text(x = 1e+05, y = -log10(0.01)+0.3, "Q = 0.01 cutoff", col = "grey")
abline(h = -log10(0.05), col = "grey", lwd = 2, lty = "dashed")
text(x = 1e+05, y = -log10(0.05)+0.3, "Q = 0.05 cutoff", col = "grey")


# Add points for the test set overlap
points(x = as.numeric(paste(Indices)), y = -log10(as.numeric(paste(Q))), pch = 16, 
       col = "red", cex = 1)

# Add text in different locations to avoid overlap
text(x = as.numeric(paste(Indices))[c(1, 5, 8, 11:12)], 
     y = -log10(as.numeric(paste(Q)))[c(1, 5, 8, 11:12)], 
     rownames(test_set)[c(1, 5, 8, 11:12)], 
     cex = .7, pos = 2, offset = 2.1)
text(x = as.numeric(paste(Indices))[2], y = -log10(as.numeric(paste(Q)))[2] + 0.1, 
     rownames(test_set)[2], cex = .7, pos = 2, offset = 1.25)
text(x = as.numeric(paste(Indices))[7], y = -log10(as.numeric(paste(Q)))[7] + 0.1, 
     rownames(test_set)[7], cex = .7, pos = 2, offset = 1.25)
text(x = as.numeric(paste(Indices))[10], y = -log10(as.numeric(paste(Q)))[10] + 0.15, 
     rownames(test_set)[10], cex = .7, pos = 2, offset = 1.25)
text(x = as.numeric(paste(Indices))[3], y = -log10(as.numeric(paste(Q)))[3] + .025, 
     rownames(test_set)[3], cex = .7, pos = 2, offset = .25)
text(x = as.numeric(paste(Indices))[6], y = -log10(as.numeric(paste(Q)))[6] + 0.15, 
     rownames(test_set)[6], cex = .7, pos = 2, offset = .25)
text(x = as.numeric(paste(Indices))[9], y = -log10(as.numeric(paste(Q)))[9] + 0.15, 
     rownames(test_set)[9], cex = .7, pos = 2, offset = 4)
text(x = as.numeric(paste(Indices))[4], y = -log10(as.numeric(paste(Q)))[4] + .025, 
     rownames(test_set)[4], cex = .7, pos = 2, offset = 4)

# Find the top couple of points
# topPoints <- CombineOrdered[(nrow(CombineOrdered) - 10):nrow(CombineOrdered), ]
# 
# points(x = as.numeric(paste(rownames(topPoints))), 
#        y = -log10(as.numeric(paste(topPoints$qvalues))), 
#        pch = 16, col = "blue", cex = 1)
# text(x = as.numeric(paste(rownames(topPoints)))[c(1, 7, 8, 11)], 
#      y = -log10(as.numeric(paste(topPoints$qvalues)))[c(1, 7, 8, 11)], 
#      topPoints$Gene_Region[c(1, 7, 8, 11)], cex = .7, pos = 2, offset = 2)
# text(x = as.numeric(paste(rownames(topPoints)))[c(3, 9)], 
#      y = -log10(as.numeric(paste(topPoints$qvalues)))[c(3,9)], 
#      topPoints$Gene_Region[c(3, 9)], cex = .7, pos = 2, offset = .3)
# text(x = as.numeric(paste(rownames(topPoints)))[c(5, 10)], 
#      y = -log10(as.numeric(paste(topPoints$qvalues)))[c(5, 10)], 
#      topPoints$Gene_Region[c(5, 10)], cex = .7, pos = 2, offset = 4.5)

legend("topleft", inset = 0.05, title = "DMGRs", c("Test Set"), 
       fill = c("red"), cex = .7, bty = "n")
dev.off()
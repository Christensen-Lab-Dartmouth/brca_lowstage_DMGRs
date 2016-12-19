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
# this will load custom CGHeatmap and heatmap.3 function
source("III.DMGR_analysis/Scripts/Functions/make_heatmaps.R") 
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
adjustedQ <- read.table("V.Validation/Data/Validation_qvalues_adjusted.csv", 
                        sep = ",", stringsAsFactors = F)
Val.Deltas <- read.table("V.Validation/Data/Validation_delta.csv", sep = ",", 
                         stringsAsFactors = F)

# Combine them into one variable
adjustedQ <- cbind(adjustedQ, Val.Deltas)

# Load extended annotation file. This file does not ignore multiple genes for a single cpg
annotation <- read.table("I.Data_Processing/Files/Expanded_annotationfile.csv", 
                         stringsAsFactors = F, row.names = 1, 
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
rownames(test_set) <- test_set[[1]]
test_set[[1]] <- NULL
test_set <- as.data.frame(test_set)

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
# Load Raw Beta Values
################################
# Load Betas
Betas <- read_tsv("V.Validation/Data/Validation_Betas.tsv")
rownames(Betas) <- Betas[[1]]
Betas[[1]] <- NULL
Betas <- as.data.frame(Betas)

# Load Annotation File with Gene Regions
annotationGR <- read.table("I.Data_Processing/Files/Expanded_annotationfile.csv", stringsAsFactors = F, row.names = 1, 
                           header = T, sep = ",", nrows = 1000000, comment.char = "")

# Load Annotation File with Map Info
annotationMap <- read_csv("I.Data_Processing/Files/HumanMethylation450K_Annotation_File.csv", skip = 7)
annotationMap <-  as.data.frame(annotationMap)

# Load covariate file
covariates <- read.table("V.Validation/Tables/ValidationSampleID.txt", sep = ",", row.names = 1, 
                         header = T, stringsAsFactors = F)
covariates <- covariates[colnames(Betas), ]

################################
# Build twelve DMGR info in validation set
################################
# We only need tumor vs. control info
cov <- data.frame(covariates$Sample_Group)
rownames(cov) <- rownames(covariates)
colnames(cov) <- "sample.type"
cancer <- rownames(cov)[cov$sample.type == "IDC"]
control <- rownames(cov)[cov$sample.type == "control"]

# Get tumor and control beta values
cancer <- Betas[ ,cancer]
control <- Betas[ ,control]

# Build the table
ValidationTable <- c()
for (DMGR in rownames(test_set)) {
  Iso_Results <- CombineOrdered[CombineOrdered$Gene_Region %in% DMGR, ]
  minQ <- min(Iso_Results$qvalues)
  minCG <- Iso_Results[Iso_Results$qvalues == minQ, 1]
  TotalCGs <- nrow(Iso_Results)
  # Subset to only significant CpGs (< 0.05 Q Value)
  Iso_Results <- Iso_Results[Iso_Results$qvalues < 0.05, ]
  if (nrow(Iso_Results) < 1) {
    rowbuild <- c('NA', minQ, 'NA', 'NA', TotalCGs, minCG)
  } else {
    # Subset to only unique CpGs
    Iso_Results <- Iso_Results[!duplicated(Iso_Results$TargetID), ]
    medB <- median(Iso_Results$beta)
    medD <- median(Iso_Results$Delta)
    sigCG <- c()
    sign <- c()
    for (cg in Iso_Results$TargetID) {
      cancerSub <- mean(as.numeric(paste(cancer[cg, ])))
      controlSub <- mean(as.numeric(paste(control[cg, ])))
      if (cancerSub > controlSub) {
        sign <- paste(sign, '+', sep = ' ')
      } else {
        sign <- paste(sign, '-', sep = ' ')
      }
      sigCG <- paste(sigCG, cg, sep = ';')
    }
    rowbuild <- c(sign, minQ, medB, medD, TotalCGs, sigCG)
  }
  ValidationTable <- rbind(ValidationTable, rowbuild)
}

rownames(ValidationTable) <- rownames(test_set)
colnames(ValidationTable) <- c('sign', 'min Q', 'med B', 'med D', 
                               'denominator', 'CpGs')

write.table(ValidationTable, 'V.Validation/Tables/ValidationDMGRs.csv',
            sep = ',', col.names = NA)

################################
# Load Common Overlaps
################################
CommonOverlaps <- read.csv("III.DMGR_analysis/Tables/commonLowStageOverlaps_FullAnnotation_extended.csv", sep = ",", 
                           row.names = 1, header = T, stringsAsFactors = F)

################################
# Get data ready for heatmap
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
source("III.DMGR_analysis/Scripts/Functions/make_heatmaps.R") 
png("V.Validation/Figures/Heatmap_Validation.png", height = 800, width = 900)
CGHeatmap("Validation Set", lowcgs, regions, annotationGR, annotationMap, 
          Covariates = covariates, Betas, Validation = T)
dev.off()

# Output Line Graph visualization
png("V.Validation/Figures/LineGraph_Validation.png", height = 700, width = 800)
par(mar = c(6, 6, 3, 9))
plot(x = 1:nrow(plotReady), y = -log10(as.numeric(paste(plotReady$qvalues))), 
     xlab = "", ylab = "", col = "black", 
     cex = 1.3, cex.lab = 2, xaxt = 'n', yaxt = 'n',
     pch = 16, bty = "n")

title(ylab = "-log10 Q value", line = 3.5, cex.lab = 3)
title(xlab = "Q Rank", line = 4, cex.lab = 3)

# Extend the length of the x axis
axis(1, at = c(0, 1e5, 2e5, 3e5, 4e5, 5e5, 6e5), col = "black", cex.lab = 2, cex.axis = 2, lwd = 2)
axis(2, col = "black", cex.lab = 2, cex.axis = 2, lwd = 2)

# Add dashed lines indicating cutoff
abline(h = -log10(0.01), col = "grey", lwd = 3, lty = "dashed")
text(x = 1e+05, y = -log10(0.01) + 0.3, "Q = 0.01", col = "grey", cex = 2.5)
abline(h = -log10(0.05), col = "grey", lwd = 3, lty = "dashed")
text(x = 1e+05, y = -log10(0.05) + 0.3, "Q = 0.05", col = "grey", cex = 2.5)

# Add points for the test set overlap
points(x = as.numeric(paste(Indices)), y = -log10(as.numeric(paste(Q))), pch = 16, 
       col = "red", cex = 1.3)

# Add text in different locations to avoid overlap
textcex = 1.7
text(x = as.numeric(paste(Indices))[1],
     y = -log10(as.numeric(paste(Q)))[1] + .15,
     rownames(test_set)[1],
     cex = textcex, pos = 2, offset = 2.1)  # HDAC4 Body

text(x = as.numeric(paste(Indices))[c(5, 8)], 
     y = -log10(as.numeric(paste(Q)))[c(5, 8)], 
     rownames(test_set)[c(5, 8)], 
     cex = textcex, pos = 2, offset = 2.1)
text(x = 600000, y = 0.2, 
     rownames(test_set)[12], cex = textcex, pos = 2, offset = 1.25) # FCGRT TSS1500
text(x = as.numeric(paste(Indices))[2] - 10000, y = -log10(as.numeric(paste(Q)))[2] + 0.05, 
     rownames(test_set)[2], cex = textcex, pos = 2, offset = 1.25)  # GNAI2 Body
text(x = as.numeric(paste(Indices))[11], y = -log10(as.numeric(paste(Q)))[11] + 0.06, 
     rownames(test_set)[11], cex = textcex, pos = 2, offset = 1.25)  # RPTOR Body
text(x = as.numeric(paste(Indices))[7], y = -log10(as.numeric(paste(Q)))[7] + 0.2, 
     rownames(test_set)[7], cex = textcex, pos = 2, offset = 1.25)  # LRP5 Body
text(x = as.numeric(paste(Indices))[10], y = -log10(as.numeric(paste(Q)))[10] + 0.2, 
     rownames(test_set)[10], cex = textcex, pos = 2, offset = 1.25)  # ANKRD11 5UTR
text(x = as.numeric(paste(Indices))[3] - 40000, y = -log10(as.numeric(paste(Q)))[3] + 7, 
     rownames(test_set)[3], cex = textcex, pos = 2, offset = .25)
text(x = as.numeric(paste(Indices))[6] - 205000, y = -log10(as.numeric(paste(Q)))[6] + 7, 
     rownames(test_set)[6], cex = textcex, pos = 2, offset = .25)
text(x = as.numeric(paste(Indices))[9] + 30000, y = -log10(as.numeric(paste(Q)))[9] + 0.17, 
     rownames(test_set)[9], cex = textcex, pos = 2, offset = 4)  # CMIP Body
text(x = as.numeric(paste(Indices))[4], y = -log10(as.numeric(paste(Q)))[4] + .025, 
     rownames(test_set)[4], cex = textcex, pos = 2, offset = 4)

legend("topleft", inset = 0.05, title = "DMGRs", c("Test Set"), 
       fill = c("red"), cex = 2, bty = "n")
dev.off()

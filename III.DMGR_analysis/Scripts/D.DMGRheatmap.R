#####################################################################
# ~~~~~~~~~~~~~~~~~~
# Tumor subtype and cell type independent DNA methylation alterations 
# associated with stage progression in invasive breast carcinoma 
# ~~~~~~~~~~~~~~~~~~
# Way, G., Johnson, K., Christensen, B. 2015
#
# Plot Heatmaps of a given list of CpGs
#####################################################################

################################
# Load Libraries
################################
#install.packages('Hmisc')
library(readr)
library(Hmisc)
# load custom CGHeatmap and heatmap.3 function
source("III.DMGR_analysis/Scripts/Functions/make_heatmaps.R") 

################################
# Load Data
################################
# Stages and subtypes of interest
low <- c("Stage I", "Stage IA", "Stage IB", "Stage II", "Stage IIA", "Stage IIB")
high <- c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC", "Stage IV")
stages <- list(low = low, high = high)
subtypes <-  c("Basal", "Her2", "LumA", "LumB", "Normal")

# Load Betas
Betas <- read_tsv("I.Data_Processing/Data/TCGA_BRCA_Betas.tsv")
rownames(Betas) <- Betas[[1]]
Betas[[1]] <- NULL
Betas <- as.data.frame(Betas)

# Load Covariates
covariates <- read.table("I.Data_Processing/Files/BRCAtarget_covariates.csv", row.names = 1, 
                         header = T, sep = ",", stringsAsFactors = F)

# Load Annotation File with Gene Regions
annotationGR <- read.table("I.Data_Processing/Files/Expanded_annotationfile.csv", 
                           stringsAsFactors = F, row.names = 1, header = T, sep = ",", 
                           nrows = 1000000, comment.char = "")

# Load Annotation File with Map Info
#annotationMap <- read_csv("I.Data_Processing/Files/HumanMethylation450_15017482_v.1.1.csv")
annotationMap <- read_csv("I.Data_Processing/Files/HumanMethylation450K_Annotation_File.csv", skip = 7)
annotationMap <- as.data.frame(annotationMap)

# The colnames for the beta file have an "X" appended to the beginning of each basename, remove it
rownames(covariates) <- covariates$Basename

################################
# Process Data
################################
# Subset the covariate data to only the samples in the beta file
covariates <- covariates[intersect(rownames(covariates), colnames(Betas)), ]
covariates <- covariates[covariates$PAM50.RNAseq != "" | covariates$sample.type == "Solid Tissue Normal", ]

# Remove Subtype Assignment of "Solid Tissue Normal" samples
covariates[covariates$sample.type == "Solid Tissue Normal", ]$PAM50.RNAseq <- ""

# Interested in "low" vs "high"
for (i in 1:length(stages)) {
  subset <- stages[[i]]
  for (j in 1:length(subset)) {
    covariates$pathologic_stage[covariates$pathologic_stage == subset[j]] <- names(stages)[i]
  }
}

# Make sure the "tumor adjacent" samples are marked in this column
covariates$pathologic_stage[covariates$sample.type == "Solid Tissue Normal"] <- "normal"

# Only accept samples that have high or low assignments
covariates <- covariates[covariates$pathologic_stage == "low" | covariates$pathologic_stage == "high" | covariates$pathologic_stage == "normal",]

# Subset Betas to those samples with PAM50 data and stage of interest
Betas <- Betas[ ,rownames(covariates)]

rm(high, i, j, low, subset, stages)

################################
# Load Common Overlaps
################################
CommonOverlaps <- read.csv("III.DMGR_analysis/Tables/commonLowStageOverlaps_FullAnnotation_extended.csv", 
                           sep = ",", row.names = 1, header = T, stringsAsFactors = F)

################################
# Prepare data for heatmaps
################################
# Extract Genes and Gene Regions
Genes <- c()
Regions <- c()
for (i in 1:nrow(CommonOverlaps)) {
  Genes <- c(Genes, unlist(strsplit(rownames(CommonOverlaps)[i], " "))[1])
  Regions <- c(Regions, unlist(strsplit(rownames(CommonOverlaps)[i], " "))[2])
}

################################
# Gene Specific Heatmap
################################
# Write all Heatmaps to file
#for (i in 12:nrow(CommonOverlaps)) {
for (i in 1:nrow(CommonOverlaps)) {
  # Function in make_heatmaps.R
  #cgs <- ExtractCommonCGs(Genes[i], CommonOverlaps)
  cgs <- ExtractCommonCGs(rownames(CommonOverlaps[i,]), CommonOverlaps)
  cgs <- unique(cgs)
  regions <- c(rep(Regions[i], length(cgs)))
  if(length(cgs) > 1){
  # Get all genespecific heatmaps and save to disk
  png(filename = paste("III.DMGR_analysis/Figures/heatmaps/", Genes[i], "_", Regions[i], 
                       "CommonHeatmap.png", sep = ""), width = 900, height = 700)
  CGHeatmap(Genes[i], cgs, regions, annotationGR, annotationMap, Betas, covariates)
  dev.off()
  }
}

################################
# Total CpG Heatmap
################################
# Create Large Heatmap with all CpGs
lowcgs <- c()
regions <- c()
for (i in 1:length(Genes)) {
  # Get all the cgs
  cg <- ExtractCommonCGs(Genes[i], CommonOverlaps = CommonOverlaps)
  lowcgs <- c(lowcgs, unique(cg))
  
  # Get all the regions
  tmpregion <- c(rep(Regions[i], length(unique(cg))))
  regions <- c(regions, tmpregion)
}

# Save to disk
png(filename = "III.DMGR_analysis/Figures/heatmaps/Heatmap_LowOverlappingCpGs_0.01.png", 
    width = 1800, height = 1400)
CGHeatmap("Common Low Stage CpGs\n (q value cutoff = 0.01)", lowcgs, regions, annotationGR, 
          annotationMap, Betas, covariates, stages = c("low", "normal"), Full = T, Validation = F)
dev.off()

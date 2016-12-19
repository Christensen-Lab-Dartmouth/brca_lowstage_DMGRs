#####################################################################
# ~~~~~~~~~~~~~~~~~~
# Tumor subtype and cell type independent DNA methylation alterations 
# associated with stage progression in invasive breast carcinoma 
# ~~~~~~~~~~~~~~~~~~
# Way, G., Johnson, K., Christensen, B. 2015
#
# Examine whether DNA methylation of selected CpGs is associated with 
# gene expression
#####################################################################

################################
#Load Libraries and the Plot Function
################################
library(readr)
library(plyr)
source("IV.Genomic_analysis/Scripts/Functions/MethRNAseq_functions.R") 
source("III.DMGR_analysis/Scripts/Functions/make_heatmaps.R") #Use the scripts in here to extract the common CpGs

################################
# Load and Subset Data
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

# Load TCGA BRCA Normal RNAseq Data
NormalRNAseq <- read_tsv("IV.Genomic_analysis/Data/unc.edu_BRCA_IlluminaHiSeq_RNASeqV2.geneExp.whitelist_normal")

rownames(NormalRNAseq) <- NormalRNAseq[[1]]
NormalRNAseq[[1]] <- NULL
NormalRNAseq <- as.data.frame(NormalRNAseq)

NormalRNAseq <- NormalRNAseq[-grep("[?]", laply(rownames(NormalRNAseq), 
                                                function (x) {unlist(strsplit(x, "[|]"))[1]})), ]


NormalRNAseq <- NormalRNAseq[-grep("SLC35E2", laply(rownames(NormalRNAseq), 
                                                    function (x) {unlist(strsplit(x, "[|]"))[1]})), ]
colnames(NormalRNAseq) <- substr(colnames(NormalRNAseq), 1, 15)
colnames(NormalRNAseq) <- gsub("-", ".", colnames(NormalRNAseq))
rownames(NormalRNAseq) <- laply(rownames(NormalRNAseq), function (x) {unlist(strsplit(x, "[|]"))[1]})

# Load annotation file
annotation <- read_csv("I.Data_Processing/Files/HumanMethylation450K_Annotation_File.csv", skip = 7)
annotation <- as.data.frame(annotation)

# Load Covariates
covariates <- read.table("I.Data_Processing/Files/BRCAtarget_covariates.csv", row.names = 1, 
                         header = T, sep = ",", stringsAsFactors = F)

# The colnames for the beta file have an "X" appended to the beginning of each basename, remove it
rownames(covariates) <- covariates$Basename

# Subset the covariate data to only the samples in the beta file and then to only the ones with PAM50 data
covariates <- covariates[intersect(rownames(covariates), colnames(Betas)), ]
covariates <- covariates[covariates$PAM50.RNAseq != "", ]

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

################################
# Run Function
################################
# Load Common Overlaps
CommonOverlaps <- read.csv("III.DMGR_analysis/Tables/commonLowStageOverlaps_FullAnnotation_extended.csv", 
                           row.names = 1, header = T, stringsAsFactors = F)

# Get all the genes in common to low stage tumors
Genes <- laply(rownames(CommonOverlaps), function(x){unlist(strsplit(x, " "))[1]})

# What are the number of comparisons made here? Bonferroni adjusted p value required.
num_unique_cpgs <- 0
for (gene in 1:length(Genes)) {
  CpGs <- unique(ExtractCommonCGs(Genes[gene], CommonOverlaps))
  num_cpgs <- length(CpGs)
  num_unique_cpgs <- num_unique_cpgs + num_cpgs
}

# 101 Unique CpGs, made for 6 comparisons (5 subtypes + all)
# Bonferroni adjustment should be made for 6 * 101 = 
alpha <- 0.05 / (6 * num_unique_cpgs)

# Loop over all genes to output several plots investigating methylation influencing gene expression
significantCor <- c()
#for (gene in 17:length(Genes)) {
for (gene in 1:length(Genes)) {
  # Extract the CGs associated with a specific gene
  #CpGs <- unique(ExtractCommonCGs(Genes[gene], CommonOverlaps))
  CpGs <- unique(ExtractCommonCGs(rownames(CommonOverlaps)[gene], CommonOverlaps))
  
  for (i in 1:length(CpGs)) {
    # Create and save all of the plots for each combination of CpGs and Genes
    png(paste("IV.Genomic_analysis/Figures/GeneExprs/", Genes[gene], "_", CpGs[i], ".png", sep = ""), 
        height = 400, width = 400)
    corTable <- methSeqPlot(gene = Genes[gene], betas = Betas, cg = CpGs[i], covariates = covariates,
                            method = 'spearman', stages = "low", subtypes = subtypes, 
                            normalExprs = NormalRNAseq)
    dev.off()
    
    # Output the Correlation Analysis to File as well
    write.table(corTable, paste("IV.Genomic_analysis/Tables/GeneExprs/", Genes[gene], "_", CpGs[i], ".csv", sep = ""),
                row.names = T, col.names = NA, sep = ",")
    
    # Test if there are any siginificant findings
    if (length(corTable[corTable[ , 3] <= alpha, ]) > 0 | length(corTable[corTable[ , 4] <= alpha, ]) > 0) {
      for (sig_row in 1:nrow(corTable)) {
        sigHit <- paste("IV.Genomic_analysis/Tables/GeneExprs/", Genes[gene], "_", CpGs[i], ".csv", sep = "")
        if (corTable[sig_row, 3] <= 0.05 | corTable[sig_row, 4] <= 0.05 ) {
          hitinclude <- c(sigHit, corTable[sig_row, ], rownames(corTable)[sig_row])
          significantCor <- rbind(significantCor, hitinclude)
        }
      }
    }
  }
}

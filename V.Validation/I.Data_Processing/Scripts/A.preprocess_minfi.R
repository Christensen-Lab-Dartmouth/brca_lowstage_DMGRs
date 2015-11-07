#####################################################################
# ~~~~~~~~~~~~~~~~~~
# Tumor subtype and cell type independent DNA methylation alterations 
# associated with stage progression in invasive breast carcinoma 
# ~~~~~~~~~~~~~~~~~~
# Way, G., Johnson, K., Christensen, B. 2015
#
# The script will curate useable samples and CpG sites and preprocess with Funnorm
# Using the package minfi as introducted by Aryee et al 2014 Workflow as described in the manuscript
#####################################################################

################################
# Load Libraries
################################
library(minfi) #v1.12.0
library(IlluminaHumanMethylation450kmanifest) #v0.4.0
library(IlluminaHumanMethylation450kanno.ilmn12.hg19) #v0.2.1
source("I.Data_Processing/Scripts/Functions/read.450k.exp2.R")

# Load Command Arguments
args <- commandArgs(trailingOnly = T)
outfile <- args[1]

##########################################
# Step 1: Load IDAT files of interest
##########################################
# Location of IDAT files for study (Place the IDAT files in a designated file location)
idat <- args[2]

# Be sure to have a .csv file of covariates in the same location of the idat files (See README)
targets <- read.450k.sheet(idat)

# Append a basename of the idat file location to the target datasheet
# The manifest file is downloaded with the TCGA data from cBioPortal
manifest <- read.csv("I.Data_Processing/Files/file_manifest.csv", skip = 2)

# Subset file to only level 1 data
manifest <- manifest[manifest[ ,"n.a"] == 1, ]

# Match the basenames to the sample names
targets$Basename <- manifest$jhu.usc.edu_BRCA.HumanMethylation450.1.9.0.sdrf.txt[match(targets$X_SAMPLE_ID, manifest$selected_samples)]  
targets$Basename <- substr(targets$Basename, 1, nchar(as.character(paste(targets$Basename))) - 9)

# Only consider the samples with idat files
targets <- targets[!is.na(targets$Basename), ]

# Write this file to the directory to use in subset script
write.table(targets, file = "I.Data_Processing/Files/BRCAtarget_covariates.csv", 
            row.names = T, col.names = NA, sep = ",")

# Reads the 'targets' like a dataframe
# The function is a custom function replacing the read.450k.exp() traditional function	
RGsetEx <- read.450k.exp2(base = idat, targets = targets)

# Access the phenotype data associated with an experiment  
pheno <- pData(RGsetEx)	

# Perform a QC report if the user supplies a "summary" argument
if (args[3] == "summary") {
  # Produce a PDF QC report                                      
  qcReport(RGsetEx, sampNames = pheno$X_SAMPLE_ID, 
           sampGroups = pheno$sample.type, pdf = "I.Data_Processing/Files/qcReport.pdf") 
  # This will append a plot to the qc report file
  mdsPlot(RGsetEx, sampNames = pheno$X_SAMPLE_ID, sampGroups = pheno$sample.type)
  # Produce a density plot of beta values       
  densityPlot(RGsetEx, sampGroups = pheno$sample.type, main = "Beta", xlab = "Beta")
}

# Filtering Probes:
# Probes with a detection p-value above signifigance (0.00001) should approached with caution
detect.p = detectionP(RGsetEx, type = "m+u")
# Begin Filtering Steps at a Detection P value of 0.00001
failed <- detect.p > 0.00001
# The detection p value must be above 0.00001 for 25% of all samples to remove the probe
cg.f <- rowMeans(failed) > 0.25
# Subset the cgs used 
cgUse <- rownames(detect.p)[!cg.f] # Remove 2,932 cgs
cat("Remove", (nrow(detect.p) - length(cgUse)), "due to high detection p values\n")

# Remove Chen and Sex Specific Probes
# Load Annotation File
annotation <- read.table("I.Data_Processing/Files/annotationfile.tsv", stringsAsFactors = F, 
                         row.names = 2, header = T, sep = "\t", nrows = 500000, comment.char = "")[ ,-1]

# Retrieve all non-Chen and non-sex specific probes and SNP probes to subset beta file
use <- rownames(annotation[annotation$excludeChen == 0 & annotation$SNPinProbe == 0 & annotation$Sex == 0, ])
cat("Remove", (length(cgUse) - length(use)), "due to Chen Probes, SNP Probes, and Sex Probes\n")

# Consider the intersection of the cgs given in the beta file (detection P value filtered) 
# and the excluded cgs
filtered.p  <- intersect(cgUse, use)

rm(manifest, idat, targets, detect.p)
#################################################                
# STEP 2: Convert intensities to methylation signal and Normalize              
################################################# 
# Functional Normalization is a between-array normalization. Removes unwanted variation by regressing
# out variability explained by the control probes present on the array. 
Mset_FunNorm <- preprocessFunnorm(RGsetEx, nPCs = 2, sex = pheno$gender)

# Use the intersect function to eliminate bad probes from FunNorm data.
# This is done using the p-value matrix from the previous steps.
intersection <- intersect(rownames(Mset_FunNorm), filtered.p)
# How many probes remain?  
length(intersection)  # 388,772

# subset the FunNorm results
Mset_FunNorm_sub = Mset_FunNorm[intersection, ]

rm(RGsetEx, pheno, Mset_FunNorm)
#################################################                
# STEP 3: Extract Beta-values              
################################################# 
# Get Beta-Values
Betas <- getBeta(Mset_FunNorm_sub)

# Write table
beta.file <- paste(outfile, "_Betas.tsv", sep = "")
write.table(Betas, file = beta.file, sep = "\t", row.names = T, col.names = NA)
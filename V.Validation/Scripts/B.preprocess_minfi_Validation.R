#####################################################################
# ~~~~~~~~~~~~~~~~~~
# Tumor subtype and cell type independent DNA methylation alterations 
# associated with stage progression in invasive breast carcinoma 
# ~~~~~~~~~~~~~~~~~~
# Way, G., Johnson, K., Christensen, B. 2015
#
# The script will curate useable samples from GSE60185 and CpG sites and preprocess with Funnorm
# Fleischer et al. 2015 http://www.ncbi.nlm.nih.gov/pubmed/25146004 Validation Set
# Using the package minfi as introducted by Aryee et al 2014 Workflow as described in the man
#####################################################################

################################
# Load Libraries
################################
#source("https://bioconductor.org/biocLite.R")
#biocLite("minfi")
#biocLite("IlluminaHumanMethylation450kmanifest")
#biocLite("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library(minfi)  # v1.12.0
library(plyr)
library(IlluminaHumanMethylation450kmanifest) #v0.4.0
library(IlluminaHumanMethylation450kanno.ilmn12.hg19) #v0.2.1
source("I.Data_Processing/Scripts/Functions/read.450k.exp2.R")

#args <- commandArgs(trailingOnly = T)
##########################################
# Step 1: Load IDAT files of interest
##########################################
# Location of IDAT files for study (Place the IDAT files in a designated file location)
# idat <- "../../Documents/IDAT"
#idat <- args[1]
idat <- "../../Documents/IDAT"

# Read a .csv file of covariates in the same location of the idat files (Be sure the other .csv file was previously removed from the IDAT file location)
manifest <- read.metharray.sheet(idat)

# Subset the Manifest File to only Low Stage Breast Tumors and Normal Samples
manifest <- manifest[manifest[,"Sample_Group"] == "IDC" | manifest[,"Sample_Group"] == "control", ]

# Confirm that the IDAT Files match the Basename in the manifest file
idat.files <- list.files(idat)
base <- unique(laply(idat.files, function (x) {paste(unlist(strsplit(x, "_"))[2], 
                                                   unlist(strsplit(x, "_"))[3], sep = "_")}))
manifest.base <- unique(laply(manifest$Basename, function (x) {paste(unlist(strsplit(x, "_"))[2], 
                                                                   unlist(strsplit(x, "_"))[3], 
                                                                   sep = "_")}))

length(intersect(manifest.base, base)) #232
# Therefore, all samples in manifest file have IDAT files

# Reads the 'targets' like a dataframe; the function is a custom function replacing the read.450k.exp() traditional function	
RGsetEx <- read.450k.exp2(base = idat, targets = manifest)

# Access the phenotype data associated with an experiment  
pheno <- pData(RGsetEx)	
write.table(pheno, file = "V.Validation/Tables/ValidationSampleID.txt", sep = ",", 
            row.names = T, col.names = NA)

# Filtering Probes:
# Probes with a detection p-value above signifigance (0.00001) should approached with caution
detect.p = detectionP(RGsetEx, type = "m+u")
# Begin Filtering Steps at a Detection P value of 0.00001
failed <- detect.p > 0.00001
# The detection p value must be above 0.00001 for 25% of all samples to remove the probe
cg.f <- rowMeans(failed) > 0.25
# Subset the cgs used 
cgUse <- rownames(detect.p)[!cg.f]  # Remove 703 cgs
cat("Remove", (nrow(detect.p) - length(cgUse)), "due to high detection p values\n")

# Remove Chen Probes and Sex Specific Probes
# Load Annotation File
annotation <- read.table("I.Data_Processing/Files/annotationfile.tsv", stringsAsFactors = F, 
                         row.names = 2, header = T, sep = "\t", nrows = 500000, comment.char = "")[ ,-1]

# Retrieve all non-Chen and non-sex specific probes and SNP probes to subset beta file
use <- rownames(annotation[annotation$excludeChen == 0 & annotation$SNPinProbe == 0 & annotation$Sex == 0, ])

cat("Remove", (length(cgUse) - length(use)), "due to Chen Probes, SNP Probes, and Sex Probes\n") #94485
# Consider the intersection of the cgs given in the beta file (detection P value filtered) and the excluded cgs
filtered.p  <- intersect(cgUse, use)

rm(manifest, idat, targets, detect.p)
#################################################                
# STEP 2: Convert intensities to methylation signal and Normalize              
################################################# 
# Functional Normalization is a between-array normalization. Removes unwanted variation by regressing
# out variability explained by the control probes present on the array. 
# Note that we do not have sex information for each sample
Mset_FunNorm <- preprocessFunnorm(RGsetEx, nPCs = 2, sex = NULL)

# Use the intersect function to eliminate bad probes from FunNorm data.
# This is done using the p-value matrix from the previous steps.

intersection <- intersect(rownames(Mset_FunNorm), filtered.p)

# How many probes remain?  
length(intersection) #390253

# subset the FunNorm results
Mset_FunNorm_sub = Mset_FunNorm[intersection, ]

rm(RGsetEx, pheno, Mset_FunNorm)
#################################################                
# STEP 3: Extract Beta-values              
################################################# 
# Get Beta-Values
Betas <- getBeta(Mset_FunNorm_sub)

# Write table
write.table(Betas, file = "V.Validation/Data/Validation_Betas.tsv", sep = "\t", row.names = T, 
            col.names = NA)

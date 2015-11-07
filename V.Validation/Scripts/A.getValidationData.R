#####################################################################
# ~~~~~~~~~~~~~~~~~~
# Tumor subtype and cell type independent DNA methylation alterations 
# associated with stage progression in invasive breast carcinoma 
# ~~~~~~~~~~~~~~~~~~
# Way, G., Johnson, K., Christensen, B. 2015
#
# Download Fleischer et al 2015 Validation Set
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60185
#####################################################################

################################
# Load commandArgs
################################
args <- commandArgs(trailingOnly = T)
destination <- args[1]

# Create the directory if it does not already exist (will not rewrite)
dir.create(destination, recursive = T, showWarnings = F)

################################
# Load Libraries
################################
library(GEOquery)

################################
# Get Data
################################
# load series and platform data from GEO, the IDAT files are located in the supplement
getGEOSuppFiles("GSE60185", baseDir = "V.Validation/Data/") 

# Get Phenotype Data
gse <- getGEO("GSE60185")
pheno <- pData(gse[[1]])
write.table(pheno, "V.Validation/Tables/Validation_PhenoData_GSE58999.csv", sep = ",", 
            row.names = T, col.names = NA)

rm(gse, pheno)
################################
# Unzip the GEO data
################################
# First, untar the folder download into an external directory
untar("V.Validation/Data/GSE60185/GSE60185_RAW.tar", exdir = destination)

# Next, gunzip the files in this newly untarred directory replacing all files
files <- list.files("~/Documents/IDAT/")[grep(".gz", list.files("~/Documents/IDAT/"))]
for (i in 1:length(files)) {
  gunzip(paste(destination, files[i], sep = ""), overwrite = T)
}

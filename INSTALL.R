#####################################################################
# ~~~~~~~~~~~~~~~~~~
# Tumor subtype and cell type independent DNA methylation alterations 
# associated with stage progression in invasive breast carcinoma 
# ~~~~~~~~~~~~~~~~~~
#
# Way, G., Johnson, K., Christensen, B. 2015
#
# Run This script to install all the R libraries required to run the
# analysis pipeline.
#####################################################################

#####################
# CRAN libraries
#####################
cran <- c("plyr", "reshape2", "ggplot2", "qvalue", "gridExtra", "readr", "Hmisc",
          "VennDiagram")
install.packages(cran, repos='http://cran.us.r-project.org')

#####################
# Bioconductor libraries
#####################
#First, install bioconductor
source("http://bioconductor.org/biocLite.R")
biocLite()

bioC <- c("minfi", "IlluminaHumanMethylation450kmanifest", "IlluminaHumanMethylation450kanno.ilmn12.hg19",
          "RefFreeEWAS", "isva", "limma", "Homo.sapiens", "Gviz", "GEOquery")

biocLite(bioC)


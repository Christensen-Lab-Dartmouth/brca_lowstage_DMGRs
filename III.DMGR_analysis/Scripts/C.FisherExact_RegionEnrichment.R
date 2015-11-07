#####################################################################
# ~~~~~~~~~~~~~~~~~~
# Tumor subtype and cell type independent DNA methylation alterations 
# associated with stage progression in invasive breast carcinoma 
# ~~~~~~~~~~~~~~~~~~
# Way, G., Johnson, K., Christensen, B. 2015
#
# Investigate enrichment of specific regions observed to be significantly 
# differentialy methylated
#####################################################################

################################
# Load Libraries
################################
library(readr)
library(plyr)
source("III.DMGR_analysis/Scripts/Functions/DMcgs_functions.R")

################################
# Load Constants
################################
UseRegions <- c("Body", "1stExon", "3'UTR", "5'UTR", "TSS1500", "TSS200")

################################
# Load Data
################################
# Load extended annotation file. This file does not ignore multiple genes for a single cpg
annotation <- read.table("I.Data_Processing/Files/Expanded_annotationfile.csv", 
                         stringsAsFactors = F, row.names = 1, header = T, sep = ",", 
                         nrows = 1000000, comment.char = "")

# Append a column of only the first gene and first region to the annotation file
annoFirst <- apply(annotation, 1, function (x) {splitGene(x[12])[1]})
annoReg <- apply(annotation, 1, function (x) {splitGene(x[14])[1]})
annotation <- cbind(annotation, paste(annoFirst, annoReg))
colnames(annotation)[ncol(annotation)] <- "GeneRegion"

# Load low stage overlaps
common <- read_csv("III.DMGR_analysis/Tables/commonLowStageOverlaps_FullAnnotation_extended.csv")

################################
# Summarize Region Enrichment 
################################
# Observe Total Counts
TotalReg <- table(annotation$UCSC_RefGene_Group)[UseRegions]

# Get differentially methylated regions in common
ObsRegions <- table(laply(common[ ,1], function (x) {unlist(strsplit(x, " "))[2]}))

# Get the data in the proper format
Other <- c(nrow(annotation[annotation$UCSC_RefGene_Group != "", ]) - as.numeric(TotalReg["Body"]), 
           as.numeric(TotalReg["Body"]))
Body <- c(nrow(common) - as.numeric(ObsRegions["Body"]), as.numeric(ObsRegions["Body"]))

# Obtain Contingency Table and perform Fisher's Exact Test
Contingency <- cbind(Other, Body)
fisher <- fisher.test(Contingency)

# extract info
CI <- paste(round(fisher$conf.int[1], 2), "-", round(fisher$conf.int[2], 2))
Odds <- round(fisher$estimate,2)
Pval <- format(fisher$p.value,scientific = T)
Results <- c(Odds, CI, Pval)
names(Results) <- c("OR", "95% CI", "P")

# Write to file
write.table(Results, file = "III.DMGR_analysis/Tables/GeneBodyEnrichmentTest.csv", 
            sep = ",", row.names = T, col.names = NA)

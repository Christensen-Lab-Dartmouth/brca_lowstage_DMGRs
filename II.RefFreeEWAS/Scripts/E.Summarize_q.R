#####################################################################
# ~~~~~~~~~~~~~~~~~~
# Tumor subtype and cell type independent DNA methylation alterations 
# associated with stage progression in invasive breast carcinoma 
# ~~~~~~~~~~~~~~~~~~
# Way, G., Johnson, K., Christensen, B. 2015
#
# Summarize q values resulting from RefFree
#####################################################################

################################
# Load Libraries
################################
library(plyr)

################################
# Constants
################################
subtypes <- c("Basal", "Her2", "LumA", "LumB", "Normal")
stages <- c("low", "high")

################################
# Load Functions
################################
# Script with custom functions
source("II.RefFreeEWAS/Scripts/Functions/doRefFree_functions.R")
source("III.DMGR_analysis/Scripts/Functions/DMcgs_functions.R")

################################
# Load and Subset Data
################################
# Load model specific q values associated with single CpGs
f.location <- "II.RefFreeEWAS/Data/"
files <- list.files(f.location)

# Obtain all the adjusted file locations
qfiles <- c()
for (i in 1:length(subtypes)) {
  for (j in 1:length(stages)) {
    qfile <- paste(f.location, files[grepl(paste(subtypes[i], "_", stages[j], "_", 
                                                 "qvalues_adjusted", sep = ""), files)], sep = "")
    qfiles <- c(qfiles, qfile)
  }
}

# Read in all the adjusted q values and store in a master list
QFileList <- list()
for (i in 1:length(qfiles)) {
  QFileList[[i]] <- read.table(qfiles[i], stringsAsFactors = F, row.names = 1, header = T, sep = ",")
  # The identifier is built into the filename; extract it
  name <- unlist(strsplit(unlist(strsplit(unlist(strsplit(qfiles[i], "[.]"))[2], "/"))[3], "_"))[1:2]
  names(QFileList)[i] <- paste(name[1], name[2], sep = "_")
}

# Load Extended Annotation File
annotation <- read.table("I.Data_Processing/Files/Expanded_annotationfile.csv", stringsAsFactors = F, 
                         row.names = 1, header = T, sep = ",", nrows = 1000000, comment.char = "")

# Append a column of only the first gene and first region to the annotation file
annoFirst <- apply(annotation, 1, function(x){splitGene(x[12])[1]})
annoReg <- apply(annotation, 1, function(x){splitGene(x[14])[1]})
annotation <- cbind(annotation, paste(annoFirst, annoReg))
colnames(annotation)[ncol(annotation)] <- "GeneRegion"

################################
# Combine annotation file and CpG specific q values for all models
################################
use <- rownames(QFileList[[1]])
AnnoQList <- list()

for (i in 1:length(QFileList)) {
  # get the QFile ready to match to the annotation file
  tmp <- QFileList[[i]][match(annotation$TargetID, use), ]
  # combined file 
  AnnoQList[[i]] <- cbind(annotation, tmp)
  names(AnnoQList)[i] <- names(QFileList)[i]
}






rm(annotation, tmp, annoFirst, annoReg, f.location, files, i, j, name, qfile, QFileList, qfiles, use)
################################
# Plot Various Q Value Cutoffs
################################
qrange <- seq(0, 0.1, by =.01)

for (i in 1:length(stages)) {
  lengthCommon <- c()
  for (j in 1:length(qrange)) {
    goodq_val <- findGoodQ(qcut = qrange[j], AnnoQList, subtypes = subtypes[1:4], stages[i], returning = "Values")
    lengthCommon <- c(lengthCommon, length(goodq_val))
  }
  
  # Get title of plot ready
  if (stages[i] == 'low') {
    stage_name = 'Early'
  } else {
    stage_name = 'Late'
  }
  
  # Make plot and save figures
  png(paste("II.RefFreeEWAS/Figures/QValueThreshold_", stage_name, ".png", sep = ""), height = 400, width = 550) 
  par(mar = c(4,6,4,4))
  plot(x = qrange, y = lengthCommon, pch = 19, xlab = "Q Value Cutoffs", xlim = c(-.01, .11) ,
       ylab = "Number of Overlapping\nGene Regions", 
       main = paste("Q value Threshold\nOverlapping", stage_name, "Stage"))
  abline(h = lengthCommon[2], lty = 2, col = "red")  # lengthCommon[2] = 0.01 q cut
  abline(h = lengthCommon[6], lty = 2, col = "blue")  # lengthCommon[6] = 0.05 q cut
  text(qrange, lengthCommon, paste(qrange, "-", lengthCommon), cex = 0.75, pos = 4)
  dev.off()
}

################################
# Extract Top 400 Overlapping Genes in Low Stage Overlaps
################################
top <- findGoodQ(qcut = 0.1, AnnoQList, subtypes = subtypes[1:4], "low", returning = "All")
genes <- lapply(rownames(top), function(x){unlist(strsplit(x, " "))[1]})
restrictTop <- unique(genes)[1:400]
write.table(restrictTop, "II.RefFreeEWAS/Tables/Top400SimilarGenes.txt", sep = ",", row.names = F)


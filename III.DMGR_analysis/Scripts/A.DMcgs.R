#####################################################################
# ~~~~~~~~~~~~~~~~~~
# Tumor subtype and cell type independent DNA methylation alterations 
# associated with stage progression in invasive breast carcinoma 
# ~~~~~~~~~~~~~~~~~~
# Way, G., Johnson, K., Christensen, B. 2015
#
# The script will subset the identified differentially methylated gene 
# regions (DMGRs) and describe direction of aberrent methylation for a 
# given q value cutoff.
#####################################################################

################################
# Load Libraries
################################
library(minfi)
library(readr)

# Script with custom functions
source("II.RefFreeEWAS/Scripts/Functions/doRefFree_functions.R")
source("III.DMGR_analysis/Scripts/Functions/DMcgs_functions.R")

################################
# Load Data
################################
# Load total covariate file
covariates_full <- read.table("I.Data_Processing/Files/BRCAtarget_covariates.csv", row.names = 1,
                         header = T, sep = ",", stringsAsFactors = F)
rownames(covariates_full) <- covariates_full$Basename

# Load extended annotation file. This file does not ignore multiple genes for a single cpg
annotation <- read.table("I.Data_Processing/Files/Expanded_annotationfile.csv", stringsAsFactors = F, 
                         row.names = 1, header = T, sep = ",", nrows = 1000000, comment.char = "")

# Append a column of only the first gene and first region to the annotation file
annoFirst <- apply(annotation, 1, function(x){splitGene(x[12])[1]})
annoReg <- apply(annotation, 1, function(x){splitGene(x[14])[1]})
annotation <- cbind(annotation, paste(annoFirst, annoReg))
colnames(annotation)[ncol(annotation)] <- "GeneRegion"

################################
# Analysis
################################
#s <- c('low', 'high')
#t <- c('Basal', 'Her2', 'LumA', 'LumB', 'Normal')

s <- c('low')
#s <- c('high')

t <- c('Basal')
#t <- c('Her2')
#t <- c('LumA')
#t <- c('LumB')
#t <- c('Normal')

# for debugging only 
#model_stage <- s
#model_subtype <- t

for (model_stage in s) {
  for (model_subtype in t) {
    args <- c(model_stage, model_subtype, "0.01")
    print(args)
    stage <- args[1]
    subtype <- args[2]
    qvalcut <- as.numeric(paste(args[3]))
    
    ################################
    # Load and Subset Data
    ################################
    f.location <- "I.Data_Processing/Data/"
    files <- list.files(f.location)
    beta.file <- paste(f.location, files[grep(subtype, files)], sep = "")
    cat(beta.file, "\n")
    
    # Load the file holding beta values corresponding to the specific subtype
    beta2 <- read_tsv(beta.file, col_names = T)
    rownames(beta2) <- beta2[['X1']]
    beta2[ ,1] <- NULL
    beta2 <- as.data.frame(beta2)
    
    # Subset the covariate data to only the samples in the beta file
    covariates <- covariates_full[intersect(rownames(covariates_full), colnames(beta2)), ]
    
    # Label samples
    cancer <- intersect(colnames(beta2), rownames(covariates[grepl("Tumor", covariates$sample.type),]))
    normal <- intersect(colnames(beta2), rownames(covariates[grepl("Normal", covariates$sample.type),]))
    
    ################################
    # Obtain q values with cutoffs
    ################################
    # args may state to get differentially methylated cpgs on the entire list of samples, or, 
    # only for a specific stage. The output of the qvalList() function are two data frames 
    # in a list that indicate CpG q values for adjusted and unadjusted models
    q <- qvalList(beta = beta2, subtype = subtype, stage = stage, covariates = covariates, 
                  qvalcut = qvalcut, filelist = "II.RefFreeEWAS/Data/")
    
    # This loop will combine q values with annotation cgs and gene regions
    anno.sub <- list()
    for (i in 1:length(q[[1]])) {

      if (nrow(q[[1]][[i]]) != 0) {
        # Extract the q value information
        frame <- q[[1]][[i]]

        # name it
        frame <- cbind(rownames(frame), frame)

        # subset the annotation file to only those with significant q values
        sub <- annotation[annotation$TargetID %in% rownames(frame), ]

        # only consider cpgs that have gene information
        sub <- sub[sub[ , 12] != "", ]

        # subset q value info
        frame.sub <- frame[match(sub$TargetID, rownames(frame)), ]

        # combine q value info with annotation
        anno.sub[[i]] <- cbind(sub, frame.sub)
        names(anno.sub)[i] <- names(q[[1]])[i]

      } else {
        cat("No Significant q value cgs")
      }

    }
    
    ################################
    # Compile DMGRs, location, direction and write to file
    ################################
    # Apply the functions in "III.DMGR_analysis/Scripts/Functions/DMcgs_functions.R" 
    # to output specific matrix and write to file
    for (j in 1:length(anno.sub)) {
      # for debugging only
      #anno.sub[[j]] <- anno.sub[[j]][1:100, ]
                  
      # Subset Beta File
      betaSub <- beta2[unique(anno.sub[[j]]$TargetID), rownames(q[[2]])]
      
      # Initialize a new dataframe to store gene info
      newFrame <- c()
      for (i in 1:nrow(anno.sub[[j]])) {
        
        row <- anno.sub[[j]][i, ]
        
        # use function getGeneInfo to extract information for each row
        info <- getGeneInfo(row)
        info <- c(paste(info[1], info[2]), paste(info[1], info[3]), info[4:length(info)])
        newFrame <- rbind(newFrame, info)
      }
      
      # Extract the significantly differentially methylated cpgs
      cpgs <- newFrame[,6]
      cat("Number of total:", length(cpgs), "\n", "Number of Unique:", length(unique(cpgs)), "\n")
      
      # Subset to cancer and normal samples
      Csamps <- intersect(colnames(betaSub), cancer)
      Nsamps <- intersect(colnames(betaSub), normal)
      
      # Observe the direction of the differential methylation
      methdir <- c()
      for (c in 1:length(cpgs)) {
        cancerval <- mean(as.numeric(paste(betaSub[cpgs[c], Csamps])), na.rm = T)
        normalval <- mean(as.numeric(paste(betaSub[cpgs[c], Nsamps])), na.rm = T)
        
        if (cancerval > normalval) {
          dir <- "+"
        } else if(cancerval < normalval) {
          dir <- "-"
        } else {
          dir <- "."
        }
        methdir <- c(methdir, dir)
      }
      
      # Combine it together
      newFrame <- cbind(newFrame, methdir)
      
      # Collapse info into writable dataframe
      methylatedFrame <- collapseInfo(newFrame, 1, annotation)
      write.table(methylatedFrame, file = paste("III.DMGR_analysis/Results/DMR_", names(anno.sub)[j], 
                                                qvalcut, "_with_UCSC_Region.csv", sep = ""), row.names = F, sep = ",")
    }
  }
}

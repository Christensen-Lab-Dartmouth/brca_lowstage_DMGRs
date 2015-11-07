#####################################################################
# ~~~~~~~~~~~~~~~~~~
# Tumor subtype and cell type independent DNA methylation alterations 
# associated with stage progression in invasive breast carcinoma 
# ~~~~~~~~~~~~~~~~~~
# Way, G., Johnson, K., Christensen, B. 2015

# The script will expand the annotation file into considering multiple 
# genes, regions, and transcripts for all single CpG probes
#####################################################################
################################
# Load Libraries
################################
library(plyr)

################################
# Load Data
################################
# Load and order the annotation file
# This annotation file was curated by Andy Houseman
annotation <- read.table("I.Data_Processing/Files/annotationfile.tsv", stringsAsFactors = F, 
                         row.names = 1, header = T, sep = "\t", nrows = 500000, comment.char = "")

################################
# Load Function
################################
# Function will split the cg if it is explanatory of multiple genes and regions
# The function is written to be used with apply(), so, the function is applied to each row
expandAnnot <- function(row) {
  # Make sure the variable is of character type
  row <- as.character(row)
  
  # If once you split the gene name and there is more than 1 gene, perform the following:
  # Note- The annotation file has the UCSC_refGene in column 12
  if (length(unlist(strsplit(row[12], ";"))) > 1) {
    # Split the gene name
    name <- unlist(strsplit(row[12], ";"))
    # Split the gene transcript
    reg <- unlist(strsplit(row[13], ";"))
    # Split the gene region
    group <- unlist(strsplit(row[14], ";"))
    # Make sure that the splits are all the same
    if (length(name) == length(reg) & length(name) == length(group) & length(reg) == length(group)) {
      # Build a character matrix the same dimension as the amount of genes
      build <- c()
      for (i in 1:length(name)) {
        build <- rbind(build, row)
      }
      # Place the results of the split into the appropriate positions in the build matrix
      build[ ,12] <- name
      build[ ,13] <- reg
      build[ ,14] <- group
      # The columns are named the same as the annotation file
      colnames(build) <- colnames(annotation)
      # Return this information if a single cpg maps to multiple elements
      return(build)
      
      # If the gene name, regions, and transcripts are not all the same length, 
      # then return the same row in matrix form without the split
    } else {
      tmp <- matrix(NA, 1, length(row))
      tmp[1, ] <- row
      colnames(tmp) <- colnames(annotation)
      return(tmp)
    }
    
    # If the CpG only describes one element, return it in matrix form
  } else {
    tmp <- matrix(NA, 1, length(row))
    tmp[1, ] <- row
    colnames(tmp) <- colnames(annotation)
    return(tmp)
  }
}

################################
# Expand Annotation File
################################
# Apply the expandAnnot() function to the rows of the annotation file; this results in a list
annotation_Full <- apply(annotation, 1, expandAnnot)

# rbind all elements of the list to make the newAnnotation file
newAnnotation <- do.call(rbind, annotation_Full)

# Assign sequential rownames to the new annotation file
rownames(newAnnotation) <- seq(1, nrow(newAnnotation), by = 1)

# Make it a data.frame to get ready for writing to disk
writeReady <- data.frame(newAnnotation)

# write file
write.table(writeReady, file = "I.Data_Processing/Files/Expanded_annotationfile.csv", sep = ",", 
            row.names = T, col.names = NA)
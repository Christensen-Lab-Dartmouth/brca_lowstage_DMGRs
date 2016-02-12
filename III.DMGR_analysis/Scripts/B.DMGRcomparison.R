#####################################################################
# ~~~~~~~~~~~~~~~~~~
# Tumor subtype and cell type independent DNA methylation alterations 
# associated with stage progression in invasive breast carcinoma 
# ~~~~~~~~~~~~~~~~~~
# Way, G., Johnson, K., Christensen, B. 2015
#
# Output venn diagrams of overlapping gene-regions and tables describing 
# low and high stage overlaps
#####################################################################

################################
# Load Libraries
################################
library(limma)
library(VennDiagram)

################################
# Load Functions
################################
# This function will determine which, if any CpGs overlap in the set
countcgs <-  function(x) {
  tmp <- intersect(unlist(strsplit(x[5], ";")), 
                   intersect(unlist(strsplit(x[10], ";")), 
                             intersect(unlist(strsplit(x[15], ";")), 
                                       unlist(strsplit(x[20], ";")))))
  return(tmp)
}

# This function will count the number of unique cg hits
countUniquecgs <-  function(x) {
  tmp <- c(unlist(strsplit(x[5], ";")), unlist(strsplit(x[10], ";")), 
           unlist(strsplit(x[15], ";")), unlist(strsplit(x[20], ";")))
  return(length(unique(tmp)))
}

# Use to output subtype specific information for plotting Venn diagram
getgenevenn <- function(venngenes, subtype_column) {
  venngenes <- as.data.frame(venngenes)
  return_venn_count <- c()
  getcolumn <- venngenes[ ,subtype_column]
  for (g in 1:length(getcolumn)) {
    if (getcolumn[g] == 1) {
      return_venn_count <- c(return_venn_count, g)
    }
  }
  return(return_venn_count)
}

################################
# Load Data
################################
file.loc <- "III.DMGR_analysis/Results/"
masterFiles <- list.files(file.loc)

# Load covariate tumor purity summary file
covariatesummary <- read.table("I.Data_Processing/Tables/covariateTumorPuritySummary.csv", 
                               sep = ",", row.names = 1, header = T)

# Load Adjusted Files
AdjustedFiles <- masterFiles[grepl("_adjusted", masterFiles)]
AdjustedList <- list()
for (i in 1:length(AdjustedFiles)) {
  AdjustedList[[i]] <- read.table(paste(file.loc,AdjustedFiles[[i]], sep = ""), header = T, 
                                  sep = ",", stringsAsFactors = F)
  cat(nrow(AdjustedList[[i]]), ": ", AdjustedFiles[i], "\n")
  names(AdjustedList)[i] <- AdjustedFiles[i]
}

# Load Unadjusted Files
UnadjustedFiles <- masterFiles[grepl("_unadjusted", masterFiles)]
UnadjustedList <- list()
for (i in 1:length(UnadjustedFiles)) {
  UnadjustedList[[i]] <- read.table(paste(file.loc, UnadjustedFiles[[i]], sep = ""), 
                                    header = T, sep = ",", stringsAsFactors = F)
  cat(nrow(UnadjustedList[[i]]), ": ", UnadjustedFiles[i], "\n")
  names(UnadjustedList)[i] <- UnadjustedFiles[i]
}

#Load Full Covariate File
covariates <- read.table("I.Data_Processing/Files/BRCAtarget_covariates.csv", row.names = 1, 
                         header = T, sep = ",", stringsAsFactors = F)

################################
# Load Constants
################################
# the index where the gene regions are held in the annotation file
ind <- 1
subtypes <- c("Basal", "Her2", "LumA", "LumB", "Normal")

################################
# Within Subtype Comparison
################################
# Initialize several lists and characters
summaryModelList <- list()
TotalRegions <- c()
MasterList <- list()
for (i in 1:length(subtypes)) {
  # Get the files in the adjusted list according to the given subtype
  files <- grep(subtypes[i], names(AdjustedList))  
  
  # Loop over the number of files in the AdjustedList that correspond to the given subtype
  sheets <- list()
  for (j in 1:length(files)) {
    sheets[[j]]<- AdjustedList[[files[j]]]
    names(sheets)[j] <- names(AdjustedList)[files[j]]
    
    # rename the column with gene regions
    colnames(sheets[[j]])[ind] <- c("Gene:Region")
    cat(nrow(sheets[[j]]), "\n")
  }
  
  # loop over the adjusted files (high vs. low subtype) to make a character string holding region information
  regions <- c()
  for (k in 1:length(sheets)) {
    if (k == 1) {
      regions <- sheets[[k]][ ,ind]
    } else {
      regions <- c(regions, sheets[[k]][ ,ind])
    }  
  } 
  
  # get all the unique gene regions
  regions <- unique(regions)
  
  # build a background of Total regions
  TotalRegions <- c(TotalRegions, regions)
  cat(subtypes[i], ": Number of unique gene regions:", length(regions), "\n")
  
  # Create a matrix to visualize Gene:Region overlaps within subtype across stage
  sumModel <- matrix(NA, nrow = length(regions), ncol = length(sheets))
  for (k in 1:length(sheets)) {
    vector <- c()
    for (l in 1:length(regions)) {
      if (regions[l] %in% sheets[[k]][ ,ind]) {
        tmp <- 1
      } else {
        tmp <- 0
      }
      vector <- c(vector, tmp)
    }
    sumModel[ ,k] <- vector
  }
  rownames(sumModel) <- regions
  colnames(sumModel) <- c(paste(subtypes[i], "high", sep = "_"), paste(subtypes[i], "low", sep = "_"))
  
  # store info in the summary model list  
  summaryModelList[[i]] <- sumModel
  names(summaryModelList)[i] <- subtypes[i]
  
  # Look at overlapping Regions (high vs. low vs. overlap)
  lowindex <- c()
  highindex <- c()
  overlapindex <- c()
  for (p in 1:nrow(sumModel)) {
    row <- sumModel[p,]
    if (row[1] == 1 & row[2] == 1) {
      overlapindex <- c(overlapindex, p)  
    } else if (row[1] == 1 & row[2] == 0) {
      highindex <- c(highindex, p)
    } else {
      lowindex <- c(lowindex, p)
    }
  }
  
  both.h <- c(overlapindex, highindex)
  both.l <- c(overlapindex, lowindex)
  
  # get gene regions included in specific sets
  High <- sheets[[1]][sheets[[1]][,ind] %in% rownames(sumModel[both.h,]),]
  Low <- sheets[[2]][sheets[[2]][,ind] %in% rownames(sumModel[both.l,]),]
  HighOnly <- sheets[[1]][sheets[[1]][,ind] %in% rownames(sumModel[highindex,]),]
  LowOnly <- sheets[[2]][sheets[[2]][,ind] %in% rownames(sumModel[lowindex,]),] 
  OverlapHigh <- sheets[[1]][sheets[[1]][,ind] %in% rownames(sumModel[overlapindex,]),]
  OverlapLow <- sheets[[2]][sheets[[2]][,ind] %in% rownames(sumModel[overlapindex,]),]
  
  MasterList[[i]] <- list("High" = High, "Low" = Low, "High-Only" = HighOnly, "Low-Only" = LowOnly, 
                          "Overlap" = list("OverlapHigh" = OverlapHigh, "OverlapLow" = OverlapLow))
  names(MasterList)[i] <- subtypes[i]
}

TotalRegions <- unique(TotalRegions)
rm(AdjustedFiles, AdjustedList, files, highindex, i, j, k, l, lowindex, masterFiles, overlapindex, p, row, sheets, tmp, UnadjustedFiles, UnadjustedList, vector, High, Low, OverlapHigh, OverlapLow, regions, sumModel, HighOnly, LowOnly, both.h, both.l)
################################
# Venn Diagrams of overlap within subtype
################################
for (i in 1:length(summaryModelList)) {
  # Prepare Venn Diagram from Summary Model List
  venn <- vennCounts(summaryModelList[[i]])
  
  # Get text info (tumor purity and sample size) to paste onto venn diagram
  subs <- covariatesummary[grepl(subtypes[i], rownames(covariatesummary)),]
  
  # this will hold the text information to paste onto the plot
  texting <- c()
  for (j in 1:nrow(subs)) {
    string <- paste("n =", subs[j,1], " Tumor Mean =", subs[j,2], "\n        Tumor SD =", subs[j,3])
    texting <- c(texting, string)
  }
  
  # Write png figure to disk
  png(filename = paste("III.DMGR_analysis/Figures/", names(summaryModelList)[i], "_FullAnnotation.png", 
                        sep = ""), width = 700, height = 700)
  vennDiagram(venn, main = paste("\nOverlapping Gene Regions across Stage: \nFull Annotation -", subtypes[i], 
                                 sep = ""))
  #text(x = c(-1.8,1.8), y = c(-1.7,-1.7),  labels = texting)
  dev.off()
}

rm(i, j, string, texting, venn)
################################
# View "High" Subtype Gene Region Overlap
################################
Intersections <- c()
for (i in 1:length(MasterList)) {
  
  builder <- c()
  for (j in 1:length(MasterList)) {
    tmp <- length(intersect(MasterList[[i]][[1]][ ,ind], MasterList[[j]][[1]][ ,ind]))
    builder <- c(builder, tmp)
  }
  Intersections <- rbind(Intersections, builder)
}
colnames(Intersections) <- rownames(Intersections) <- subtypes

################################
# View "Low" Subtype Gene Region Overlap
################################
IntersectionsLow <- c()
for (i in 1:length(MasterList)) {
  
  builder <- c()
  for (j in 1:length(MasterList)) {
    tmp <- length(intersect(MasterList[[i]][[2]][ ,ind], MasterList[[j]][[2]][ ,ind]))
    builder <- c(builder, tmp)
  }
  
  IntersectionsLow <- rbind(IntersectionsLow, builder)
}
colnames(IntersectionsLow) <- rownames(IntersectionsLow) <- subtypes

################################
# Venn Diagram of Gene Region Overlaps across subtypes within stage
################################
stage <- c("high", "low")
VennLabels <- list()
for (k in 1:length(stage)) {
  # Subset covariate summary by stage to extract text for plot
  subs <- covariatesummary[grepl(stage[k], rownames(covariatesummary)),]

  # initialize a venn matrix to get ready for overlaps 
  # (-1 subtype because we are not considering Normal-Like for this analysis)
  Venn <- matrix(NA, nrow = length(TotalRegions), ncol = (length(subtypes)-1)) 
  rownames(Venn) <- TotalRegions
  colnames(Venn) <- subtypes[1:4]
  for (i in 1:(length(MasterList) - 1)) {
    vector <- c()
    for (j in 1:length(TotalRegions)) {
      if(TotalRegions[j] %in% MasterList[[i]][[k]][ ,ind]) {
        tmp <- 1
      } else {
        tmp <- 0
      }
      vector <- c(vector, tmp)
    }
    Venn[ ,i] <- vector
  }
  
  VennLabels[[k]] <- Venn
  names(VennLabels)[k] <- stage[k]
  
  basal_venn <- getgenevenn(Venn, 1)
  her2_venn <- getgenevenn(Venn, 2)
  luma_venn <- getgenevenn(Venn, 3)
  lumb_venn <- getgenevenn(Venn, 4)
  venn.plot <- venn.diagram(x = list('Basal' = basal_venn,
                                     'Lum B' = lumb_venn,
                                     'Her2' = her2_venn,
                                     'Lum A' = luma_venn),
                            filename = paste("III.DMGR_analysis/Figures/Venn_", stage[k], ".png", sep = ''),
                            height = 3000, width = 2150,
                            fill = c("red", "cyan", "pink", "blue"),
                            cat.cex = rep(1.6, 4),
                            margin = 0.02, cex = 1.2,
                            main.cex = 2)
}

################################
# Specifically Investigating Overlapping Gene Regions associated with Low Stage
################################
#These are the regions in common for low stage (across all four PAM50 subtypes)
CommonLow <- rownames(VennLabels[[2]][VennLabels[[2]][,1] == 1 & VennLabels[[2]][,2] == 1 & VennLabels[[2]][,3] == 1 & VennLabels[[2]][,4] ==  1,])

commonSummary <- c()
for (i in 1:(length(MasterList)-1)) {
  
  # subset the master list subtype specific significant DMGRs with the common low gene regions
  tmp <- MasterList[[i]][[2]][MasterList[[i]][[2]][,ind] %in% CommonLow, ]
  
  # extract information from this subset
  Q <- tmp$medQval
  B <- tmp$medBetaCoef
  D <- tmp$medDelta
  Sign <- tmp$Sign
  cgs <- tmp$cpgs
  denom <- tmp$denominator
  commonSummary <- cbind(commonSummary, Sign, Q, B, D, cgs)
}

# Name the elements of the common summary
commonSummary <- cbind(commonSummary, denom)
rownames(commonSummary) <- tmp[,1]
colnames(commonSummary) <- c("Basal_sign", "Basal_q", "Basal_B", "Basal_D", "Basal_cpgs",
                             "Her2_sign", "Her2_q", "Her2_B", "Her2_D", "Her2_cpgs", 
                             "LumA_sign", "LumA_q", "LumA_B", "LumA_D","LumA_cpgs", 
                             "LumB_sign", "LumB_q", "LumB_B", "LumB_D", "LumB_cpgs", "Denominator")

# Perform this for every row in the commonSummary matrix
int <- unlist(apply(commonSummary, 1, countcgs))
uni <- unlist(apply(commonSummary, 1, countUniquecgs))

# combine this to the commonsummary table
commonSummary <- cbind(commonSummary, int[match(rownames(commonSummary), names(int))], 
                       uni[match(rownames(commonSummary), names(uni))])
colnames(commonSummary)[22] <- "intersecting_cgs"
colnames(commonSummary)[23] <- "unique_cgs"

# Write to file
write.table(commonSummary, file = "III.DMGR_analysis/Tables/commonLowStageOverlaps_FullAnnotation_extended.csv", sep = ",", row.names = T, col.names = NA)

################################
# Common Low Stage DMGRs also associated with High Stage
################################
# Are these genes also present and where in high stage samples?
commonSummaryHigh <- c()
for (i in 1:(length(MasterList) - 1)) {
  #Subset the Master list to high stage models
  tmp1 <- MasterList[[i]][[1]][MasterList[[i]][[1]][ ,ind] %in% CommonLow,]
  
  # Extract information from the Master List
  what <- tmp1[match(CommonLow, tmp1[,ind]),]
  Q <- what$medQval
  B <- what$medBetaCoef
  D <- what$medDelta
  Sign <- what$Sign
  cgs <- what$cpgs
  demon <- what$denominator
  commonSummaryHigh <- cbind(commonSummaryHigh, Sign, Q, B, D, cgs)
}

commonSummaryHigh <- cbind(commonSummaryHigh, denom)
rownames(commonSummaryHigh) <- tmp[,1]
colnames(commonSummaryHigh) <- colnames(commonSummary)[c(-22, -23)]

write.table(commonSummaryHigh, file = "III.DMGR_analysis/Tables/commonLowStageOverlaps_HighStage.csv", sep = ",", row.names = T, col.names = NA)

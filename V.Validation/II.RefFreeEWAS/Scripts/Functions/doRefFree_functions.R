#####################################################################
# ~~~~~~~~~~~~~~~~~~
# Tumor subtype and cell type independent DNA methylation alterations 
# associated with stage progression in invasive breast carcinoma 
# ~~~~~~~~~~~~~~~~~~
# Way, G., Johnson, K., Christensen, B. 2015
#
# This script will store all of the custom functions that are used in 
# the doRefFree analyses
#####################################################################

################################
# SubsetStage Function
################################
# This function will subset stage for the input dataframe of covariates
subsetStage <- function (data, stage) {
  tmpframe <- c()
  
  # Separate the data according to the given stage assignments
  for (i in 1:length(stage)) {
    tmp <- data[data$pathologic_stage == stage[i], ]
    # Accept only primary tumors
    tmp <- tmp[tmp$sample.type == "Primary Tumor", ]
    # Combine the info together
    tmpframe <- rbind(tmpframe, tmp)
  }
  
  # Also, take all the solid tissue normal samples in this step as well
  tmpframe <- rbind(tmpframe, data[data$sample.type == "Solid Tissue Normal", ])
  return(tmpframe)
}

################################
# customRefFree Function
################################
# This function will make the design matrices, run RefFreeEWAS, and save output files
# k is the number of dimensions; it can be prespecified if you do not wish to run EstDimRMT
customRefFree <- function (covariates, betas, age = T, findRMT = T, bootstraps = 5, 
                           save = F, iPOI = 2, k) {
  
  # subset betas
  newBeta <- betas[ ,rownames(covariates)]
  
  # Null Model
  mod0 <- model.matrix(~1, data = covariates)
  
  if (age == T) {
    # Model adjusting for age as a continuous variable
    mod <- model.matrix(~ covariates$sample.type + covariates$age.Dx)
  } else {
    mod <- model.matrix(~ covariates$sample.type)
  }

  # Obtain Bstar
  tmpBstar <- as.matrix(newBeta[ , ]) %*% mod %*% solve(t(mod)%*%mod)
  
  # Get Residual Matrix
  tmpResid <- as.matrix(newBeta[ , ])-tmpBstar%*%t(mod)
  
  # Use the function EstDimRMT in the package 'isva'
  if (findRMT == T) {
    require("isva")
    RMT <- EstDimRMT(as.matrix(tmpResid), plot = F)
    k <- RMT$dim
  } else {
    k = k
  }
  
  # Perform RefFreeEWAS
  # Generate reffree model
  rf0 <- RefFreeEwasModel(as.matrix(newBeta[ , ]), mod, k) 
  
  # Use model to bootstrap
  rfb <- BootRefFreeEwasModel(rf0, bootstraps) 
  
  # Store in list
  reffreeList <- list(rfb, rf0)
  
  ################################
  # Extract RefFreeEWAS results
  ################################
  # Summarize bootstrap results
  rfbSummary <- summary(rfb)
  
  # Extract standard error for all cpgs
  SE <- rfbSummary[ , , , 2]
  
  # Get standard errors for Beta and Bstar
  seBeta <- apply(rfb[ , ,"B", ], 1:2, sd)
  seBstar <- apply(rfb[ , ,"B*", ], 1:2, sd)
  
  # Degrees of freedom
  denDegFree <- dim(mod)[1] - dim(mod)[2]
  residDf <-denDegFree
  
  # Extract deltas
  Delta <- rf0$Bstar[ ,2] - rf0$Beta[ ,2]
  seDelta <- apply(rfb[ , , "B*", ]-rfb[ , , "B", ], 1:2, sd)
  pvDelta <- 2 * pt(-abs(Delta)/seDelta[ ,2], denDegFree)

  # Get Adjusted and Unadjusted P-values
  pvAdj <- 2 * pt( -abs(rf0$Beta[ ,iPOI]/SE[ ,iPOI, 1]) , residDf)
  pvUnadj <- 2 * pt( -abs(rf0$Bstar[ ,iPOI]/SE[ ,iPOI, 2]) , residDf)
  
  # Get adjusted and unadjusted model coefficients
  coefAdj<-rf0$Beta[ ,iPOI]
  coefUnadj<-rf0$Bstar[ ,iPOI]
  
  # Create a matrix and write to file
  results <- cbind(coefAdj, coefUnadj, pvAdj, pvUnadj, Delta, pvDelta)
  
  # Get the list ready to return
  returnlist <- list(reffreeList, results, k)
  return(returnlist)
}

################################
# findGoodQ Function
################################
# This function will return a set of genes, regions, gene regions, or cgs according to a given q value cutoff, stage, or subtype
findGoodQ <- function (qcut, AnnotationList, stage, subtypes, consider = "GeneRegion", 
                       returning = "Values") {
  # Subset the annotation list to only accept the given stages
  AnnotationList <- AnnotationList[grepl(stage, names(AnnotationList))]
  # Loop over the given subtypes and reassign the Annotation List
  WhatList <- list()
  for (i in 1:length(subtypes)) {
    tmp <- AnnotationList[grepl(subtypes[i], names(AnnotationList), fixed = T)]
    WhatList[[i]] <- as.data.frame(tmp)
    colnames(WhatList[[i]]) <- colnames(tmp[[1]])
    names(WhatList)[i] <- paste(subtypes[i], "_", stage, sep = "")
  }
  # Reassign annotation list
  AnnotationList <- WhatList
  
  # Ask if we are returning values only, or the entire dataframe; if all, then initialize a list
  if (returning == "All") {
    AllFrames <- list()
  }
  
  # Get all the significant Q values
  sigQ <- list()
  for (i in 1:length(AnnotationList)) {
    tmp <- AnnotationList[[i]]
    # Remove cgs that were previously filtered
    tmp <- tmp[!is.na(tmp$qvalues), ]
    # only accept cgs that have met the q value cutoff
    tmp <- tmp[tmp$qvalues < qcut, ]
    # assign the unique gene regions, or Genes, or CpGs, to the internal list
    if (consider == "UCSC_RefGene_Name") {
      sigQ[[i]] <- unique(tmp$UCSC_RefGene_Name)
    } else if(consider == "TargetID") {
      sigQ[[i]] <- unique(tmp$TargetID)
    } else if(consider == "UCSC_RefGene_Group") {
      sigQ[[i]] <- unique(tmp$UCSC_RefGene_Group)
    } else {
      sigQ[[i]] <- unique(tmp$GeneRegion)
    }
    # store subtype specific information in the AllFrames list
    if (returning == "All") {
      AllFrames[[i]] <- tmp
    }
  }
  
  # Get all the intersects of the significant hits according to the q val cut
  for (j in 1:length(sigQ)) {
    if (j == 1) {
      compare <- sigQ[[j]]
    } else {
      compare <- intersect(compare, sigQ[[j]])
    }
  }
  # Remove "NA NA" which are CpGs that are unmapped to gene:regions
  compare <- compare[compare != "NA NA"]
  
  if (returning != "All") {
    # Return which were in common
    return(as.character(compare))
  } else {
    FilteredFrame <- data.frame()
    for (i in 1:length(AllFrames)) {
      compareMedians <- c()
      for (j in 1:length(compare)) {
        subset <- AllFrames[[i]][AllFrames[[i]][ ,consider] %in% compare[j], ]
        compareMedians <- c(compareMedians, median(subset$qvalues))
      }
      if (i == 1) {
        FilteredFrame <- compareMedians
      } else {
        FilteredFrame <- cbind(FilteredFrame, compareMedians)
      } 
    }
    colnames(FilteredFrame) <- paste(subtypes, "medianQ", sep = "-")
    rownames(FilteredFrame) <- compare
    TotalMedian <- apply(FilteredFrame, 1, median)
    FilteredFrame <- cbind(FilteredFrame, TotalMedian)
    FilteredFrame <- FilteredFrame[order(FilteredFrame[ ,ncol(FilteredFrame)], decreasing = F), ]
  }
}

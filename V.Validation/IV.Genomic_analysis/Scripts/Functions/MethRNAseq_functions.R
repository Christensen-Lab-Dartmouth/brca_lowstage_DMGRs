#####################################################################
# ~~~~~~~~~~~~~~~~~~
# Tumor subtype and cell type independent DNA methylation alterations 
# associated with stage progression in invasive breast carcinoma 
# ~~~~~~~~~~~~~~~~~~
# Way, G., Johnson, K., Christensen, B. 2015
#
# Scripts to visualize each DMGRs effect on gene expression
#####################################################################

################################
# methSeqPlot Function
################################
# Function that plots betas against RNA seq data for the same sample for a specific gene
# Gene can be a string of multiple strings
methSeqPlot <- function (gene, betas, cg, covariates, stages, subtypes, normalExprs, method = "spearman") {
  # initialize connection to cbp
  require("cgdsr")
  cbiop <- CGDS("http://www.cbioportal.org/public-portal/")
  
  # We are interested in the TCGA brca_tcga_pub2015 Study Data
  study <- getCancerStudies(cbiop)[grep("brca_tcga_pub2015", getCancerStudies(cbiop)[,1]), 1:2][1]
  cases <- unlist(strsplit(getCaseLists(cbiop, study$cancer_study_id)[1,5], " "))
  
  # Case list ID and Profile data hold the same brca_tcga_pub2015_rna_seq_v2_mrna name
  case.listID <- getCaseLists(cbiop, study$cancer_study_id)[grep("brca_tcga_pub2015_rna_seq_v2_mrna", getCaseLists(cbiop, study$cancer_study_id)[ ,1]), 1]
  profileData <- getGeneticProfiles(cbiop, study)[grep("brca_tcga_pub2015_rna_seq_v2_mrna", getGeneticProfiles(cbiop, study)[, 1]), 1][2]
  
  # Get the RNA seq data for the specific gene, or genes, of interest
  if (length(gene) == 1) {
    brca.express <- getProfileData(cbiop, gene, profileData, case.listID)
    
    # Append Normal Info
    NormalGene <- t(normalExprs[grep(paste("\\<", gene, "\\>", sep = ""), rownames(normalExprs)), ])
    brca.express <- rbind(brca.express, NormalGene)
    brca.express <- cbind(rownames(brca.express), brca.express)
  } else {
    brca.express <- list()
    for (i in 1:length(gene)) {
      brca.express[[i]] <- getProfileData(cbiop, gene[i], profileData, case.listID)
      
      # Append Normal Info
      NormalGene <- t(normalExprs[grep(paste("\\<", gene, "\\>", sep = ""), rownames(normalExprs)), ])
      brca.express[[i]] <- rbind(brca.express[[i]], NormalGene)
    }
  }
  
  # Also, this must be cpg specific
  covariates$X_SAMPLE_ID <- gsub("-", ".", covariates$X_SAMPLE_ID)
  
  # Get common samples
  common <- intersect(covariates$X_SAMPLE_ID, rownames(brca.express))
  
  # Subset expression data to common samples
  brca.sub <- brca.express[rownames(brca.express) %in% common, ]
  
  # Subset covariate file with common samples
  covSub <- covariates[covariates$X_SAMPLE_ID %in% rownames(brca.sub), ]
  
  # Samples must be in same order
  covSub <- covSub[match(covSub$X_SAMPLE_ID, rownames(brca.sub)), ]
  
  # Rename covSub to "Normal-Like" tumors
  covSub$PAM50.RNAseq[covSub$sample.type == "Primary Tumor" & covSub$PAM50.RNAseq == "Normal"] <- "Normal-Like"
  covSub$PAM50.RNAseq[covSub$sample.type == "Solid Tissue Normal"] <- "Normal"
  
  # Subset the beta file
  betaSub <- t(betas[cg, rownames(covSub)])
  
  # Compile a dataframe with the information you will need to plot
  plotFrame <- as.data.frame(cbind(covSub$PAM50.RNAseq, covSub$pathologic_stage, brca.sub[ ,2], betaSub[match(rownames(covSub), rownames(betaSub)), ]), 
                             stringsAsFactors = F)
  rownames(plotFrame) <- rownames(brca.sub)
  colnames(plotFrame) <- c("Subtype", "Stage", "RNAseq", "Beta")
  plotFrame$Stage[plotFrame$Subtype == "Normal"] <- "Normal"
  
  # Reorder the factor
  plotFrame$Stage <- factor(plotFrame$Stage, levels = c("low", "high", "Normal"))
  
  # Get correlations ready for output (Remove Normal-Like Tumors)
  outputCorReadyFrame <- plotFrame[plotFrame$Subtype != "Normal-Like", ]
  
  useSubtypes <- c("LumA", "LumB", "Her2", "Basal")
  
  correlations <- samplesize <- probability <- matrix(NA, nrow = length(unique(outputCorReadyFrame$Subtype)), ncol = 2)
  for (i in 1:length(useSubtypes)) {
    for (j in 1:2) {
      if (useSubtypes[i] != "Normal") {
        
        # Subset the plot frame according to the subtype and stage
        subset <- outputCorReadyFrame[outputCorReadyFrame$Subtype == useSubtypes[i] & outputCorReadyFrame$Stage == as.character(paste(unique(outputCorReadyFrame$Stage)[j])), ]
        
        # Add the correlations to the matrix
        correlations[i, j] <- round(cor(as.numeric(paste(subset$RNAseq)), as.numeric(paste(subset$Beta)), method = method), 2)
        
        # Get the probabilities
        pr <- cor.test(as.numeric(paste(subset$RNAseq)), as.numeric(paste(subset$Beta)), method = method)
        
        # Sore in tables
        probability[i, j] <- round(pr$p.value, 2)
        samplesize[i, j] <- nrow(subset)
      }
    }
    # Subset the info
    subset <- outputCorReadyFrame[outputCorReadyFrame$Subtype == "Normal", ]
    
    # Add Normal sample correlations, probabilities and sample size
    correlations[5, 1] <- round(cor(as.numeric(paste(subset$RNAseq)), as.numeric(paste(subset$Beta)), method = method), 2)
    pr <- cor.test(as.numeric(paste(subset$RNAseq)), as.numeric(paste(subset$Beta)), method = method)
    probability[5, 1] <- round(pr$p.value, 2)
    samplesize[5, 1] <- nrow(subset)
    
    allTumor <- plotFrame[plotFrame$Stage != "Normal", ]
    
    # Add Full correlations, probabilities and sample size
    correlations[5, 2] <- round(cor(as.numeric(paste(allTumor$RNAseq)), as.numeric(paste(allTumor$Beta)), 
                                    method = method), 2)
    pr <- cor.test(as.numeric(paste(allTumor$RNAseq)), as.numeric(paste(allTumor$Beta)), method = method)
    probability[5, 2] <- round(pr$p.value, 2)
    samplesize[5, 2] <- nrow(allTumor)
  }
  colnames(correlations) <- colnames(samplesize) <- colnames(probability) <- c("low", "high")
  rownames(correlations) <- rownames(samplesize) <- rownames(probability) <- c(useSubtypes, "Normal_and_All")
 
  outputFrame <- cbind(correlations, probability, samplesize)
  colnames(outputFrame) <- c("low_cor", "high_cor", "low_p", "high_p", "low_n", "high_n")
  
  # Get All data except Normal samples (And Normal-Like) ready to plot in bottom grid
  modelData <- data.frame(outputCorReadyFrame[outputCorReadyFrame$Subtype != "Normal", ], stringsAsFactors = F)
  
  # Get rest of data together (Normal Samples and a full set)
  allData <- modelData
  allData$Subtype <- allData$Stage <- "All Tumor"
  allData <- rbind(outputCorReadyFrame[outputCorReadyFrame$Subtype == "Normal", ], allData)
  
  # Setup plotting window
  par(mfrow=c(3, 2), oma = c(5,4,2,2), mar = c(1.8, 1.8, 2, 2))
  
  # Loop over model data and output a scatter plot for each subtype
  colorSubtypes <- c("blue", "green", "purple", "red")
  otherData_category <- c("All Tumor", "Normal")
  otherColor <- c("orange", "black")
  
  for (new_total_plots in 1:length(otherData_category)) {
    plotData <- allData[allData$Subtype == otherData_category[new_total_plots], ]
    fitPlot <- lm(as.numeric(paste(plotData$Beta)) ~ as.numeric(paste(plotData$RNAseq)))
    plot(plotData$RNAseq, plotData$Beta, pch = 16, col = otherColor[new_total_plots], xlab = "", ylab = "", cex.axis = 1.5)
    abline(fitPlot, col = otherColor[new_total_plots], lwd = 2)
  }
  
  for (new_plot in 1:length(useSubtypes)) {
    plotData <- modelData[modelData$Subtype == useSubtypes[new_plot] & modelData$Stage == "low", ]
    fitPlot <- lm(as.numeric(paste(plotData$Beta)) ~ as.numeric(paste(plotData$RNAseq)))
    plot(plotData$RNAseq, plotData$Beta, pch = 16, col = colorSubtypes[new_plot], xlab = "", ylab = "", cex.axis = 1.5)
    abline(fitPlot, col = colorSubtypes[new_plot], lwd = 2)
  }

  title(xlab = "RNAseq Expression", ylab = "Methylation Beta Values", outer = T, line = 1.5, cex.lab = 2)
  title(main = paste(gene, "Expression vs.", cg, "Methylation"), line = 0, outer = T, cex.main = 2)
  
  return(outputFrame)
}
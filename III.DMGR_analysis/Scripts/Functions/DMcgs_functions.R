#####################################################################
# ~~~~~~~~~~~~~~~~~~
# Tumor subtype and cell type independent DNA methylation alterations 
# associated with stage progression in invasive breast carcinoma 
# ~~~~~~~~~~~~~~~~~~
# Way, G., Johnson, K., Christensen, B. 2015
#
# This script will store all of the custom functions that are used in 
# the DMGR analyses
#####################################################################

################################
# SubsetStage Function
################################
# This function will subset stage for the input dataframe of covariates
subsetStage <- function (data, stage) {
  tmpframe <- c()
  for (i in 1:length(stage)) {
    tmp <- data[data$pathologic_stage == stage[i], ]
    tmp <- tmp[tmp$sample.type == "Primary Tumor", ]
    tmpframe <- rbind(tmpframe, tmp)
  }
  tmpframe <- rbind(tmpframe, data[data$sample.type == "Solid Tissue Normal", ])
  return(tmpframe)
}

################################
# Obtain qvalues Function
################################
qvalList <- function (beta, subtype, stage, covariates, qvalcut = qvalcut, betacut = "No", filelist, type = "Test") {
  if (type == "Test") {
    # Subset covariate data based on stage assignment
    if (stage == "low") {
      stageSub <- c("Stage I", "Stage IA", "Stage IB", "Stage II", "Stage IIA", "Stage IIB")
    } else if (stage == "high") {
      stageSub <- c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC", "Stage IV")
    } else if (stage == "stageI") {
      stageSub <- c("Stage I", "Stage IA", "Stage IB")
    } else if (stage == "stageII") {
      stageSub <- c("Stage II", "Stage IIA", "Stage IIB")
    } else if (stage == "stageIII") {
      stageSub <- c("Stage IIIA", "Stage IIIB", "Stage IIIC")
    } else if (stage == 'stageIV'){
      stageSub <- c("Stage IV")
    } else if (stage == 'all'){
      stageSub <- c("Stage I", "Stage IA", "Stage IB", "Stage II", "Stage IIA", "Stage IIB", "Stage IIIA", "Stage IIIB", "Stage IIIC", "Stage IV")
    } else {
      print("enter valid stage indicator"); stop()
    }
    
    # This function is taken from "doRefFree_functions.R"
    stageCov <- subsetStage(covariates, stageSub)
    
    # subset betas
    newBeta <- beta[ ,rownames(covariates)]
    
    ################################
    # Load qvalue data
    ################################
    # load qvalue
    pattern <- paste(subtype, "_", stage, "_qvalues*", sep = "")
    
    # load delta
    pattern_d <- paste(subtype, "_", stage, "_delta*", sep = "")
    files <- list.files(filelist)
    qfiles <- files[grepl(pattern, files)]
    dfiles <- files[grepl(pattern_d, files)]
    
    # read in delta values
    tmp_d <- read.table(paste(filelist, dfiles, sep = ""), stringsAsFactors = F, sep = ",")
    
    # Combine q values and subset according to cutoffs
    qlist <- list()
    for (i in 1:length(qfiles)) {
      # read in q values, beta values
      tmp <- read.table(paste(filelist, qfiles[i], sep = ""), stringsAsFactors = F, sep = ",")
      
      # Perform Filtering Steps
      tmp <- tmp[tmp[ ,2] <= as.numeric(paste(qvalcut)), ]
      if (betacut != "No") {
        tmp <- tmp[abs(tmp[ ,3]) >= as.numeric(paste(betacut)), ]
      }
      # subset tmp_d with the significant cgs
      tmp_dsub <- tmp_d[rownames(tmp), ]
      # combine into a single dataframe
      tmp <- cbind(tmp, tmp_dsub)
      
      qlist[[i]] <- tmp
      names(qlist)[i] <- unlist(strsplit(qfiles[i], "[.]"))[1]
    }
    bhInfo <- list(qlist, stageCov)
    return(bhInfo)
    
    # However, the function can also handle the Validation set analysis
  } else if(type == "Validation"){
    
    # Load model specific q values associated with single CpGs
    unadjustedQ <- read.table("V.Validation/Data/Validation_set_qvalues_unadjusted.csv", sep = ",", stringsAsFactors = F)
    adjustedQ <- read.table("V.Validation/Data/Validation_set_qvalues_adjusted.csv", sep = ",", stringsAsFactors = F)
    
    # Load model specific delta values associated with single CpGs
    Val.Deltas <- read.table("V.Validation/Data/Validation_set_delta.csv", sep = ",", stringsAsFactors = F)
    
    # Combine deltas to adjusted and unadjusted models
    unadjustedQ <- cbind(unadjustedQ, Val.Deltas)
    adjustedQ <- cbind(adjustedQ, Val.Deltas)
    
    qlist <- list("unadjusted" = unadjustedQ, "adjusted" = adjustedQ)
    
    filter.qlist <- list()
    for (i in 1:length(qlist)) {
      tmp <- qlist[[i]]
      
      # Perform Filtering Steps
      tmp <- tmp[tmp[ ,2] <= as.numeric(paste(qvalcut)), ]
      
      # assign to filter.qlist
      filter.qlist[[i]] <- tmp
      names(filter.qlist)[i] <- names(qlist)[i]
    }
    return(filter.qlist)
  }
}

################################
# splitGene() Function
################################
# This function will split the USCS_RefGene_Name by ";" and will store the first and last genes listed for the specific cpg
splitGene <- function (Gene) {
  g <- unlist(strsplit(Gene, ";"))
  result <- c(g[1], g[length(g)])
  return(result)
}

################################
# getGeneInfo() Function
################################
# This function will take the previous function and apply it to each row extracting a vector of first gene name, 
# gene region of methylation, the associated q value, and the first and last gene names
getGeneInfo <- function (row){
  use <- c(splitGene(row$UCSC_RefGene_Name)[1], splitGene(row$UCSC_RefGene_Group)[1], 
           row$q, row$beta, row$Delta, row$TargetID)
  return(use)
}

################################
# parseMethDir() Function
################################
# This function will parse together the direction of the methylation change for each given cpg
parseMethDir <- function(methdir) {
  len <- length(methdir)
  if (len == 1) {
    this <- "methdir[1]"
  } else {
    #Build the paste input depending on how long it is
    for (i in 1:len) {
      if (i == 1) {
        this <- paste("paste(", "methdir", "[",i,"]", ",", sep = "")
      } else if(i != len){
        this <- paste(this, paste("methdir", "[", i, "]", ",", sep = ""))
      } else {
        this <- paste(this, paste("methdir", "[", i, "])", sep = ""))
      }
    }
  }
  return(eval(parse(text=this)))
}

################################
# collapseInfo() Function
################################
# This function will input a matrix of genes and the associated subsettable column and output a matrix of specific genes associated with cpgs and methylated regions with associated q values
collapseInfo <- function (genes, column, annotation) {
  newF <- c()
  for (i in 1:length(unique(genes[ ,column]))) {
    cons <- genes[grepl(unique(genes[ ,column])[i], genes[ ,column]),]
    
    if (class(cons) == "matrix") {
      # get median q value for matrix
      medq <- median(as.numeric(paste(cons[ ,2])))
      # get median beta coefficient
      medb <- median(as.numeric(paste(cons[ ,3])))
      # get median delta
      medd <- median(as.numeric(paste(cons[ ,4])))
      len <- nrow(cons)
      mdir <- cons[ ,6]
      cpgs <- cons[ ,5]
      denom <- length(unique(annotation[annotation$GeneRegion == cons[1, 1], ]$TargetID))
      cpgUse <- c()
      for (j in 1:length(cpgs)) {
        if (j == 1) {
          cpgUse <- cpgs[j]
        } else {
          cpgUse <- paste(cpgUse, cpgs[j], sep = ";")
        }    
      }
      
      use <- c(cons[1, ][1], len, parseMethDir(mdir), medq, medb, medd, cpgUse, denom)
    } else {
      medq <- as.numeric(paste(cons[2]))
      medb <- as.numeric(paste(cons[3]))
      medd <- as.numeric(paste(cons[4]))
      len <- 1
      mdir <- cons[6]
      cpgUse <- cons[5]
      denom <- length(unique(annotation[annotation$GeneRegion == cons[1], ]$TargetID))
      use <- c(cons[1], len, parseMethDir(mdir), medq, medb, medd, cpgUse, denom)
    }
    
    newF <- rbind(newF, use)
  }
  colnames(newF) <- c("Gene Region", "nCpG", "Sign", "medQval", "medBetaCoef", "medDelta","cpgs", "denominator")
  rownames(newF) <- 1:nrow(newF)
  return(newF)
}
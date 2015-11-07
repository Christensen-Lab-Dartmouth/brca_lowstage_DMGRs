#####################################################################
# ~~~~~~~~~~~~~~~~~~
# Tumor subtype and cell type independent DNA methylation alterations 
# associated with stage progression in invasive breast carcinoma 
# ~~~~~~~~~~~~~~~~~~
# Way, G., Johnson, K., Christensen, B. 2015

# This script will summarize samples and inclusion as well as relevant CpGs
#####################################################################
################################
# Load Data
################################
# Load Annotation File: This annotation file was curated by Andy Houseman
annotation <- read.table("I.Data_Processing/Files/annotationfile.tsv", stringsAsFactors = F, row.names = 2, header = T, sep = "\t", nrows = 500000, comment.char = "")[,-1]

# Load total covariate file
covariatesFull <- read.table("I.Data_Processing/Files/BRCAtarget_covariates.csv", row.names = 1, header = T, sep = ",", stringsAsFactors = F)

################################
# Load Function
################################
# This function will count the percentages of a prespecified covariate name and 
# which variables of the covariate name are to be counted
countingCov <- function (covariateName, variable) {
  count <- c()
  for (i in 1:length(variable)) {
    build <- c()
    # if the variable is prespecified before hand, (E.g. if a list is given for multiple values)
    if (class(variable) == "list") {
      # begin loop over list elements
      for (j in 1:length(variable[[i]])) {
        subset <- covariates[covariates[ ,covariateName] == variable[[i]][[j]], ]
        if(j == 1) {
          # build a variable of all the subsets that belong to the specific variable
          build <- subset
        } else {
          build <- rbind(build, subset)
        }
      }
      # Begin Counting
      
      # how many samples are in the build
      size <- nrow(build)
      # find the percentage of the total samples within the covariate file
      percent <- (size / nrow(covariates)) * 100
      # get the results in the proper format
      add <- paste(size, " (", round(percent, 0), "%)", sep = "")
      # build the variable to accomodate all the elements in the variable list
      count <- c(count, add)
      
    } else { # if the variable is not a list and just a string
      
      # perform subset step and get the data in the proper format
      subset <- covariates[covariates[ ,covariateName] == variable[i], ]
       
      # Begin Counting
      size <- nrow(subset)
      percent <- (size / nrow(covariates)*100)
      add <- paste(size, " (", round(percent, 0), "%)", sep = "")
      count <- c(count, add)
    }
  }
  
  # name the variables
  if (class(variable) == "list") {
    names(count) <- names(variable)
  } else {
    names(count) <- variable
  }
  return(count)
}

################################
# Load Constants
################################
# initialize several constants; subtypes and covariate files to consider.
pam50 <- c("Basal", "Her2", "LumA", "LumB", "Normal", "MissingPAM50", "Total")
cov <- c("Age", "Stage", "sample.type")

# stages will be a list of elements that all describe a given stage number
stages <- list(stageI = c("Stage I", "Stage IA", "Stage IB"), 
               stageII = c("Stage II", "Stage IIA", "Stage IIB"), 
               stageIII = c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"), 
               stageIV = c("Stage IV"))

# We know what we want our table 1 to look like, so begin building it
Total <- matrix(NA, nrow = 10, ncol = length(pam50))
for (subtype in 1:length(pam50)) {
  # start by subsetting the covariateFull file into the subtype specific covariate file
  if (pam50[subtype] == "Total") {
    # consider total set
    covariates <- covariatesFull
    
    # Describe samples with Missing PAM50 assignments
  } else if(pam50[subtype] == "MissingPAM50") {
    # consider missing data
    covariates <- covariatesFull[covariatesFull$PAM50.RNAseq == "", ]
    
  } else {
    # consider specific PAM50 subtype
    covariates <- covariatesFull[covariatesFull$PAM50.RNAseq == pam50[subtype], ]
  }
  
  # Loop over the covariate information with the prespecified data
  for (i in 1:length(cov)) {
    # get the total count
    totCount <- nrow(covariates)
    Total[5, subtype] <- totCount
    
    # loop through the covariate information given and get data in the right format for output
    if (cov[i] == "Age") {
      mean <- round(mean(covariates$age.Dx, na.rm = T), 1)
      sd <- round(sd(covariates$age.Dx, na.rm = T), 1)
      # Add to the matrix along the way
      Total[1, subtype] <- paste(mean, " (", sd, ")", sep = "")
      
      } else if (cov[i] == "Stage") {
      stageInfo <- countingCov("pathologic_stage", stages)
      Total[6:9, subtype] <- t(stageInfo)
      missingStage <- nrow(covariates[covariates$pathologic_stage == "" | covariates$pathologic_stage == "Stage X" | covariates$pathologic_stage == "[Discrepancy]",])
      Total[10, subtype] <- paste(missingStage, " (", round((missingStage / nrow(covariates)*100), 0), "%)", sep = "")             
   
       } else if (cov[i] == "sample.type") {
      typeInfo <- countingCov("sample.type", c("Primary Tumor", "Metastatic", "Solid Tissue Normal"))
      Total[2:4, subtype] <- t(typeInfo)
    }
  }
}

# Total File row and column names
colnames(Total) <- pam50
rownames(Total) <- c("Age\nMean (SD)", "Primary Tumor", "Metastatic", "Solid Tissue Normal", "Total Samples",
                     "Stage I", "Stage II", "Stage III", "Stage IV","Stage Missing")

# write file to disk
write.table(Total, "I.Data_Processing/Tables/Table1_SampleInfo.csv", row.names = T, col.names = NA, sep = ",")

################################
# Build Table for CpGs
################################
# Read in an example dataset that has all of the cpgs that we used in the experiment
used <- read.table("I.Data_Processing/Data/BRCAmethSubset_PAM50.RNAseq_Normal.tsv", 
                   row.names = 1, header = T, sep = "\t", stringsAsFactors = F, 
                   nrows = 400000, comment.char = "")

# The CpGs that we removed
removed <- setdiff(rownames(annotation), rownames(used))

# Get Specific Information
TotalCGs <- nrow(annotation) # 485,512
TotalRemoved <- length(removed) # 96,740

# How many were removed from High Detection P Values (CpG must fail in more than 25% of samples with a p value greater than 0.00001)
DetectionP <- 2932

# Also remove Chen Probes and Sex specific Probes
Other <- TotalRemoved - DetectionP

# Prepare a vector to write to file
CGRemoval <- c(TotalCGs, nrow(used), TotalRemoved, DetectionP, Other)
names(CGRemoval) <- c("Total CpGs", "CpGs after Filtering", "Total Removed", "High Detection P Value", 
                      "Chen and Sex Specific Probes")

# write CpG file to disk
write.table(t(CGRemoval), "I.Data_Processing/Tables/Table1_CpGsRemoved.csv", sep = ",", row.names = F)
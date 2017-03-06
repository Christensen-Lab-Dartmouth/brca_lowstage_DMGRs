#####################################################################
# ~~~~~~~~~~~~~~~~~~
# Tumor subtype and cell type independent DNA methylation alterations 
# associated with stage progression in invasive breast carcinoma 
# ~~~~~~~~~~~~~~~~~~
# Way, G., Johnson, K., Christensen, B. 2015
# 
# The script will plot the TCGA tumor purity estimate stratified by 
# tumor stage and PAM50 classification
#####################################################################
################################
# Load and Subset Data
################################
library(reshape2)
library(ggplot2)

# Load and subset covariate file
covariates <- read.table("I.Data_Processing/Files/BRCAtarget_covariates.csv", row.names = 1, 
                         header = T, sep = ",", stringsAsFactors = F)

# We are only interested in summarizing Primary Tumors
covariates <- covariates[covariates$sample.type == "Primary Tumor",]
covariates <- covariates[covariates$pathologic_stage != "" & covariates$pathologic_stage != "[Discrepancy]" & covariates$pathologic_stage != "Stage X" & covariates$PAM50.RNAseq != "",]
# covSubset <- covariates[ ,c(2,5,8,9,10)]
covSubset <- covariates[ ,c(3,6,9,10,11)]

# Remove normal tumors
covSubset <- covSubset[covSubset$PAM50.RNAseq != "Normal", ]

# Make a list of low and high stage identifiers
low <- c("Stage I", "Stage IA", "Stage IB", "Stage II", "Stage IIA", "Stage IIB")
high <- c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC", "Stage IV")
stage <- list("low" = low, "high" = high)

# Rename the covariate subset frame
for(i in 1:length(stage)) {
  for(j in 1:length(stage[[i]])) {
    covSubset$pathologic_stage[covSubset$pathologic_stage == stage[[i]][j]] <- names(stage)[i]
  }
}

################################
# Create Summary Table
################################
# constants we are interested in plotting
interest <- c("percent_tumor_cells_TOP", "percent_stromal_cells_TOP", "percent_normal_cells_TOP")
subtypes <- c("Basal", "Her2", "LumA", "LumB")
stagesub <- c("high", "low")

# Output a summary score for the constants of interest to get ready for writing to file
summary <- c()
for (i in 1:length(subtypes)) {
  # Subset the covariate file by PAM50 subtype
  cov <- covSubset[covSubset$PAM50.RNAseq == subtypes[i], ]
  stageInfo <- c()
  
  # For each stage, get respective information
  for (j in 1:length(stagesub)) {
    covsub <- cov[cov$pathologic_stage == stagesub[j], ]
    n <- nrow(covsub)
    
    # Get the mean and sd of the percentages of interest
    mean <- c()
    sd <- c()
    for (k in 1:length(interest)) {
      pm <- paste(round(mean(covsub[ ,interest[k]], na.rm = T),2), "%", sep = "")
      psd <- paste(round(sd(covsub[ ,interest[k]], na.rm = T),2), "%", sep = "")
      
      mean <- c(mean, pm)
      sd <- c(sd, psd)
    }
    
    # Compile together into a usable format
    tmp <- c(n, mean, sd)
    stageInfo <- rbind(stageInfo, tmp)
  }
  
  summary <- rbind(summary, stageInfo)
}

# reorder the summary file
summary <- cbind(summary[,1], summary[,2], summary[,5], summary[,3], summary[,6], summary[,4], summary[,7])

rownames(summary) <- c(paste(rep(subtypes[1],2), stagesub), paste(rep(subtypes[2], 2), stagesub), 
                       paste(rep(subtypes[3],2), stagesub), paste(rep(subtypes[4],2), stagesub))
colnames(summary) <- c("n", "Tumor Cells (Mean)", "Tumor Cells (SD)", "Stromal Cells (Mean)", "Stromal Cells (SD)", "Normal Cells (Mean)", "Normal Cells (SD)")

write.table(summary, file = "I.Data_Processing/Tables/covariateTumorPuritySummary.csv", sep = ",", row.names = T, col.names = NA)

################################
# Build BarChart
################################
# First, build a dataframe to melt to get ready for ggplot
# ggplotCov <- covariates[ ,c(2,5,8,9,10)]
ggplotCov <- covariates[ ,c(3,6,9,10,11)]

# only consider samples with complete information
ggplotCov <- ggplotCov[complete.cases(ggplotCov), ]
ggplotCov <- ggplotCov[ggplotCov$PAM50.RNAseq != "Normal", ]

# Load Barchart constants
low <- c("Stage I", "Stage IA", "Stage IB", "Stage II", "Stage IIA", "Stage IIB")
high <- c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC", "Stage IV")
stage <- list("Low" = low, "High" = high)

# Rename all classifiers either "low" or "high" stage
for (i in 1:length(stage)) {
  for (j in 1:length(stage[[i]])) {
    ggplotCov$pathologic_stage[ggplotCov$pathologic_stage == stage[[i]][j]] <- names(stage)[i]
  }
}

# plot ggplot (Supplemental Figure S1)
ggplot(data = melt(ggplotCov), aes(x = variable, y = as.numeric(paste(value)))) + 
  geom_boxplot(aes(fill = PAM50.RNAseq, colour = pathologic_stage)) + 
  scale_x_discrete(labels = c("Tumor Cells", "Stromal Cells", "Normal Cells")) + 
  scale_colour_manual(values = c("gold", "green"), labels = c('Late', 'Early')) + 
  scale_fill_manual(values = c("red", "pink", "blue", "cyan"),
                    labels = c('Basal-Like', 'Her2', 'Luminal A', 'Luminal B')) +
  xlab("") + ylab("Percentage") + 
  ggtitle("Summary of Tumor Purity") + theme_bw()

ggsave(filename = paste("I.Data_Processing/Figures/covariateSummaryBarPlot.png", sep = ""), width = 10, height = 8)

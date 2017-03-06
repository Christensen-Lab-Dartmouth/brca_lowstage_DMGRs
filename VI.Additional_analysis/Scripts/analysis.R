#####################
# CNV analysis
#####################
dir = 'C:\\Users\\atitus\\Documents\\github\\brca_lowstage_DMGRs'
setwd(dir)
data = read.csv('VI.Additional_analysis\\Data\\CNA_Genes_data.csv')
data_sub = data[grep('^1p', data$Cytoband), ]
data_sub = data_sub[order(data_sub$Cytoband, decreasing = T), ]
data_sub$X. = NULL
write.csv(data_sub, file = "VI.Additional_analysis\\Data\\CNA_Genes_data_1p.csv")


#####################
# Prepare data for the Manhattan Plot
#####################
#install.packages('qqman')
library(qqman)
library(data.table)
results_dir = "II.RefFreeEWAS/Data/"
files = list.files(results_dir)
results_files = files[grep("_low_qvalues_adjusted", files)]
results_files = results_files[-(grep('^Normal', results_files))]
basal = read.csv(paste(results_dir, results_files[grep("^Basal", results_files)], sep = "/"))
her2 = read.csv(paste(results_dir, results_files[grep("^Her", results_files)], sep = "/"))
luma = read.csv(paste(results_dir, results_files[grep("^LumA", results_files)], sep = "/"))
lumb = read.csv(paste(results_dir, results_files[grep("^LumB", results_files)], sep = "/"))
full = cbind(basal$pvalues, her2$pvalues, luma$pvalues, lumb$pvalues)
rownames(full) = rownames(basal)
full = cbind(full, rowMeans(full))
colnames(full) = c('BasalP', "Her2P", "LumAP", "LumBP", 'MeanP')
temp_data = as.data.frame(full[, "MeanP"])
colnames(temp_data) <- "MeanP"
#temp_data$negLogP = -log(temp_data$MeanP)
temp_data$IlmnID = rownames(temp_data)
annot_data = as.data.frame(fread('I.Data_Processing/Files/Reduced_annotationfile.csv'))
annot_data = annot_data[order(annot_data, decreasing = T), ]
merged_data = merge(annot_data, temp_data, by= "IlmnID")
#merged_data$MeanP = NULL
colnames(merged_data) = c('SNP', 'CHR', 'BP', 'P')
merged_data = merged_data[merged_data$CHR != "MULTI", ]
merged_data$CHR = as.numeric(merged_data$CHR)
merged_data$BP = as.numeric(merged_data$BP)
#####################
# Create the Manhattan Plot
#####################
png("VI.Additional_analysis/Figures/EWAS_manhattan_lowStage.png", width = 12, height = 8, units = 'in', res = 300)
manhattan(merged_data)
dev.off()
#####################
# Create the qq-plot Plot
#####################
#png("VI.Additional_analysis/Figures/EWAS_qqplot_lowStage.png", width = 12, height = 8, units = 'in', res = 300)
#qq(merged_data$P, main = "Q-Q plot of EWAS p-values")
#dev.off()


#####################
# Create bar plot subtype by -log(Q)
#####################
data = read.csv('VI.Additional_analysis/Data/DMGR_Q_Positions.csv')
data$nLogQ = -log(data$AvgMedQ)
library(ggplot2)
png("VI.Additional_analysis/Figures/medQ_by_genomic_position_all.png", width = 12, height = 8, units = 'in', res = 300)
ggplot(data, aes(fill=Subtype, y=nLogQ, x=Position)) +
  geom_bar(position="dodge", stat="identity") + 
  scale_x_continuous(breaks = unique(data$Position), labels=unique(data$DMGR)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, max(temp_data$nLogQ))) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "DMGR spaced by genomic position", y= "-log(Q)") +
  geom_hline(aes(yintercept = 2), color = "grey") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
dev.off()

temp_data = data[data$Position > 15000, ]

png("VI.Additional_analysis/Figures/medQ_by_genomic_position.png", width = 12, height = 8, units = 'in', res = 300)
ggplot(temp_data, aes(fill=Subtype, y=nLogQ, x=Position)) +
  geom_bar(position="dodge", stat="identity") + 
  scale_x_continuous(breaks = unique(temp_data$Position), labels=unique(temp_data$DMGR)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, max(temp_data$nLogQ))) + 
  theme(axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1)) +
  labs(x = "DMGR spaced by genomic position", y= "-log(Q)") +
  geom_hline(aes(yintercept = 2), color = "grey") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
dev.off()

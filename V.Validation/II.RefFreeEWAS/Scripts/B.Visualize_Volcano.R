#####################################################################
# ~~~~~~~~~~~~~~~~~~
# Tumor subtype and cell type independent DNA methylation alterations 
# associated with stage progression in invasive breast carcinoma 
# ~~~~~~~~~~~~~~~~~~
# Way, G., Johnson, K., Christensen, B. 2015

# The script will output volcano plots for all models adjusted and unadjusted
#####################################################################

################################
# Load Libraries
################################
library(qvalue)
library(ggplot2)
library(gridExtra)

################################
# Load Constants
################################
subtype <- c("Basal", "Her2", "LumA", "LumB", "Normal")
stage <- c("low", "high")

# Loop over each subtype and stage combo to output volcano plots
for (subt in subtype) {
  for (sta in stage) {
    
    # Read in data
    ad_file <- paste(subt, sta, "qvalues", "adjusted.csv", sep = "_")
    un_file <- paste(subt, sta, "qvalues", "unadjusted.csv", sep = "_")
    deltas <- paste(subt, sta, "delta.csv", sep = "_")
    
    adjusted <- read.table(paste("II.RefFreeEWAS/Data/", ad_file, sep = ""), sep = ",")
    unadjusted <- read.table(paste("II.RefFreeEWAS/Data/", un_file, sep = ""), sep = ",")
    deltas <- read.table(paste("II.RefFreeEWAS/Data/", deltas, sep = ""), sep = ",")
    
    ################################
    # Visualize RefFreeEWAS results
    ################################
    # Append adjusted Q values to the results file
    qtmp <- qvalue(adjusted$pvalues, fdr.level = 0.01)
    qtmp_high <- qvalue(adjusted$pvalues, fdr.level = 0.05)
    
    # Get Q value cutoffs to draw lines in plots
    tmp1 <- adjusted[qtmp$significant == T, ]
    tmp2 <- adjusted[qtmp_high$significant == T, ]
    
    # Of these cutoffs, what is the lowest p value?
    qcut1 <- min(-log10(as.numeric(paste(tmp1[ ,1])))) 
    qcut2 <- min(-log10(as.numeric(paste(tmp2[ ,1]))))
    
    # Find max and min for plotting margins
    maxX <- max(c(max(as.numeric(paste(adjusted[ ,3]))),as.numeric(paste(max(unadjusted[ ,3])))))
    maxY <- max(c(-log10(as.numeric(paste(unadjusted[ ,1])))[!is.infinite(-log10(as.numeric(paste(unadjusted[ ,1]))))], 
                  -log10(as.numeric(paste(adjusted[ ,1])))[!is.infinite(-log10(as.numeric(paste(adjusted[ ,1]))))]))
    
    # Plotting margin for adjusted volcano plot
    maxY1 <- max(-log10(as.numeric(paste(adjusted[ ,1])))[!is.infinite(-log10(as.numeric(paste(adjusted[ ,1]))))])
    
    # Store the plots to then arrange them on a grid to save as png
    # Unadjusted Plot
    un <- ggplot(unadjusted, aes(as.numeric(paste(unadjusted[ ,3])), -log10(as.numeric(paste(unadjusted[ ,1]))))) + 
      geom_point(aes(colour = deltas[ ,1])) + 
      scale_color_gradient2(low = "blue", mid="grey", high = "red") + 
      labs(list(x = "Beta Coefficient", y = "-log10 p Value", title = paste("Unadjusted\n", subt, sta), color = "Delta"), size = 18) + 
      geom_hline(yintercept = qcut1, color = "red", lintetype = "longdash") + 
      geom_hline(yintercept = qcut2, acolor = "black", lintetype = "longdash") +
      xlim((-1 * maxX), maxX) + ylim(0, maxY) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), panel.border = element_blank(), 
            axis.line = element_line(), axis.text = element_text(size = rel(1.4), color = "black"),
            axis.title = element_text(size = rel(1.6)), plot.title = element_text(size = rel(2)),
            legend.text = element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.3)),
            axis.ticks = element_line(size = rel(1.4), color = "black"),
            axis.line = element_line(size = rel(1.4), color = "black"))
    
    # Adjusted plot
    ad <- ggplot(adjusted, aes(as.numeric(paste(adjusted[ ,3])), -log10(as.numeric(paste(adjusted[ ,1]))))) + 
      geom_point(aes(colour = deltas[ ,1])) + 
      scale_color_gradient2(low = "blue", mid="grey", high = "red") + 
      labs(list(x = "Beta Coefficient", y = "-log10 p Value", title = paste("Reference Free Adjusted\n", subt, sta), color = "Delta")) + 
      geom_hline(yintercept = qcut1, color = "red", lintetype = "longdash") + 
      geom_hline(yintercept = qcut2, acolor = "black", lintetype = "longdash") +
      xlim((-1 * maxX), maxX) + ylim(0, maxY) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), panel.border = element_blank(), 
            axis.line = element_line(), axis.text = element_text(size = rel(1.4), color = "black"),
            axis.title = element_text(size = rel(1.6)), plot.title = element_text(size = rel(2)),
            legend.text = element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.3)),
            axis.ticks = element_line(size = rel(1.4), color = "black"),
            axis.line = element_line(size = rel(1.4), color = "black"))
    
    # Adjusted plot with different max y value
    ad1 <- ggplot(adjusted, aes(as.numeric(paste(adjusted[,3])), -log10(as.numeric(paste(adjusted[,1]))))) + 
      geom_point(aes(colour = deltas[ ,1])) + 
      scale_color_gradient2(low = "blue", mid="grey", high = "red") + 
      labs(list(x = "Beta Coefficient", y = "-log10 p Value", title = paste("Reference Free Adjusted \n", subt, sta, "(Resized)"), color = "Delta")) + 
      geom_hline(yintercept = qcut1, color = "red", lintetype = "longdash") + 
      geom_hline(yintercept = qcut2, acolor = "black", lintetype = "longdash") +
      xlim((-1 * maxX), maxX) + ylim(0, maxY1) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), panel.border = element_blank(), 
            axis.line = element_line(), axis.text = element_text(size = rel(1.4), color = "black"),
            axis.title = element_text(size = rel(1.6)), plot.title = element_text(size = rel(2)),
            legend.text = element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.3)),
            axis.ticks = element_line(size = rel(1.4), color = "black"),
            axis.line = element_line(size = rel(1.4), color = "black"))
    
    png(paste("II.RefFreeEWAS/Figures/", subt, "_", sta, "_volcano.png", sep = ""), width = 1000, height = 450)
    grid.arrange(un, ad, ad1, ncol = 3, nrow = 1)
    dev.off()
  }
}
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
#install.packages("gridExtra")
#install.packages("R.utils")
library(qvalue)
library(ggplot2)
library(grid)
library(gridExtra)
library(R.utils)

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
    
    # Generate the theme of each ggplot
    overall_theme <-  theme(panel.grid.major = element_blank(), 
                            panel.grid.minor = element_blank(), 
                            panel.background = element_blank(), 
                            panel.border = element_blank(),
                            axis.line = element_line(size = rel(4), color = "black"), 
                            axis.text = element_text(size = rel(6), color = "black"),
                            axis.ticks = element_line(size = rel(5), color = "black"),
                            axis.ticks.length = unit(0.5, "cm"),
                            axis.title = element_text(size = rel(6.8)),
                            axis.title.y = element_text(vjust = 4.5),
                            axis.title.x = element_text(vjust = -4.5),
                            plot.title = element_text(size = rel(7.5)),
                            legend.key.size = unit(2, "cm"),
                            legend.text = element_text(size = rel(4.5)), 
                            legend.title = element_text(size = rel(4.7)),
                            plot.margin = unit(c(1.2, 2, 2, 1.5), 'cm'))
    
    # Get information about titles
    if (subt == 'Basal') {
      subtype_name <- 'Basal-like'
    } else {
      subtype_name <- subt
    }
    
    if (sta == 'low') {
      stage_name <- 'Early'
    } else {
      stage_name <- 'Late'
    }
    
    
    # Store the plots to then arrange them on a grid to save as png
    # Unadjusted Plot
    un <- ggplot(unadjusted, aes(as.numeric(paste(unadjusted[ ,3])), 
                                 -log10(as.numeric(paste(unadjusted[ ,1]))))) + 
      geom_point(aes(colour = deltas[ ,1]), size = 9) + 
      scale_color_gradient2(low = "blue", mid = "grey", high = "red") + 
      labs(list(x = "Beta Coefficient", 
                y = "-log10 p Value", 
                title = paste("Unadjusted\n", subtype_name, "-", stage_name), 
                color = "Delta")) + 
      geom_hline(yintercept = qcut1, color = "red", linetype = "dashed", size = 1.8) + 
      geom_hline(yintercept = qcut2, color = "black", linetype = "dashed", size = 1.8) +
      xlim((-1 * maxX), maxX) + ylim(0, maxY) + overall_theme
    
    # Adjusted plot
    ad <- ggplot(adjusted, aes(as.numeric(paste(adjusted[ ,3])), 
                               -log10(as.numeric(paste(adjusted[ ,1]))))) + 
      geom_point(aes(colour = deltas[ ,1]), size = 9) + 
      scale_color_gradient2(low = "blue", mid = "grey", high = "red") + 
      labs(list(x = "Beta Coefficient", 
                y = "-log10 p Value", 
                title = paste("Reference Free Adjusted\n", subtype_name, "-", stage_name), 
                color = "Delta")) + 
      geom_hline(yintercept = qcut1, color = "red", linetype = "dashed", size = 1.8) + 
      geom_hline(yintercept = qcut2, color = "black", linetype = "dashed", size = 1.8) +
      xlim((-1 * maxX), maxX) + ylim(0, maxY) + overall_theme
    
    # Adjusted plot with different max y value
    ad1 <- ggplot(adjusted, aes(as.numeric(paste(adjusted[,3])), 
                                -log10(as.numeric(paste(adjusted[,1]))))) + 
      geom_point(aes(colour = deltas[ ,1]), size = 9) + 
      scale_color_gradient2(low = "blue", mid="grey", high = "red") + 
      labs(list(x = "Beta Coefficient", 
                y = "-log10 p Value", 
                title = paste("Adjusted (Resized)\n", subtype_name, "-", stage_name), 
                color = "Delta")) + 
      geom_hline(yintercept = qcut1, color = "red", linetype = "dashed", size = 1.8) + 
      geom_hline(yintercept = qcut2, color = "black", linetype = "dashed", size = 1.8) +
      xlim((-1 * maxX), maxX) + ylim(0, maxY1) + overall_theme
    
    png(paste("II.RefFreeEWAS/Figures/", subtype_name, "_", stage_name, "_volcano.png", sep = ""), 
        width = 4000, height = 1800)
    grid.arrange(un, ad, ad1, ncol = 3, nrow = 1)
    dev.off()
  }
}

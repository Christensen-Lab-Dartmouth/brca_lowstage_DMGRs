#####################################################################
# ~~~~~~~~~~~~~~~~~~
# Tumor subtype and cell type independent DNA methylation alterations 
# associated with stage progression in invasive breast carcinoma 
# ~~~~~~~~~~~~~~~~~~
# Way, G., Johnson, K., Christensen, B. 2015
#
# Visualize genome location and differentially methylated CpGs for 
# all specific genes
#####################################################################

################################
# Load Libraries
################################
#source("https://bioconductor.org/biocLite.R")
#biocLite("Gviz")
#biocLite("Homo.sapiens")
library(Gviz)
library(Homo.sapiens)
library(readr)
library(plyr)

################################
# Load Functions
################################
# This function will pull all the CpGs that are associated with a given gene region in the overlap file
ExtractCommonCGs <- function (gene) {
  # Get the gene info for the specific gene
  geneInfo <- CommonOverlaps[gene,]
  # Get all of the cpgs of interest regardless of the subtype
  CpGs <- c(unlist(strsplit(geneInfo[ ,5], ";")), unlist(strsplit(geneInfo[ ,10], ";")), 
            unlist(strsplit(geneInfo[ ,15], ";")), unlist(strsplit(geneInfo[ ,20], ";")))
  return(CpGs)
}

PlotGvizTracks <- function (gene, toggle = c(100,00)) {
  # Get the gene symbol and granges object for the input gene
  wh <- genesymbol[gene]
  wh <- range(wh, ignore.strand = T)
  
  # Extract all of the differentially  methylated CpGs from the given gene input
  CpGs <- unique(ExtractCommonCGs(gene))
  
  # Subset the annotation file
  Anno.Gene <- annotation[annotation$Name %in% CpGs, ]
  
  # Get the location of the significant CpGs
  CpG.loc <- Anno.Gene$MAPINFO
  
  # Find the chromosome number of the gene of interest
  chrom <- paste("chr", Anno.Gene$CHR[1], sep = "")
  
  # Identify where to start and where to end
  from <- wh@ranges@start - toggle[1]
  to <- wh@ranges@start+ wh@ranges@width + toggle[2]
  
  #######################
  # Tracks
  #######################
  # Get the ideogram track
  itr <- IdeogramTrack(genome="hg19", chromosome=chrom)
  
  # Obtain the actual genomic location id numbers
  gtr <- GenomeAxisTrack(genome="hg19", chromosome=chrom)
  
  # Pull information from UCSC 
  # Get all known genes in that range
  knownGenes <- UcscTrack(genome = "hg19", chromosome = chrom, track = "knownGene", from = from, to = to, 
                          trackType = "GeneRegionTrack", rstarts = "exonStarts", rends = "exonEnds", gene = "name", 
                          symbol = "name", transcript = "name", strand = "strand", fill = "orange", name = gene)
  
  # Get CpG island information
  cpgIslands <- UcscTrack(genome = "hg19", chromosome = chrom, track = "cpgIslandExt", from = from, to = to, 
                          trackType = "AnnotationTrack", start = "chromStart", end = "chromEnd", id = "name", 
                          shape = "box", fill = "#006400", name = "CpG Islands")
  
  # Get DNase Sensitivity information
  dnase <- UcscTrack(genome = "hg19", chromosome = chrom, track = "DNase Clusters", from = from, to = to, 
                     trackType = "DataTrack", strand = "*", name = "DNase", start="start", end="end", data="score", 
                     type="hist", window=-1, windowSize=1500, fill.histogram="blue",
                     col.histogram="black")
  
  # Draw lines where the specific CpGs are
  ht1 <- HighlightTrack(genome = "hg19", trackList = list(knownGenes, cpgIslands, dnase), start = CpG.loc, width = 1, chromosome = chrom, 
                        col = "red", inBackground = F)
  
  # Plot the tracks on the screen
  plotTracks(trackList = list(itr, gtr, ht1), from = from, to = to, sizes = c(1, 2, 4, 2, 2.5), 
                                              cex = 0.85, cex.sampleNames = 2, showTitle = T)
}

################################
# Load Data
################################
# Load extended annotation file. This file does not ignore multiple genes for a single cpg
annotation <- read_csv("I.Data_Processing/Files/HumanMethylation450K_Annotation_File.csv", skip = 7)
annotation <- as.data.frame(annotation)

# Load Common Overlaps
CommonOverlaps <- read.csv("III.DMGR_analysis/Tables/commonLowStageOverlaps_FullAnnotation_extended.csv", 
                           row.names = 1, header = T, stringsAsFactors = F)

# Get all the genes in common to low stage tumors
Genes <- laply(rownames(CommonOverlaps), function (x) {unlist(strsplit(x, " "))[1]})

# Holds information about genomic location for a given gene
data(genesymbol, package = "biovizBase")

################################
# Make Plots
################################
for (i in 1:length(Genes)) {
  png(paste("IV.Genomic_analysis/Figures/GeneLocation/GenomicLocation_", Genes[i], ".png", sep = ""), height = 320, width = 510)
  PlotGvizTracks(Genes[i])
  dev.off()
}

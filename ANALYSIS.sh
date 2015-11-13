#!bin/bash

#####################################################################
# ~~~~~~~~~~~~~~~~~~
# Identification and validation of tumor subtype and cell type independent DNA
# methylation alterations in breast carcinoma
# ~~~~~~~~~~~~~~~~~~
#
# Way, G., Johnson, K., Christensen, B. 2015
#
# The following bash script will perform a cell proportion analysis of 
# solid tumor TCGA breast cancer samples. We will use this correction 
# to analyze DNA methylation changes of low stage TCGA breast tumors 
# stratefied by PAM50 Subtype.
#####################################################################

#Install R Packages
R --no-save < INSTALL.R

####################
# Step I: Data Processing
####################
# Data can be downloaded from TCGA
# https://tcga-data.nci.nih.gov/tcga/dataAccessMatrix.htm

#~~~~GET DATA
# Data Type: DNA Methylation; Batch Number: All; Data Level: Level 1
# Availibility: Available; Preservation: All; Center/Platform: JHU_USC (HumanMethylation450)
# Access date: March 2015

# Download Clinical Data from UCSC Cancer Genome Browser: 
# https://genome-cancer.ucsc.edu/proj/site/hgHeatmap/

# NOTE: It is important to place the IDAT files and the Clinical data matrix in the same folder
IDAT_loc="../../Documents/mdata/TCGAbreast_idat/" # Place download folder location here

#~~~~~PREPROCESS DATA
# Run preprocessing script (if third argument is "summary", then the script will output a processing summary)
R --no-save --args "I.Data_Processing/Data/TCGA_BRCA" $IDAT_loc "nosummary" < I.Data_Processing/Scripts/A.preprocess_minfi.R

# Subset data to spectific subtypes; tumor adjacent samples are included in all subsets
R --no-save --args PAM50.RNAseq Basal Her2 LumA LumB Normal < I.Data_Processing/Scripts/B.subsetBetas.R

#~~~~~EXPAND ANNOTATION FILE
# This file will expand the annotation file to include multiple cgs that are associated with multiple gene:regions
R --no-save < I.Data_Processing/Scripts/C.editAnnotationFile.R # We will use this expanded annotation file in downstream analyses

# Describe the data you are using (Generate Table 1)
# To get information about CpG values, you need to have run a RefFreeEWAS already to understand which CpGs were actually used
R --no-save < I.Data_Processing/Scripts/D.DescribeTableI.R

# Describe TCGA tumor purity estimates
R --no-save < I.Data_Processing/Scripts/E.TumorPurity.R

####################
# Step II: Reference-Free cell type adjustment modeling
####################
# NOTE- This was performed on the Discovery cluster at Dartmouth College
# Apply a model correcting for age, and cell proportion heterogeneity estimates with 1000 bootstraps
# The python script runs the R script "II.RefFreeEWAS/Scripts/A.doRefFree.R"
python II.RefFreeEWAS/Scripts/Functions/RefFree_Discovery.py

# Output Volcano Plots
R --no-save < II.RefFreeEWAS/Scripts/B.Visualize_Volcano.R

# Summarize dimensionality
# First get the dimensionality of the full dataset as well as all of the tumor adjacent samples
R --no-save < II.RefFreeEWAS/Scripts/C.FullDataDimensionality.R

# Next, get all dimensionality estimates into a table
R --no-save < II.RefFreeEWAS/Scripts/D.Summarize_Dimensions.R

# Following the results of the Reference Free modeling approach, investigate and summarize q value cutoffs
R --no-save < II.RefFreeEWAS/Scripts/E.Summarize_q.R  # The script will output the 400 genes used in the DAVID Analysis

####################
# Step III: Identifying Differentially Methylated Gene Regions (DMGRs)
####################
# NOTE- This step was performed in the Discovery cluster at Dartmouth College
python III.DMGR_analysis/Scripts/Functions/runDM.py  # Script will run A.DMcgs.R

# Compare DMGRs within PAM50 across stage and across PAM50 within Stage (Output venn diagrams comparing overlaps) and low/high stage overlap tables
R --no-save < III.DMGR_analysis/Scripts/B.DMGRcomparison.R

# Perform Gene:Region Enrichment Test
R --no-save < III.DMGR_analysis/Scripts/C.FisherExact_RegionEnrichment.R

# Generate HeatMaps for all Gene:Regions and a heatmap in common to all low stage overlapping gene region associated CpGs
R --no-save < III.DMR_analysis/Scripts/D.DMGRheatmap.R 

####################
# Step IV: DMGR Gene Expression and Genomic Locations
####################
# Input the 400 genes as given by II.RefFreeEWAS/Scripts/D.Summarize_q.R [http://david.abcc.ncifcrf.gov/tools.jsp] to run a pathways analysis

# Get RNAseq data for TCGA BRCA data from Synapse
# To install synapse follow instructions listed here: https://www.synapse.org/#!Wiki:syn1768504/ENTITY
synapse get syn1446185
# Move synapse data to appropriate location
mv unc.edu* IV.Genomic_analysis/Data

# For the 12 common low stage DMGRs, explore CpG association with gene expression
R --no-save < IV.Genomic_analysis/Scripts/A.MethRNAseq.R

# This script will plot the genomic locations of all genes and the associated significant CpG methylation differences
R --no-save < IV.Genomic_analysis/Scripts/B.GenomicLocations.R

####################
# Step V: Validation
####################
# NOTE: The argument supplied will hold the IDAT files for download
IDAT_validation="../../Documents/IDAT/"

# Get the Fliescher et al 2015 Data from GEO GSE60185 
R --no-save --args $IDAT_validation < V.Validation/Scripts/A.getValidationData.R

# Remove the other .csv file that is automatically downloaded with the GEO profile (Not removing this file will cause errors in preprocessing)
rm $IDAT_validation/*.csv

# Copy the manifest file (Located in "V.Validation/Data/GSE60185_manifest.csv") to the location of the IDAT files
cp V.Validation/Data/GSE60185_manifest.csv $IDAT_validation

# Use minfi to preprocess the Fliescher Data
R --no-save --args  $IDAT_validation < V.Validation/Scripts/B.preprocess_minfi_Validation.R

# Run RefFreeEWAS on the Validation Set using 1000 bootstraps
# NOTE- This was performed on the Discovery cluster at Dartmouth College
R --no-save --args 1000 < V.Validation/Scripts/C.doRefFree_Validation.R

# Output a table of dimensionality estimation for the validation set
R --no-save < V.Validation/Scripts/D.Dimensionality_Validation.R

# Compare Validation DMGRs to Test Set DMGRs (output Scatter plot and Heatmap)
R --no-save < V.Validation/Scripts/E.Validation_DMGRComparison.R

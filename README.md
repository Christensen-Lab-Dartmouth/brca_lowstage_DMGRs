# README #

#####################################################################
# Analytical Code for "Identification and validation of tumor subtype and cell type independent DNA methylation alterations in breast carcinoma" 
Way, G., Johnson, K., Christensen, B. 2015, Manuscript in preparation
#####################################################################
[![DOI](https://zenodo.org/badge/18957/gwaygenomics/brca_lowstage_DMGRs.svg)](https://zenodo.org/badge/latestdoi/18957/gwaygenomics/brca_lowstage_DMGRs)
#######################
# SUMMARY
#######################
The following repository contains all scripts required to reproduce an analysis of low stage invasive
breast carcinoma investigating similarities in DNA methylation measured by The Cancer Genome Atlas on
the Illumina 450k platform. At the core of the analysis is a reference free adjustment of cell type
proportion (see Houseman, Molitor, and Marsit. Bioinformatics 2014) on each PAM50 subtype stratified by 
low and high stages. After the adjustment, we observe key gene regions of differential methylation in 
common to all PAM50 subtypes in the low stage. We also validate these findings in a Validation set (see 
Yang et al. Genome Biology 2015).

#######################
# CONTACT
#######################
Please report all bugs and direct scripting questions to:
Gregory.Way.GR@dartmouth.edu.

Questions regarding the analysis or other correspondance should be directed to:
Brock.C.Christensen@dartmouth.edu

#######################
# ANALYSIS
#######################
All scripts are intended to be run in order, as defined in ANALYSIS.sh. The bash script can be run 
directly to reproduce the pipeline but this is not recommended since some steps are computationally
intensive. Instead, to successfully reproduce this analysis, we recommend following the bash script line
by line. The folder structure is labelled to facilitate easy step-wise navigation, as are the Scripts/ 
folders in each parent folder. All scripts should be run from the top directory in the repository. For example: 

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Place downloaded IDAT files in this folder
IDAT_loc="../../Documents/mdata/TCGAbreast_idat/"

# Run preprocessing script (if third argument is "summary", then the script will output a processing summary)
R --no-save --args "I.Data_Processing/Data/TCGA_BRCA" $IDAT_loc "nosummary" < I.Data_Processing/Scripts/A.preprocess_minfi.R
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#######################
# DATA
#######################
Data is not stored directly in this repository and should be downloaded according to the ANALYSIS.sh script.

#######################
# DEPENDENCIES
#######################
To install all required packages from CRAN and Bioconductor please run the INSTALL.R script

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
R --no-save < INSTALL.R
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

########################
# VERSION CONTROL
########################

# R Version
* R 3.1.2

# R Packages
* BiocGenerics_0.14.0 
* Biobase_2.28.0
* Biostrings_2.36.1
* bumphunter_1.8.0
* DBI_0.3.1
* fastICA_1.2-0 
* foreach_1.4.2
* Formula_1.2-1
* GenomicFeatures_1.20.1
* GenomeInfoDb_1.4.0
* GenomicRanges_1.20.4
* GEOquery_2.34.0
* ggplot2_1.0.1
* GO.db_3.1.2
* gridExtra_2.0.0
* Gviz_1.12.0
* Hmisc_3.16-0
* Homo.sapiens_1.1.2
* IlluminaHumanMethylation450kanno.ilmn12.hg19_0.2.1
* IlluminaHumanMethylation450kmanifest_0.4.0
* IRanges_2.2.3
* isva_1.8
* iterators_1.0.7
* lattice_0.20-31
* limma_3.24.7
* locfit_1.5-9.1
* minfi_1.14.0
* org.Hs.eg.db_3.1.2
* OrganismDBi_1.10.0
* plyr_1.8.3
* readr_0.1.1
* RefFreeEWAS_1.3
* reshape2_1.4.1
* RSQLite_1.0.0
* S4Vectors_0.6.0
* survival_2.38-1
* TxDb.Hsapiens.UCSC.hg19.knownGene_3.1.2
* qvalue_2.0.0
* XVector_0.8.0

# Python Version
* Python 2.7.6

# Operating System
* Ubuntu 14.04.2 LTS
* MAC OSX

########################
# ACKNOWLEDGEMENTS
########################
This work was supported by the Institute for Quantitative Biomedical Sciences and two grants P20GM104416 and R01DE02277 (BCC)

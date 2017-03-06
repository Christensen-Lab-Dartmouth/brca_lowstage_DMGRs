# Cell-type (and subtype) independent DNA methylation alterations in breast cancer 

Titus, A.J., Way, G., Johnson, K., Christensen, B. 2016, Manuscript in preparation

[![DOI](https://zenodo.org/badge/18957/gwaygenomics/brca_lowstage_DMGRs.svg)](https://zenodo.org/badge/latestdoi/18957/gwaygenomics/brca_lowstage_DMGRs)

## Summary 
The following repository contains all scripts required to reproduce an analysis
of low stage invasive breast carcinoma investigating similarities in DNA
methylation measured by The Cancer Genome Atlas on the Illumina 450k platform.
At the core of the analysis is a reference free adjustment of cell type
proportion ([Houseman, Molitor, and Marsit. Bioinformatics
2014](10.1093/bioinformatics/btu029)) on each PAM50 subtype stratified by low
and high stages. After the adjustment, we observe key gene regions of
differential methylation in common to all PAM50 subtypes in the low stage. We
also validate these findings in a Validation set
([Yang et al. Genome Biology 2015](10.1186/s13059-015-0699-9)).

## Contact 

For all code related questions please file a [GitHub
issue](https://github.com/gwaygenomics/brca_lowstage_DMGRs/issues)

Questions regarding the analysis or other correspondance should be directed to:
Brock.C.Christensen@dartmouth.edu

## Analysis

All scripts are intended to be run in order, as defined in `run_pipeline.sh`.
The bash script can be run directly to reproduce the pipeline but this is not
recommended since some steps are computationally intensive. Instead, to
successfully reproduce this analysis, we recommend following the bash script
line by line.

The folder structure is labelled to facilitate easy step-wise navigation, as are
the `Scripts/` folders in each parent folder. All scripts should be run from the
top directory in the repository. For example: 

```sh
# Place downloaded IDAT files in this folder
IDAT_loc="../../Documents/mdata/TCGAbreast_idat/"

# Run preprocessing script 
R --no-save --args "I.Data_Processing/Data/TCGA_BRCA" $IDAT_loc "nosummary" < \
I.Data_Processing/Scripts/A.preprocess_minfi.R
```

## Data 

Data is not stored directly in this repository and should be downloaded
according to the `run_pipeline.sh` script.

## Acknowledgements 

This work was supported by the Institute for Quantitative Biomedical Sciences
and two grants P20GM104416 and R01DE02277 (BCC)

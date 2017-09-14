# Deconvolution of DNA methylation identifies differentially methylated gene regions on 1p36 across breast cancer subtypes

Titus, A.J., Way, G., Johnson, K., Christensen, B. Scientific Reports 2017 (DOI: doi:10.1038/s41598-017-10199-z)

[![DOI](https://zenodo.org/badge/45754471.svg)](https://zenodo.org/badge/latestdoi/45754471)


## Summary 

**Update: This work is [now published](https://www.nature.com/articles/s41598-017-10199-z) in Scientific Reports!**

The following repository contains all scripts required to reproduce an analysis
of early stage invasive breast carcinoma investigating similarities in DNA
methylation as measured by The Cancer Genome Atlas on the Illumina 450k platform.
At the core of the analysis is a reference-free adjustment of cell type
proportion ([Houseman _et al._ BMC Bioinformatics 2014](https://doi.org/10.1186/s12859-016-1140-4))
on each PAM50 subtype stratified by early and late stages. After the adjustment,
we observe key differentially methylated gene regions (DMGRs) in common to all
early stage tumors regardless of PAM50 subtype. We also validate these findings
in a Validation set ([Yang _et al._ Genome Biology 2015](10.1186/s13059-015-0699-9)).

Our findings implicate a small region localized entirely on chromosome 1p36.3
that harbors common DMGRs. This region has previously been shown to be important
for cancer initiation and prognosis.

## Contact 

For all code related questions please file a [GitHub
issue](https://github.com/Christensen-Lab-Dartmouth/brca_lowstage_DMGRs/issues)

Questions regarding the analysis or other correspondance should be directed to:
Brock.C.Christensen@dartmouth.edu

## Data 

Data is not stored directly in this repository and should be downloaded
according to the `run_pipeline.sh` script.

## Acknowledgements 

This work was supported by the Institute for Quantitative Biomedical Sciences
and two grants P20GM104416 and R01DE02277 (BCC) and a training grant T32LM012204 (AJT).

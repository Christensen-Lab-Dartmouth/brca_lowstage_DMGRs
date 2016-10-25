# Run preprocessing script (if third argument is "summary", then the script will output a processing summary)
# R --no-save --args "I.Data_Processing/Data/TCGA_BRCA" $IDAT_folder "nosummary" < I.Data_Processing/Scripts/A.preprocess_minfi.R

# Subset data to spectific subtypes; tumor adjacent samples are included in all subsets
# R --no-save --args PAM50.RNAseq Basal Her2 LumA LumB Normal < I.Data_Processing/Scripts/B.subsetBetas.R

#~~~~~EXPAND ANNOTATION FILE
# This file will expand the annotation file to include multiple cgs that are associated with multiple gene:regions
# R --no-save < I.Data_Processing/Scripts/C.editAnnotationFile.R # We will use this expanded annotation file in downstream analyses

# Describe the data you are using (Generate Table 1)
# To get information about CpG values, you need to have run a RefFreeEWAS already to understand which CpGs were actually used
# R --no-save < I.Data_Processing/Scripts/D.DescribeTableI.R

# Describe TCGA tumor purity estimates
R --no-save < I.Data_Processing/Scripts/E.TumorPurity.R

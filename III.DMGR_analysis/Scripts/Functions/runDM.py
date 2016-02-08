#####################################################################
# ~~~~~~~~~~~~~~~~~~
# Tumor subtype and cell type independent DNA methylation alterations 
# associated with stage progression in invasive breast carcinoma 
# ~~~~~~~~~~~~~~~~~~
# Way, G., Johnson, K., Christensen, B. 2015
#
# This python script will determine the differential CpGs in all models
# across each stage. It will run an R script for a prespecified q value
# cutoff. 
#####################################################################

import os

subtypes = ['Basal', 'Her2', 'LumA', 'LumB', 'Normal'] 
stages = ['low', 'high']

for i in subtypes:
	for j in stages:
		os.system('R --no-save --args '+j+' '+i+' 0.01 < III.DMGR_analysis/Scripts/A.DMcgs.R')
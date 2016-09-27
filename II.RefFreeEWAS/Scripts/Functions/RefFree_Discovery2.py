#####################################################################
# ~~~~~~~~~~~~~~~~~~
# Tumor subtype and cell type independent DNA methylation alterations 
# associated with stage progression in invasive breast carcinoma 
# ~~~~~~~~~~~~~~~~~~
# Way, G., Johnson, K., Christensen, B. 2015
#
# The custom python script will run the RefFree algorithm with 1000
# bootstraps. The script will take about two days to complete the largest
# data subset. It is possible to run on a personal machine but is not 
# recommended. The analysis was performed on the DISCOVERY cluster computer
# at Dartmouth College.
#####################################################################

import os
import sys
sys.path.insert(0, 'I.Data_Processing/Scripts/Functions')

subtypes = ['Basal', 'Her2', 'LumA', 'LumB', 'Normal']
stages = ['low', 'high']
# bootstraps = 1000
bootstraps = 1

for i in subtypes:
	for j in stages:
		os.system('R --no-save --args ' + bootstraps + ' ' + j + ' ' + i + ' < II.RefFreeEWAS/Scripts/A.doRefFree2.R')

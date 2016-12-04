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
bootstraps = 10

# for i in subtypes:
# 	for j in stages:
# 		os.system('R --no-save --args '+str(bootstraps)+' '+j+' '+i+' < II.RefFreeEWAS/Scripts/A.doRefFree2.R')

#os.system('R --no-save --args '+str(bootstraps)+' low Basal ' + '< II.RefFreeEWAS/Scripts/A.doRefFree2.R')
#os.system('R --no-save --args '+str(bootstraps)+' high Basal ' + '< II.RefFreeEWAS/Scripts/A.doRefFree2.R')
#os.system('R --no-save --args '+str(bootstraps)+' low Her2 ' + '< II.RefFreeEWAS/Scripts/A.doRefFree2.R')
#os.system('R --no-save --args '+str(bootstraps)+' high Her2 ' + '< II.RefFreeEWAS/Scripts/A.doRefFree2.R')
os.system('R --no-save --args '+str(bootstraps)+' low LumA ' + '< II.RefFreeEWAS/Scripts/A.doRefFree2.R')
#os.system('R --no-save --args '+str(bootstraps)+' high LumA ' + '< II.RefFreeEWAS/Scripts/A.doRefFree2.R')
#os.system('R --no-save --args '+str(bootstraps)+' low LumB ' + '< II.RefFreeEWAS/Scripts/A.doRefFree2.R')
#os.system('R --no-save --args '+str(bootstraps)+' high LumB ' + '< II.RefFreeEWAS/Scripts/A.doRefFree2.R')
#os.system('R --no-save --args '+str(bootstraps)+' low Normal ' + '< II.RefFreeEWAS/Scripts/A.doRefFree2.R')
#os.system('R --no-save --args '+str(bootstraps)+' high Normal ' + '< II.RefFreeEWAS/Scripts/A.doRefFree2.R')

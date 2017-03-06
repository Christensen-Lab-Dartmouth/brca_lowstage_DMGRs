#####################################################################
# ~~~~~~~~~~~~~~~~~~
# Tumor subtype and cell type independent DNA methylation alterations 
# associated with stage progression in invasive breast carcinoma 
# ~~~~~~~~~~~~~~~~~~
# Way, G., Johnson, K., Christensen, B. 2015
#
# This script will store all of the custom functions that are used in 
# the doRefFree analyses
#####################################################################

################################
# SubsetStage Function
################################
# This function will subset stage for the input dataframe of covariates
subsetStage <- function (data, stage) {
        tmpframe <- c()
        
        # Separate the data according to the given stage assignments
        for (i in 1:length(stage)) {
                tmp <- data[data$pathologic_stage == stage[i], ]
                # Accept only primary tumors
                tmp <- tmp[tmp$sample.type == "Primary Tumor", ]
                # Combine the info together
                tmpframe <- rbind(tmpframe, tmp)
        }
        
        # Also, take all the solid tissue normal samples in this step as well
        tmpframe <- rbind(tmpframe, data[data$sample.type == "Solid Tissue Normal", ])
        return(tmpframe)
}

################################
# customRefFree Function
################################
# This function will make the design matrices, run RefFreeEWAS, and save output files
# k is the number of dimensions; it can be prespecified if you do not wish to run EstDimRMT
customRefFree <- function (covariates, betas, age = T, findRMT = T, bootstraps = 5, 
                           save = F, iPOI = 2, k) {
        
        # subset betas
        newBeta <- betas[ , rownames(covariates)]
        
        # Step 0: You can start with an intial estimate of M (Optional)
        # Komen_Initialized = RefFreeCellMixInitialize(Betas_Komen, K=5)
        
        # Step 1 - 2: Alternate fixing Mu and Omega by iterating from 2 to 10 (Kmax) cell types
        DMGR_RefFree_Array <- RefFreeCellMixArray(newBeta, Klist=3:10, iters=25)
        
        # Step 3: Bootstrap method for determining the optimal number of Classes K
        RefFree_DMGR_Boots = RefFreeCellMixArrayDevianceBoots(DMGR_RefFree_Array, newBeta, R=1000, bootstrapIterations=bootstraps)
        
        # Get the list ready to return - Possibly remove if the following code runs
        returnlist <- list(DMGR_RefFree_Array, RefFree_DMGR_Boots)
        return(returnlist)
        
        
}


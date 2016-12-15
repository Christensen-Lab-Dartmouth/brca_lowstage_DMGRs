#################
#Code taken from:
#https://github.com/kasperdanielhansen/minfi/blob/master/R/read.450k.R
#################
#There was an error in the normal read.450k.exp function in the minfi package
#This function overcomes the error

read.450k.exp2 <- function(base = NULL, targets = NULL, extended = FALSE, 
                           recursive = FALSE, verbose = FALSE) {
  if(!is.null(targets)) {
    if(! "Basename" %in% names(targets))
      stop("Need 'Basename' amongst the column names of 'targets'")
    if(!is.null(base)) {
      files <- file.path(base, targets$Basename)
    } else {
      files <- targets$Basename
    }
    rgSet <- read.metharray(files, extended = extended, verbose = verbose)
    pD <- targets
    pD$filenames <- files
    rownames(pD) <- sampleNames(rgSet)
    pData(rgSet) <- pD
    return(rgSet)
  }
  ## Now we just read all files in the directory
  Grn.files <- list.files(base, pattern = "_Grn.idat$", recursive = recursive,
                          ignore.case = TRUE, full.names = TRUE)
  Red.files <- list.files(base, pattern = "_Red.idat$", recursive = recursive,
                          ignore.case = TRUE, full.names = TRUE)
  if(length(Grn.files) == 0 || length(Red.files) == 0)
    stop("No IDAT files were found")
  commonFiles <- intersect(sub("_Grn.idat$", "", Grn.files), sub("_Red.idat$", "", Red.files))
  if(length(commonFiles) == 0)
    stop("No IDAT files with both Red and Green channel were found")
  commonFiles.Grn <- paste(commonFiles, "_Grn.idat", sep = "")
  commonFiles.Red <- paste(commonFiles, "_Red.idat", sep = "")
  if(!setequal(commonFiles.Grn, Grn.files))
    warning(sprintf("the following files only exists for the green channel: %s",
                    paste(setdiff(Grn.files, commonFiles.Grn), collapse = ", ")))
  if(!setequal(commonFiles.Red, Red.files))
    warning(sprintf("the following files only exists for the red channel: %s",
                    paste(setdiff(Red.files, commonFiles.Red), collapse = ", ")))
  rgSet <- read.450k(basenames = commonFiles, extended = extended, verbose = verbose)
  rgSet
}
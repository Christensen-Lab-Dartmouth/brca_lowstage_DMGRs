# Gregory Way 2016
#
# INSTALL.R
#
# USAGE:
# Run to install all required R libraries: "Rscript INSTALL.R"

library("methods")

mirror <- 'http://cran.us.r-project.org'
install.packages("checkpoint", repos = mirror)

library("checkpoint")

dir.create(".checkpoint")
checkpoint("2016-09-01", checkpointLocation = ".")

# CRAN
cran <- c("plyr", "reshape2", "ggplot2", "qvalue", "gridExtra", "readr", "Hmisc",
          "VennDiagram")
install.packages(cran, repos='http://cran.us.r-project.org')

# Bioconductor
source("http://bioconductor.org/biocLite.R")

bioC <- c("minfi",
	  "IlluminaHumanMethylation450kmanifest",
	  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
          "RefFreeEWAS",
	  "isva",
	  "limma",
	  "Homo.sapiens",
	  "Gviz",
	  "GEOquery",
	  "qvalue")

biocLite(bioC, suppressUpdates = TRUE)

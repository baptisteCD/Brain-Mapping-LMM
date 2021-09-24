# Wrapper to extract the values of interest from the simulated BWAS results

setwd("path/to/working/directory")
source("RR_3.1_BWAS_manhathanPlot_functions.R")

library(pROC)
library(qqman)

# read arguments provided when submitting batch job
arg = commandArgs(trailingOnly=TRUE)

# Arguments are
# 1) pathToTheBWASresults
# 2) iter
# 3) NNN
# 4) output folder

summ=summaryBWAS_allModa_NULL(path = arg[1], iter=as.numeric(arg[2]), NNN = arg[3])
write.table(summ, paste0(arg[4], "/Summary_BWAS_NVertices",arg[3], "_Niter",arg[2]), col.names = T, row.names = F )

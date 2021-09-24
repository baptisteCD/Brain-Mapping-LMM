# Wrapper to extract the values of interest from the simulated BWAS results

setwd("path/to/working/directory")
source("RR_3.1_BWAS_manhathanPlot_functions.R")

library(Morpho)
library(plyr)
library(viridis)

# read arguments provided when submitting batch job
arg = commandArgs(trailingOnly=TRUE)

# Arguments are
# 1) pathToTheBWASresults/scenario
# 2) ntotIter
# 3) NNN
# 4) output folder
# 5) iter
#
signifTList=c(7.665385e-08, 5e-8, 1e-8,  5e-9, 1e-9, 5e-10, 1e-10, 5e-11, 1e-11, 5e-12, 1e-12, 5e-13, 1e-13, 5e-14, 1e-14, 5e-15, 1e-15, 5e-16, 1e-16, 5e-17, 1e-17, 5e-18, 1e-18, 5e-19, 1e-19, 5e-20, 1e-20, 5e-21, 1e-21, 5e-22, 1e-22, 5e-23, 1e-23, 1e-24, 1e-25, 1e-26,1e-27, 1e-28, 1e-29, 1e-30, 1e-31, 1e-32, 1e-33, 1e-34, 1e-35, 1e-36, 1e-37, 1e-38, 1e-39,  1e-40, 1e-42, 1e-44, 1e-46, 1e-48, 1e-50, 1e-52)

# Open list of associated vertices
tpList=read.table(paste0("01_vertexSelection/NVertices", arg[3], "_Niter", arg[5] ), stringsAsFactors = F)

if (file.exists(paste0(arg[1], "/BWAS_NVertices", arg[3], "_Niter", arg[5], ".linear"))){
bwasResult=read.table(paste0(arg[1], "/BWAS_NVertices", arg[3], "_Niter", arg[5], ".linear"), header=T, stringsAsFactors = F) } else if (file.exists(paste0(arg[1], "/BWAS_NVertices", arg[3], "_Niter", arg[5], ".moa"))){
bwasResult=read.table(paste0(arg[1], "/BWAS_NVertices", arg[3], "_Niter", arg[5], ".moa"), header=T, stringsAsFactors = F)  } else {
bwasResult=read.table(paste0(arg[1], "/BWAS_NVertices", arg[3], "_Niter", arg[5], ".moment"), header=T, stringsAsFactors = F)
}

resTotTable=data.frame(c(signifTList))
resTot=NULL


for (moda in c("area", "thickness", "thick", "LogJacs")){

# Open and format bwas summary statistics
if (moda %in% c("thick", "LogJacs")){
bwasFormatted=formatBWASsubcortical(BWASsumstat = bwasResult)
bwasFormatted=bwasFormatted[grep(bwasFormatted$Probe, pattern = moda),]
bwas=bwasAddCoordinatesSubcortical(bwasFormatted = bwasFormatted)
}

if (moda %in% c("thickness", "area")){
# Format left and right hemisphere data
bwasFormattedlh=formatBWAScortical(BWASsumstat = bwasResult, hemi = "lh", mod = moda)
bwaslh=bwasAddCoordinatesCortical(bwasFormatted = bwasFormattedlh , hemi = "lh", moda = moda)
bwasFormattedrh=formatBWAScortical(BWASsumstat = bwasResult, hemi = "rh", mod = moda)
bwasrh=bwasAddCoordinatesCortical(bwasFormatted = bwasFormattedrh , hemi = "rh", moda = moda)
# Rbind left and right
bwas=rbind(bwasrh, bwaslh)
}

for (signifThreshold in signifTList){

# Run cluster statistic
res=clusterFWER(path=arg[4], bwas = bwas, tpList = tpList, NNN = arg[3], moda = moda, signifThreshold = signifThreshold, KNearestNeighbours = 10, iter=arg[5])
colnames(res)=paste0(colnames(res), "_", moda)
resTot=rbind(resTot, res)
}

resTotTable=cbind(resTotTable, resTot)
resTot=NULL
}

# Write final output - with all iterations
write.table(resTotTable, paste0(arg[4], "/Results_clusterFWER_NVertices", arg[3],"_iter",arg[5],  "_allSignifT_extra2.txt" ))




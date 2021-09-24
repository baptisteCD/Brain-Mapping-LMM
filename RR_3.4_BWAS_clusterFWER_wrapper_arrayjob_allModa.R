# Wrapper to extract the values of interest from the simulated BWAS results

setwd("/scratch/project/genetic_data_analysis/uqbcouvy/56_BWAS")
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

signifThreshold=0.05/652283

# Open list of associated vertices
tpList=read.table(paste0("01_vertexSelection/NVertices", arg[3], "_Niter", arg[5] ), stringsAsFactors = F)

if (file.exists(paste0(arg[1], "/BWAS_NVertices", arg[3], "_Niter", arg[5], ".linear"))){
bwasResult=read.table(paste0(arg[1], "/BWAS_NVertices", arg[3], "_Niter", arg[5], ".linear"), header=T, stringsAsFactors = F) } else if (file.exists(paste0(arg[1], "/BWAS_NVertices", arg[3], "_Niter", arg[5], ".moa"))){
bwasResult=read.table(paste0(arg[1], "/BWAS_NVertices", arg[3], "_Niter", arg[5], ".moa"), header=T, stringsAsFactors = F)  } else {
bwasResult=read.table(paste0(arg[1], "/BWAS_NVertices", arg[3], "_Niter", arg[5], ".moment"), header=T, stringsAsFactors = F)
}
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

# Run cluster statistic
res=clusterFWER(path=arg[4], bwas = bwas, tpList = tpList, NNN = arg[3], moda = moda, signifThreshold = signifThreshold, KNearestNeighbours = 10, iter=arg[5])
colnames(res)=paste0(colnames(res), "_", moda)
resTot=c(resTot, res)
}

# Write final output - with all iterations
write.table(resTot, paste0(arg[4], "/Results_clusterFWER_NVertices", arg[3],"_iter",arg[5],  ".txt" ))





atlasPath="path/to/atlas/folder"

#####################
# Functions
formatBWAScortical<-function(BWASsumstat, hemi, mod){
BWASsumstat$ProbeID<-BWASsumstat$Probe

atlas<-read.table(paste0(atlasPath, "fsaverage_vertex_labels_names_order2.txt"), header=T, stringsAsFactors = F)
# adapt probe prefix to modality and hemisphere
if (hemi=="lh" & mod=="thickness"){ prefix="lht"}
if (hemi=="lh" & mod=="area"){ prefix="lha" }
if (hemi=="rh" & mod=="thickness"){ prefix="rht" }
if (hemi=="rh" & mod=="area"){ prefix="rha"}
atlas$ProbeID<-paste0(prefix, "_",substr(atlas$VertexNum, 11,nchar(atlas$VertexNum)))
BWASsumstat<-merge(BWASsumstat, atlas, by="ProbeID" )

if (hemi=="lh"){
  if ( length(which(is.na(BWASsumstat$LHLabel)) )>0){
  BWASsumstat<-BWASsumstat[-which(is.na(BWASsumstat$LHLabel)),] }
BWASsumstat<-BWASsumstat[order(BWASsumstat[,"newOrderLeft"]),]
BWASsumstat$LABEL<-BWASsumstat[,"LeftHem_Label"]
BWASsumstat$LABELNAME<-BWASsumstat[,"LHLabel"]
BWASsumstat$rgbx<-BWASsumstat[,"Lrgbx"]
BWASsumstat$rgby<-BWASsumstat[,"Lrgby"]
BWASsumstat$rgbz<-BWASsumstat[,"Lrgbz"]
}
if (hemi=="rh"){
  if ( length(which(is.na(BWASsumstat$RHLabel)) )>0){
  BWASsumstat<-BWASsumstat[-which(is.na(BWASsumstat$RHLabel)),] }
BWASsumstat<-BWASsumstat[order(BWASsumstat[,"newOrderRight"]),]
BWASsumstat$LABEL<-BWASsumstat[,"RightHem_Label"]
BWASsumstat$LABELNAME<-BWASsumstat[,"RHLabel"]
BWASsumstat$rgbx<-BWASsumstat[,"Rrgbx"]
BWASsumstat$rgby<-BWASsumstat[,"Rrgby"]
BWASsumstat$rgbz<-BWASsumstat[,"Rrgbz"]
}
BWASsumstat$xax<-1:length(BWASsumstat[,1])

return(BWASsumstat)
}


#############
formatBWASsubcortical<-function(BWASsumstat){

BWASsumstat$ProbeID<-gsub(".*[a-z][_]", "", BWASsumstat$Probe)


atlas<-read.table(paste0(atlasPath, "fsaverage_subcortical_vertices_label.txt"), header=T, stringsAsFactors = F)
BWASsumstat<-merge(BWASsumstat, atlas, by="ProbeID" )

BWASsumstat<-BWASsumstat[order(BWASsumstat$subcv),]
BWASsumstat$LABEL<-BWASsumstat[,"subcv"]
BWASsumstat$LABELNAME<-BWASsumstat[,"subcvName"]

BWASsumstat$xax<-1:length(BWASsumstat[,1])

return(BWASsumstat)
}

######################
#####################
# BWAS plot

BWASPlotSimple<-function(path , bwasFile , yMax, ntotvar, phenotypeLabel){
library(png)
library(plyr)
library(qqman)

BWASsignif=NULL
chi2All=NULL
pAll=NULL
lambda=NULL
signifLog10= -log( (0.05/652283)/ntotvar, base=10 )

# Initialise plot
png(paste0(path, "/Manhathan_", bwasFile, "_simple.png"), width = 50, height = 16, units = "cm", res = 400, type="cairo")
plot.new()
layout(matrix(c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,6,7,8), nrow =1 , ncol = 20, byrow = F))

BWASsumstatfilePath = paste0(path ,  bwasFile)
BWASsumstatfile<-  read.table(BWASsumstatfilePath, header=T, stringsAsFactors = F)
BWASsumstatfile$log10p<- -log10(BWASsumstatfile$p)

for (mod in c( "thickness","area", "thick", "LogJacs")){
  for (hemi in c("lh", "rh")){
BWASsumstat=NULL

  if (mod %in% c("area", "thickness")){
BWASsumstatHemi<-formatBWAScortical(BWASsumstat = BWASsumstatfile,hemi = hemi, mod=mod)
BWASsumstat=rbind(BWASsumstat,BWASsumstatHemi)
  }

   if (mod %in% c("LogJacs", "thick")){
BWASsumstat<-formatBWASsubcortical(BWASsumstat = BWASsumstatfile)
BWASsumstat=BWASsumstat[grep(BWASsumstat$Probe, pattern = mod),]
if (hemi == "lh"){
  BWASsumstat=BWASsumstat[which(BWASsumstat$subcv<30),]
} else{
  BWASsumstat=BWASsumstat[which(BWASsumstat$subcv>30),]
  BWASsumstat$xax=1:length(BWASsumstat$ProbeID)
}
}

# Get all chi2 and pvalues
BWASsumstat$CHI2=(BWASsumstat$b/BWASsumstat$se)**2

# Set color
laab<-unique(BWASsumstat$LABEL)
colvect<-rep(c("#56B4E9","#0072B2"),length(laab)+1)
BWASsumstat$col<-0
jjj<-1
for (iii in laab){
  BWASsumstat$col[which(BWASsumstat$LABEL==iii)]<-colvect[jjj]
  jjj<-jjj+1
}

# Get rid of unused vertices
if(length(which(is.na(BWASsumstat$log10p)))>0){
BWASsumstat<-BWASsumstat[-which(is.na(BWASsumstat$log10p)),]
}

# Add hemi-moda to text file
BWASsumstat$hemi=hemi
BWASsumstat$mod=mod

# Store all vertices
BWASsignif=rbind.fill(BWASsignif, BWASsumstat)

# Manhathan plot of association
if (hemi=="lh" & mod=="thickness"){
  par(mar=c(8,4,2,0.3))
plot(BWASsumstat$xax, BWASsumstat$log10p, col=BWASsumstat$col, pch=20, cex=1,  xlab="", xaxt="n", cex.lab=1, ylim=c(0, yMax+1), ylab="-log10(p-value)", bty="n", main="", cex.main=0.7, xaxs="i", yaxs="i", xlim=c(-10, max(BWASsumstat$xax)+10))
title("Cortical thickness \n left", line = -1, cex.main=0.7)
title(paste0(phenotypeLabel), line = 1, cex.main=0.9, adj=0)
} else{

if(mod=="area"){ moda<- "Cortical area"} else if (mod=="thickness"){ moda<- "Cortical thickness" } else if (mod=="LogJacs"){ moda<- "Subcortical area" } else if (mod=="thick"){ moda<- "Subcortical thickness" }

titl<-paste( moda , "\n",ifelse(hemi=="lh", "Left", "Right"))
par(mar=c(8,0.5,2,0.5))
plot(BWASsumstat$xax, BWASsumstat$log10p, col=BWASsumstat$col, pch=20, cex=1,  xlab="", xaxt="n",  yaxt="n", cex.lab=0.9, ylim=c(0, yMax+1), ylab="", bty="n", main="", cex.main=0.7, xaxs="i", yaxs="i", xlim=c(-10, max(BWASsumstat$xax)+10) )
title(titl, line = -1, cex.main=0.7)
}
# add ticks and labels
# Change axis
laab<-unique(BWASsumstat$LABEL)
idchange<-NULL
for (iii in laab){
 idchange<-c(idchange, BWASsumstat$xax[which(BWASsumstat$LABEL==iii)][1])
}
idchange<-c(idchange, max(BWASsumstat$xax))
tickNb<-idchange[1:(length(idchange)-1)]+(idchange[2:length(idchange)]-idchange[1:(length(idchange)-1)])/2
tickNb[-1]=tickNb[-1]-0.002*max(BWASsumstat$xax)

odd <- function(x) x%%2 != 0
axis(1, at = tickNb[odd(which(tickNb>0))], las=2, labels = unique(BWASsumstat$LABELNAME)[odd(which(tickNb>0))], lwd.ticks = 1.5, lwd = 1.5 ,cex.axis=0.8, col.axis = c("#56B4E9"), xaxs="i")
axis(1, at = tickNb[!odd(which(tickNb>0))], las=2, labels = unique(BWASsumstat$LABELNAME)[!odd(which(tickNb>0))], lwd.ticks = 1.5, lwd = 1.5 ,cex.axis=0.8, col.axis = c("#0072B2"), xaxs="i")

abline(h= signifLog10 , col="darkred", lwd=2)
    }
}
dev.off()

write.csv(BWASsignif, paste0(path, "/BWAS_fullSummary_", phenotypeLabel, ".csv" ))
write.csv(BWASsignif[which(BWASsignif$log10p > signifLog10),], paste0(path, "/BWAS_signif_", phenotypeLabel , ".csv" ))

}

####################################
# QQPlot

qqplot <- function(pvector, col = "black", add = F, ylim = NULL){
expectedP <- -log10(ppoints(length(pvector)))
  observedP <- -log10(sort(pvector, decreasing = F))
  if (add == F) {
    plot(x = expectedP, y = observedP, col = col,
         xlab = expression(Expected ~ ~-log[10](italic(p))),
         ylab = expression(Observed ~ ~-log[10](italic(p))),
         pch = 16, cex = 0.7, ylim = NULL)
    abline(0, 1, col = "red")
  }else{
    points(x = expectedP, y = observedP, col = col,
           xlab = expression(Expected ~ ~-log[10](italic(p))),
           ylab = expression(Observed ~ ~-log[10](italic(p))),
           pch = 16, cex = 0.7, ylim = NULL)
  }
}


# QQplot that superimpose results from different models
superimposedQQplot=function(path , scenarioList,  legendList, colourList, variableLabel){
library(readr)
  library(qqman)
jjj=2
png(paste0(path, "QQplotCombined_assocVertices", variableLabel, ".png"), width = 15, height = 15, units = "cm", res=400)
par(mar=c(4,4,2,1))
bwas=read_csv(paste0(path, scenarioList[1], "/BWAS_fullSummary_", variableLabel, ".csv"))
qqplot(bwas$p, col=colourList[1])
 print(round(median(bwas$CHI2,na.rm=T)/qchisq(0.5,df=1),3))
legend(x = 3.5, y = 0.75*max(bwas$log10p), legend = legendList ,pch=20, pt.cex=1.5,  col = colourList)

 for (scenario in scenarioList[-1] ){
bwas=read_csv(paste0(path, scenarioList[jjj], "/BWAS_fullSummary_", variableLabel, ".csv"))
qqplot(bwas$p, col=colourList[jjj], add=T)
 print(print(round(median(bwas$CHI2,na.rm=T)/qchisq(0.5,df=1),3)))
jjj=jjj+1
 }

dev.off()
}

##################################
# Add coordinates to formated BWAS file

bwasAddCoordinatesCortical=function(bwasFormatted, hemi, moda){
# Get vertices coordinates
coor=read.table(paste0(atlasPath, hemi, ".vertexCoordinates.asc"))
coor$V1=paste0(hemi, substr(moda, 1,1), "_", coor$V1)
colnames(coor)=c("ProbeID", "X", "Y", "Z", "V4")
# Add coordinates
bwasFormatted=merge(bwasFormatted, coor[,1:4], by="ProbeID")
bwasFormatted$X=bwasFormatted$X*(-1) # coordinates are left/right swapped in asc format
bwasFormatted$hemi=hemi
return(bwasFormatted)
}

################################
# Function to identify cluster
significantClusterAroundVertex=function(bwas, signifThreshold, truePositiveLabel, nbNearestNeighbours ){
library(Rvcg)

  if (bwas$p[which(bwas$Probe == truePositiveLabel)]> signifThreshold){
  print(paste("Vertex", truePositiveLabel, "located in", bwas$LABELNAME[which(bwas$Probe==truePositiveLabel)], "not significant" ))
  } else {

bwas$voxelStatus="nonSignificant"
bwas$voxelStatus[which(bwas$p < signifThreshold)]="bwSignificant"
bwas$voxelStatus[which(bwas$Probe == truePositiveLabel)]="truePositive"

# First iter - closest vertices
clost= vcgKDtree( target = as.matrix(bwas[,c("X", "Y", "Z")]) ,  query = as.matrix(bwas[which(bwas$Probe==truePositiveLabel),c("X", "Y", "Z")]), k=nbNearestNeighbours)$index # Get vertices number

# Mark the vertices and exclude NS ones and TP one
bwas$voxelStatus[clost]="belongInCluster"
bwas$voxelStatus[which(bwas$p > signifThreshold)]="nonSignificant"
bwas$voxelStatus[which(bwas$Probe==truePositiveLabel)]="truePositive"
nCluster=length(which(bwas$voxelStatus=="belongInCluster"))

# Loop until all cluster has been identified
nIncrementCluster=nCluster
sizeCluster=0

while(nIncrementCluster>0){
clost=NULL
clost <- vcgKDtree( target = as.matrix(bwas[,c("X", "Y", "Z")]) ,  query = as.matrix(bwas[which(bwas$voxelStatus %in% c("belongInCluster", "truePositive")),c( "X","Y", "Z")]), k=nbNearestNeighbours)$index
bwas$voxelStatus[unique(c(clost))]="belongInCluster"
bwas$voxelStatus[which(bwas$p > signifThreshold)]="nonSignificant"
nIncrementCluster=length(which(bwas$voxelStatus=="belongInCluster"))-sizeCluster
sizeCluster=length(which(bwas$voxelStatus=="belongInCluster"))
}

bwas[,paste0("cluster_", truePositiveLabel)]=ifelse(bwas$voxelStatus %in% c("belongInCluster", "truePositive"), 1 , 0 )
print(paste("Found ", sum(bwas[,paste0("cluster_", truePositiveLabel)]),"significant vertices in cluster including ", truePositiveLabel, "located in", bwas$LABELNAME[which(bwas$Probe==truePositiveLabel)] ))
bwas$voxelStatus=NULL

}

return(bwas)
}


############################################
# Extract cluster FWER metrics results
clusterIdentification=function(path, bwas, moda, signifThreshold, KNearestNeighbours, phenotypeLabel){

library(plyr)
bwas11=NULL
twoTPinCLuster=fwerProblemIncreaseK=0
bwas$signifVoxel=ifelse(bwas$p < signifThreshold, 1 , 0)

# Depending on modality - loop on hemisphere or subcortical structures
if (moda %in% c("area", "thickness")){ ROIlist=c("lh", "rh")}
if (moda %in% c("LogJacs", "thick")){ ROIlist=c(10, 11, 12, 13, 17, 18, 26, 49, 50, 51, 52, 53, 54, 58)}

# Loop on the different structures and concatenate results
  for (ROI in ROIlist){
  print(ROI)
  if (moda %in% c("LogJacs", "thick")){ bwas1=bwas[which(bwas$LABEL==ROI),] }
  if (moda %in% c("area", "thickness")){ bwas1=bwas[which(bwas$hemi==ROI),] }

  # Find significant clusters
  bwas1$inCluster=0
  while (length(which(bwas1$signifVoxel==1 & bwas1$inCluster==0))>0){
  signifVertex=bwas1$Probe[which(bwas1$signifVoxel==1  & bwas1$inCluster==0)][1]
  bwas1=significantClusterAroundVertex(bwas=bwas1, signifThreshold=signifThreshold, truePositiveLabel = signifVertex, nbNearestNeighbours = KNearestNeighbours)
  bwas1$inCluster=rowSums(bwas1[,grep(x = colnames(bwas1), pattern = "cluster_"), drop=FALSE])

    if (length(table(bwas1$inCluster))>2){
 print("WARNING: some vertices have been attributed to several clusters. You may want to check the results or you may want to modify the nearest neighbour search option k, in mcNNindex() (used in function significantClusterAroundVertex), in the meantime,  the 2 clusters are merged")

    # Merge clusters and resets variable
    TPvertex=bwas1$Probe[which(bwas1$inCluster==2 & bwas1$signifVoxel==1)][1]
    bwas1[which(bwas1[,paste0("cluster_", signifVertex)]==1),paste0("cluster_", TPvertex)]=1
    bwas1[,paste0("cluster_", signifVertex)]=NULL
    fwerProblemIncreaseK=1
    bwas1$inCluster=rowSums(bwas1[,grep(x = colnames(bwas1), pattern = "cluster_"), drop=FALSE])
    }
  }
    # Combine all ROI - hemisphere of subcortical volumes
    bwas11=rbind.fill(bwas11, bwas1)
  }
    # Extract statistics and information on the full data
    nbClusters=length(grep(x = colnames(bwas11), pattern = "cluster_"))
    maxSizeClusters=max(colSums(bwas11[,grep(x = colnames(bwas11), pattern = "cluster_"), drop=FALSE], na.rm = T), na.rm = T)
    medianSizeClusters=median(colSums(bwas11[,grep(x = colnames(bwas11), pattern = "cluster_"), drop=FALSE], na.rm = T), na.rm = T)
    minSizeClusters=min(colSums(bwas11[,grep(x = colnames(bwas11), pattern = "cluster_"), drop=FALSE], na.rm = T), na.rm = T)

    FWERcluster_ROI=0
    for (iii in grep(x = colnames(bwas11), pattern = "cluster_")){
      if( length(table(bwas11$LABELNAME[which(bwas11[,iii]==1)]))>1 ){ FWERcluster_ROI=1   }   }

     # Format results
    res=as.data.frame(t(c( nbClusters, maxSizeClusters,medianSizeClusters, minSizeClusters, FWERcluster_ROI)))
    colnames(res)=c("NumberClusters", "maxSizeCluster", "medianSizeCluster", "minSizeCluster", "FWERcluster_ROI")
    # Write outputs and bwas table with all the informations (for plots etc..)

    write.table(bwas11, file = paste0(path, "/BWAS_NVertices",phenotypeLabel , "_", moda, "_clustersAndCoordinates"), quote = F )
    return(res)
}


######################
# ADD subcortical coordinates

bwasAddCoordinatesSubcortical=function(bwasFormatted, hemi, moda){
# Get vertices coordinates
coor=read.table(paste0(atlasPath, "Atlas_coordinates_subortical.txt"), stringsAsFactors = F, header=T)
# Add coordinates
bwasFormatted=merge(bwasFormatted, coor, by="ProbeID")
return(bwasFormatted)
}



########################
# Plot subcortical

plotSubcortical=function(path, variable){
library(rgl)
  for (moda in c("thick", "LogJacs")){
 for (hemi in c("lh", "rh")){
bwasPlot=read.table(paste0(path,  "/BWAS_NVertices",variable , "_", moda, "_clustersAndCoordinates"), header=T)
bwasPlot=bwasPlot[which(bwasPlot$hemi==hemi),]
bwasPlot$X=bwasPlot$X*(-1)
bwasPlot$Y=bwasPlot$Y*(-1)

# Open phenotype file
pheno=read.table(paste0("UKB_phenotypes15K/", variable, ".phen"))

# Transform betas in correlations
bwasPlot$cor=bwasPlot$b/sqrt(var(pheno$V3, na.rm = T))

# Atttribute colors on a diverging palette
bwasPlot$colorScale <- cut(bwasPlot$cor, breaks = seq(-0.1, 0.1, len = 21),  include.lowest = TRUE)

## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
library(RColorBrewer)
cols=c(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[10:6],RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[5:1]) # Select palette colours
bwasPlot$color <- colorRampPalette(c(cols))(21)[bwasPlot$colorScale] # Make it more continuous
bwasPlot$color[which(bwasPlot$signifVoxel==0)]="darkgrey"
bwasPlot$radius=ifelse(  bwasPlot$signifVoxel==1, 1.5 ,0.8 )

# Draw plots and save screenshots
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot[,c( "Z",  "X", "Y")]), col=bwasPlot$color, radius = bwasPlot$radius)
rgl.snapshot(paste0(path, "/BWAS_", variable, "_", hemi, "_", moda , "_clustersAndCoordinates_inside.png"))
rgl.close()

bwasPlot$X=bwasPlot$X*(-1)
bwasPlot$Z=bwasPlot$Z*(-1)
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot[,c( "Z",  "X", "Y")]), col=bwasPlot$color, radius = bwasPlot$radius)
rgl.snapshot(paste0(path, "/BWAS_", variable,"_", hemi, "_", moda , "_clustersAndCoordinates_outside.png"))
rgl.close()

}}
}

############################
# Plot cortical
plotCortical=function(path, variable){
library(rgl)
 for (moda in c("area", "thickness")){

 for (hemi in c("rh", "lh")){
bwasPlot=read.table(paste0(path,  "/BWAS_NVertices",variable , "_", moda, "_clustersAndCoordinates"), header=T)
bwasPlot=bwasPlot[which(bwasPlot$hemi==hemi),]

# Open phenotype file
pheno=read.table(paste0("UKB_phenotypes15K/", variable, ".phen"))

# Transform betas in correlations
bwasPlot$cor=bwasPlot$b/sqrt(var(pheno$V3, na.rm = T))

# Atttribute colors on a diverging palette
bwasPlot$colorScale <- cut(bwasPlot$cor, breaks = seq(-0.4, 0.4, len = 21),  include.lowest = TRUE)

## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
library(RColorBrewer)
cols=c(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[10:6], RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[5:1]) # Select palette colours
bwasPlot$color <- colorRampPalette(c(cols))(21)[bwasPlot$colorScale] # Make it more continuous
bwasPlot$color[which(bwasPlot$signifVoxel==0)]="darkgrey"
bwasPlot$radius=ifelse(bwasPlot$signifVoxel==1, 2, 0.8)

#bwasPlot2=bwasPlot[sample(1:length(bwasPlot$ProbeID), size = 10000),]
#bwasPlot2=rbind(bwasPlot2, bwasPlot[which(bwasPlot$signifVoxel==1 | bwasPlot$inCluster>0 ),])
bwasPlot2=bwasPlot


# Draw plot
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot2[,c( "Y","X", "Z")]), col=bwasPlot2$color, radius = bwasPlot2$radius)
rgl.snapshot(paste0(path, "/BWAS_", variable,"_", hemi, "_", moda , "_clustersAndCoordinates_inside.png"))
rgl.close()

bwasPlot2$X=bwasPlot2$X*(-1)
bwasPlot2$Y=bwasPlot2$Y*(-1)
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot2[,c( "Y","X", "Z")]), col=bwasPlot2$color, radius = bwasPlot2$radius)
rgl.snapshot(paste0(path, "/BWAS_", variable,"_", hemi, "_", moda , "_clustersAndCoordinates_outside.png"))
rgl.close()

}}
}


##############################
# Combine cortical and subcortical plots
combineCorticalSubcorticalPlots=function(path, variable){
library(grid)
library(png)
library(ggplot2)
library(gridExtra)

# List of files
ll=NULL
for (moda in c("thickness", "area","thick", "LogJacs") ){
ll=c(ll, paste0(path, "/BWAS_", variable, "_", "lh", "_", moda, "_clustersAndCoordinates_outside.png"))
ll=c(ll, paste0(path, "/BWAS_", variable, "_","lh", "_", moda, "_clustersAndCoordinates_inside.png"))
ll=c(ll, paste0(path, "/BWAS_", variable, "_", "rh", "_", moda, "_clustersAndCoordinates_inside.png"))
ll=c(ll, paste0(path, "/BWAS_", variable, "_", "rh", "_", moda, "_clustersAndCoordinates_outside.png"))
}
ll=c(ll, paste0(path, "/../legendbar.png"))

# Crop and convert format
plots2 <- lapply(ll<-ll,function(x){
  if(x!=paste0(path, "/../legendbar.png")){
   img <- as.raster(readPNG(x)[,100:1100,])} else {
    # img <- as.raster(readPNG(x))} else {
     img <- as.raster(readPNG(x)[,,])
    }
    rasterGrob(img, interpolate = T)

})


# Lay and write png
lay <- rbind(c(1,1,3,3,5,5,7,7,9,9,11,11,13,13,15,15,17),
            c(2,2,4,4,6,6,8,8,10,10,12,12,14,14,16,16, 17))
# lay <- rbind(c(1,1,3,3,5,5,7,7,9,11,13,15,17),
# c(2,2,4,4,6,6,8,8,10,12,14,16, 17))

gs=grid.arrange(grobs = plots2, layout_matrix = lay)
ggsave(paste0(path, "/Plots_Combined", variable, ".png"),width=18, height=4, gs)


}

# Plot all without significance threshold

########################
# Plot subcortical

plotSubcortical_noSignif=function(path, variable){
library(rgl)
  for (moda in c("thick", "LogJacs")){
 for (hemi in c("lh", "rh")){
bwasPlot=read.table(paste0(path,  "/BWAS_NVertices",variable , "_", moda, "_clustersAndCoordinates"), header=T)
bwasPlot=bwasPlot[which(bwasPlot$hemi==hemi),]
bwasPlot$X=bwasPlot$X*(-1)
bwasPlot$Y=bwasPlot$Y*(-1)

# Open phenotype file
pheno=read.table(paste0("../../UKB_phenotypes15K/", variable, ".phen"))

# Transform betas in correlations
bwasPlot$cor=bwasPlot$b/sqrt(var(pheno$V3, na.rm = T))

# Atttribute colors on a diverging palette
bwasPlot$colorScale <- cut(bwasPlot$cor, breaks = seq(-0.1, 0.1, len = 21),  include.lowest = TRUE)

## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
library(RColorBrewer)
cols=c(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[10:6],RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[5:1]) # Select palette colours
bwasPlot$color <- colorRampPalette(c(cols))(21)[bwasPlot$colorScale] # Make it more continuous
#bwasPlot$color[which(bwasPlot$signifVoxel==0)]="darkgrey"
bwasPlot$radius=0.8

# Draw plots and save screenshots
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot[,c( "Z",  "X", "Y")]), col=bwasPlot$color, radius = bwasPlot$radius)
rgl.snapshot(paste0(path, "/BWAS_", variable, "_", hemi, "_", moda , "_clustersAndCoordinates_inside_noSignif.png"))
rgl.close()

bwasPlot$X=bwasPlot$X*(-1)
bwasPlot$Z=bwasPlot$Z*(-1)
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot[,c( "Z",  "X", "Y")]), col=bwasPlot$color, radius = bwasPlot$radius)
rgl.snapshot(paste0(path, "/BWAS_", variable,"_", hemi, "_", moda , "_clustersAndCoordinates_outside_noSignif.png"))
rgl.close()

}}
}

############################
# Plot cortical
plotCortical_noSignif=function(path, variable){
library(rgl)
 for (moda in c("area", "thickness")){

 for (hemi in c("rh", "lh")){
bwasPlot=read.table(paste0(path,  "/BWAS_NVertices",variable , "_", moda, "_clustersAndCoordinates"), header=T)
bwasPlot=bwasPlot[which(bwasPlot$hemi==hemi),]

# Open phenotype file
pheno=read.table(paste0("../../UKB_phenotypes15K/", variable, ".phen"))

# Transform betas in correlations
bwasPlot$cor=bwasPlot$b/sqrt(var(pheno$V3, na.rm = T))

# Atttribute colors on a diverging palette
bwasPlot$colorScale <- cut(bwasPlot$cor, breaks = seq(-0.1, 0.1, len = 21),  include.lowest = TRUE)

## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
library(RColorBrewer)
cols=c(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[10:6], RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[5:1]) # Select palette colours
bwasPlot$color <- colorRampPalette(c(cols))(21)[bwasPlot$colorScale] # Make it more continuous
#bwasPlot$color[which(bwasPlot$signifVoxel==0)]="darkgrey"
bwasPlot$radius=0.8

#bwasPlot2=bwasPlot[sample(1:length(bwasPlot$ProbeID), size = 10000),]
#bwasPlot2=rbind(bwasPlot2, bwasPlot[which(bwasPlot$signifVoxel==1 | bwasPlot$inCluster>0 ),])
bwasPlot2=bwasPlot


# Draw plot
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot2[,c( "Y","X", "Z")]), col=bwasPlot2$color, radius = bwasPlot2$radius)
rgl.snapshot(paste0(path, "/BWAS_", variable,"_", hemi, "_", moda , "_clustersAndCoordinates_inside_noSignif.png"))
rgl.close()

bwasPlot2$X=bwasPlot2$X*(-1)
bwasPlot2$Y=bwasPlot2$Y*(-1)
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot2[,c( "Y","X", "Z")]), col=bwasPlot2$color, radius = bwasPlot2$radius)
rgl.snapshot(paste0(path, "/BWAS_", variable,"_", hemi, "_", moda , "_clustersAndCoordinates_outside_noSignif.png"))
rgl.close()

}}
}


##############################
# Combine cortical and subcortical plots
combineCorticalSubcorticalPlots_noSignif=function(path, variable){
library(grid)
library(png)
library(ggplot2)
library(gridExtra)

# List of files
ll=NULL
for (moda in c("thickness", "area","thick", "LogJacs") ){
ll=c(ll, paste0(path, "/BWAS_", variable, "_", "lh", "_", moda, "_clustersAndCoordinates_outside_noSignif.png"))
ll=c(ll, paste0(path, "/BWAS_", variable, "_","lh", "_", moda, "_clustersAndCoordinates_inside_noSignif.png"))
ll=c(ll, paste0(path, "/BWAS_", variable, "_", "rh", "_", moda, "_clustersAndCoordinates_inside_noSignif.png"))
ll=c(ll, paste0(path, "/BWAS_", variable, "_", "rh", "_", moda, "_clustersAndCoordinates_outside_noSignif.png"))
}
ll=c(ll, paste0(path, "/../legengBar.png"))

# Crop and convert format
plots2 <- lapply(ll<-ll,function(x){
  if(x!=paste0(path, "/../legengBar.png")){
    img <- as.raster(readPNG(x)[100:900,100:1100,])} else {
     img <- as.raster(readPNG(x)[,5000:5800,])
    }
    rasterGrob(img, interpolate = T)
})
#

# Lay and write png
lay <- rbind(c(1,1,3,3,5,5,7,7,9,11,13,15,17),
            c(2,2,4,4,6,6,8,8,10,12,14,16, 17))
gs=grid.arrange(grobs = plots2, layout_matrix = lay)
ggsave(paste0(path, "/Plots_Combined", variable, "_noSignif.png"),width=18, height=4, gs)

}


########################
# Plot subcortical - flat for easier visualisation of subcortical volumes

plotSubcortical_flat=function(path, variable){
library(rgl)
  for (moda in c("thick", "LogJacs")){
 for (hemi in c("lh", "rh")){
bwasPlot=read.table(paste0(path,  "/BWAS_NVertices",variable , "_", moda, "_clustersAndCoordinates"), header=T)
bwasPlot=bwasPlot[which(bwasPlot$hemi==hemi),]
bwasPlot$X=bwasPlot$X*(-1)
bwasPlot$Y=bwasPlot$Y*(-1)

# Open phenotype file
pheno=read.table(paste0("../../UKB_phenotypes15K/", variable, ".phen"))

# Transform betas in correlations
bwasPlot$cor=bwasPlot$b/sqrt(var(pheno$V3, na.rm = T))

# Atttribute colors on a diverging palette
bwasPlot$colorScale <- cut(bwasPlot$cor, breaks = seq(-0.4, 0.4, len = 21),  include.lowest = TRUE)

## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
library(RColorBrewer)
cols=c(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[10:6],RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[5:1]) # Select palette colours
bwasPlot$color <- colorRampPalette(c(cols))(21)[bwasPlot$colorScale] # Make it more continuous
bwasPlot$color[which(bwasPlot$signifVoxel==0)]="darkgrey"
bwasPlot$radius=ifelse(  bwasPlot$signifVoxel==1, 2 ,0.8 )

# Change coordinates for flat plotting
if (hemi=="rh"){

bwasPlot$Xf=bwasPlot$X
bwasPlot$Yf=bwasPlot$Y
bwasPlot$Zf=bwasPlot$Z

  # Hippocampus
bwasPlot$Yf[bwasPlot$subcv %in% c(17,53)]=bwasPlot$Yf[bwasPlot$subcv %in% c(17,53)]-17

# Amygdala
bwasPlot$Yf[bwasPlot$subcv %in% c(18,54)]=bwasPlot$Yf[bwasPlot$subcv %in% c(18,54)]-17
bwasPlot$Zf[bwasPlot$subcv %in% c(18,54)]=bwasPlot$Zf[bwasPlot$subcv %in% c(18,54)]+7

# Thalamus
bwasPlot$Zf[bwasPlot$subcv %in% c(10,49)]=bwasPlot$Zf[bwasPlot$subcv %in% c(10,49)]-20

# Caudate
bwasPlot$Yf[bwasPlot$subcv %in% c(11,50)]=bwasPlot$Yf[bwasPlot$subcv %in% c(11,50)]+20

# Accumbens
bwasPlot$Zf[bwasPlot$subcv %in% c(26,58)]=bwasPlot$Zf[bwasPlot$subcv %in% c(26,58)]+12

# Pallidum
bwasPlot$Yf[bwasPlot$subcv %in% c(13,52)]=bwasPlot$Yf[bwasPlot$subcv %in% c(13,52)]-15

# Second set of coordinates
bwasPlot$Xf2=bwasPlot$Xf*(-1)
bwasPlot$Zf2=bwasPlot$Zf*(-1)
bwasPlot$Yf2=bwasPlot$Yf

# Accumbens
bwasPlot$Zf2[bwasPlot$subcv %in% c(26,58)]=bwasPlot$Zf2[bwasPlot$subcv %in% c(26,58)]-7

# Pallidum
bwasPlot$Yf2[bwasPlot$subcv %in% c(13,52)]=bwasPlot$Yf2[bwasPlot$subcv %in% c(13,52)]-8

# Amygdala
bwasPlot$Yf2[bwasPlot$subcv %in% c(18,54)]=bwasPlot$Yf2[bwasPlot$subcv %in% c(18,54)]-5

# Hippocampus
bwasPlot$Yf2[bwasPlot$subcv %in% c(17,53)]=bwasPlot$Yf2[bwasPlot$subcv %in% c(17,53)]-5

# Caudate
bwasPlot$Yf2[bwasPlot$subcv %in% c(11,50)]=bwasPlot$Yf2[bwasPlot$subcv %in% c(11,50)]-8
##########################
} else if (hemi == "lh"){
# Second set of coordinates

bwasPlot$Xf=bwasPlot$X
bwasPlot$Yf=bwasPlot$Y
bwasPlot$Zf=bwasPlot$Z

  # Hippocampus
bwasPlot$Yf[bwasPlot$subcv %in% c(17,53)]=bwasPlot$Yf[bwasPlot$subcv %in% c(17,53)]-21

# Amygdala
bwasPlot$Yf[bwasPlot$subcv %in% c(18,54)]=bwasPlot$Yf[bwasPlot$subcv %in% c(18,54)]-21
bwasPlot$Zf[bwasPlot$subcv %in% c(18,54)]=bwasPlot$Zf[bwasPlot$subcv %in% c(18,54)]+7

# Thalamus
bwasPlot$Zf[bwasPlot$subcv %in% c(10,49)]=bwasPlot$Zf[bwasPlot$subcv %in% c(10,49)]-20

# Caudate
bwasPlot$Yf[bwasPlot$subcv %in% c(11,50)]=bwasPlot$Yf[bwasPlot$subcv %in% c(11,50)]+12

# Accumbens
bwasPlot$Zf[bwasPlot$subcv %in% c(26,58)]=bwasPlot$Zf[bwasPlot$subcv %in% c(26,58)]+19

# Pallidum
bwasPlot$Yf[bwasPlot$subcv %in% c(13,52)]=bwasPlot$Yf[bwasPlot$subcv %in% c(13,52)]-22

# Second set of coordinates
bwasPlot$Xf2=bwasPlot$Xf*(-1)
bwasPlot$Zf2=bwasPlot$Zf*(-1)
bwasPlot$Yf2=bwasPlot$Yf

# Caudate
bwasPlot$Yf2[bwasPlot$subcv %in% c(11,50)]=bwasPlot$Yf2[bwasPlot$subcv %in% c(11,50)]+8

# Pallidum
bwasPlot$Yf2[bwasPlot$subcv %in% c(13,52)]=bwasPlot$Yf2[bwasPlot$subcv %in% c(13,52)]+6

# Amygdala
bwasPlot$Yf2[bwasPlot$subcv %in% c(18,54)]=bwasPlot$Yf2[bwasPlot$subcv %in% c(18,54)]+3

# Hippocampus
bwasPlot$Yf2[bwasPlot$subcv %in% c(17,53)]=bwasPlot$Yf2[bwasPlot$subcv %in% c(17,53)]+3

# Accumbens
bwasPlot$Zf2[bwasPlot$subcv %in% c(26,58)]=bwasPlot$Zf2[bwasPlot$subcv %in% c(26,58)]+6

# Thalamus
bwasPlot$Zf2[bwasPlot$subcv %in% c(10,49)]=bwasPlot$Zf2[bwasPlot$subcv %in% c(10,49)]+2

}

# Draw plots and save screenshots
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot[,c( "Zf",  "Xf", "Yf")]), col=bwasPlot$color, radius = bwasPlot$radius)
rgl.snapshot(paste0(path, "/BWAS_", variable, "_", hemi, "_", moda , "_clustersAndCoordinates_inside.png"))
rgl.close()


par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot[,c( "Zf2",  "Xf2", "Yf2")]), col=bwasPlot$color, radius = bwasPlot$radius)
rgl.snapshot(paste0(path, "/BWAS_", variable,"_", hemi, "_", moda , "_clustersAndCoordinates_outside.png"))
rgl.close()

}}
}


# Plot subcortical - flat for easier visualisation of subcortical volumes

plotSubcortical_flat_multiModel=function(path, variable, otherPaths2Models){
library(rgl)
  for (moda in c("thick", "LogJacs")){
 for (hemi in c("lh", "rh")){
bwasPlot=read.table(paste0(path,  "/BWAS_NVertices",variable , "_", moda, "_clustersAndCoordinates"), header=T)

bwasPlot=bwasPlot[which(bwasPlot$hemi==hemi),]
bwasPlot$X=bwasPlot$X*(-1)
bwasPlot$Y=bwasPlot$Y*(-1)

for (model in otherPaths2Models){

  bwasPlotOther=read.table(paste0(model,  "/BWAS_NVertices",variable , "_", moda, "_clustersAndCoordinates"), header=T)
  bwasPlotOther=bwasPlotOther[which(bwasPlotOther$hemi==hemi),]
  bwasPlotOther=bwasPlotOther[,c("Probe", "signifVoxel")]
  colnames(bwasPlotOther)[2]=paste0("signifVoxel_Other")
  bwasPlot=merge(bwasPlot, bwasPlotOther, by="Probe" )

}

# Open phenotype file
pheno=read.table(paste0("UKB_phenotypes15K/", variable, ".phen"))

# Transform betas in correlations
bwasPlot$cor=bwasPlot$b/sqrt(var(pheno$V3, na.rm = T))

# Atttribute colors on a diverging palette
bwasPlot$colorScale <- cut(bwasPlot$cor, breaks = seq(-0.1, 0.1, len = 21),  include.lowest = TRUE)

## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
library(RColorBrewer)
cols=c(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[10:6],RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[5:1]) # Select palette colours
bwasPlot$color <- colorRampPalette(c(cols))(21)[bwasPlot$colorScale] # Make it more continuous
bwasPlot$signifVoxel_All=rowSums(bwasPlot[,grep(pattern = "signifVoxel", x = colnames(bwasPlot) )] )

bwasPlot$color[which(bwasPlot$signifVoxel_All==0)]="darkgrey"
bwasPlot$radius=0.4
bwasPlot$radius[which(bwasPlot$signifVoxel_All==2)]=0.6
bwasPlot$radius[which(bwasPlot$signifVoxel_All==3)]=0.8
bwasPlot$radius[which(bwasPlot$signifVoxel_All==4)]=1.1
bwasPlot$radius[which(bwasPlot$signifVoxel_All==5)]=1.3
bwasPlot$radius[which(bwasPlot$signifVoxel_All==6)]=1.5

# Change coordinates for flat plotting
if (hemi=="rh"){

bwasPlot$Xf=bwasPlot$X
bwasPlot$Yf=bwasPlot$Y
bwasPlot$Zf=bwasPlot$Z

  # Hippocampus
bwasPlot$Yf[bwasPlot$subcv %in% c(17,53)]=bwasPlot$Yf[bwasPlot$subcv %in% c(17,53)]-17

# Amygdala
bwasPlot$Yf[bwasPlot$subcv %in% c(18,54)]=bwasPlot$Yf[bwasPlot$subcv %in% c(18,54)]-17
bwasPlot$Zf[bwasPlot$subcv %in% c(18,54)]=bwasPlot$Zf[bwasPlot$subcv %in% c(18,54)]+7

# Thalamus
bwasPlot$Zf[bwasPlot$subcv %in% c(10,49)]=bwasPlot$Zf[bwasPlot$subcv %in% c(10,49)]-20

# Caudate
bwasPlot$Yf[bwasPlot$subcv %in% c(11,50)]=bwasPlot$Yf[bwasPlot$subcv %in% c(11,50)]+20

# Accumbens
bwasPlot$Zf[bwasPlot$subcv %in% c(26,58)]=bwasPlot$Zf[bwasPlot$subcv %in% c(26,58)]+12

# Pallidum
bwasPlot$Yf[bwasPlot$subcv %in% c(13,52)]=bwasPlot$Yf[bwasPlot$subcv %in% c(13,52)]-15

# Second set of coordinates
bwasPlot$Xf2=bwasPlot$Xf*(-1)
bwasPlot$Zf2=bwasPlot$Zf*(-1)
bwasPlot$Yf2=bwasPlot$Yf

# Accumbens
bwasPlot$Zf2[bwasPlot$subcv %in% c(26,58)]=bwasPlot$Zf2[bwasPlot$subcv %in% c(26,58)]-7

# Pallidum
bwasPlot$Yf2[bwasPlot$subcv %in% c(13,52)]=bwasPlot$Yf2[bwasPlot$subcv %in% c(13,52)]-8

# Amygdala
bwasPlot$Yf2[bwasPlot$subcv %in% c(18,54)]=bwasPlot$Yf2[bwasPlot$subcv %in% c(18,54)]-5

# Hippocampus
bwasPlot$Yf2[bwasPlot$subcv %in% c(17,53)]=bwasPlot$Yf2[bwasPlot$subcv %in% c(17,53)]-5

# Caudate
bwasPlot$Yf2[bwasPlot$subcv %in% c(11,50)]=bwasPlot$Yf2[bwasPlot$subcv %in% c(11,50)]-8
##########################
} else if (hemi == "lh"){
# Second set of coordinates

bwasPlot$Xf=bwasPlot$X
bwasPlot$Yf=bwasPlot$Y
bwasPlot$Zf=bwasPlot$Z

  # Hippocampus
bwasPlot$Yf[bwasPlot$subcv %in% c(17,53)]=bwasPlot$Yf[bwasPlot$subcv %in% c(17,53)]-21

# Amygdala
bwasPlot$Yf[bwasPlot$subcv %in% c(18,54)]=bwasPlot$Yf[bwasPlot$subcv %in% c(18,54)]-21
bwasPlot$Zf[bwasPlot$subcv %in% c(18,54)]=bwasPlot$Zf[bwasPlot$subcv %in% c(18,54)]+7

# Thalamus
bwasPlot$Zf[bwasPlot$subcv %in% c(10,49)]=bwasPlot$Zf[bwasPlot$subcv %in% c(10,49)]-20

# Caudate
bwasPlot$Yf[bwasPlot$subcv %in% c(11,50)]=bwasPlot$Yf[bwasPlot$subcv %in% c(11,50)]+12

# Accumbens
bwasPlot$Zf[bwasPlot$subcv %in% c(26,58)]=bwasPlot$Zf[bwasPlot$subcv %in% c(26,58)]+19

# Pallidum
bwasPlot$Yf[bwasPlot$subcv %in% c(13,52)]=bwasPlot$Yf[bwasPlot$subcv %in% c(13,52)]-22

# Second set of coordinates
bwasPlot$Xf2=bwasPlot$Xf*(-1)
bwasPlot$Zf2=bwasPlot$Zf*(-1)
bwasPlot$Yf2=bwasPlot$Yf

# Caudate
bwasPlot$Yf2[bwasPlot$subcv %in% c(11,50)]=bwasPlot$Yf2[bwasPlot$subcv %in% c(11,50)]+8

# Pallidum
bwasPlot$Yf2[bwasPlot$subcv %in% c(13,52)]=bwasPlot$Yf2[bwasPlot$subcv %in% c(13,52)]+6

# Amygdala
bwasPlot$Yf2[bwasPlot$subcv %in% c(18,54)]=bwasPlot$Yf2[bwasPlot$subcv %in% c(18,54)]+3

# Hippocampus
bwasPlot$Yf2[bwasPlot$subcv %in% c(17,53)]=bwasPlot$Yf2[bwasPlot$subcv %in% c(17,53)]+3

# Accumbens
bwasPlot$Zf2[bwasPlot$subcv %in% c(26,58)]=bwasPlot$Zf2[bwasPlot$subcv %in% c(26,58)]+6

# Thalamus
bwasPlot$Zf2[bwasPlot$subcv %in% c(10,49)]=bwasPlot$Zf2[bwasPlot$subcv %in% c(10,49)]+2

}

# Draw plots and save screenshots
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot[,c( "Zf",  "Xf", "Yf")]), col=bwasPlot$color, radius = bwasPlot$radius, alpha=0.5)
rgl.snapshot(paste0(path, "/BWAS_", variable, "_", hemi, "_", moda , "_clustersAndCoordinates_inside_multiModel.png"))
rgl.close()


par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot[,c( "Zf2",  "Xf2", "Yf2")]), col=bwasPlot$color, radius = bwasPlot$radius)
rgl.snapshot(paste0(path, "/BWAS_", variable,"_", hemi, "_", moda , "_clustersAndCoordinates_outside_multiModel.png"))
rgl.close()

}}
}





#####################################################
#
identifyClustersBWAS=function(pathToTheBWASresults, variable, outputFolder, variableLabel, scenario){

library(Morpho)
library(plyr)
library(viridis)

# read arguments provided when submitting batch job
# arg = commandArgs(trailingOnly=TRUE)
arg=c(pathToTheBWASresults, variable, outputFolder, variableLabel, scenario)

# Arguments are
# 1) pathToTheBWASresults/scenario
# 2) variable
# 3) output folder
# 4) variable label

signifThreshold=0.05/652283/5
# 1.5e-8
# Open list of associated vertices

if (file.exists(paste0(arg[1], "/BWAS_", arg[2], ".linear"))){
bwasResult=read.table(paste0(arg[1], "/BWAS_", arg[2], ".linear"), header=T, stringsAsFactors = F) } else if (file.exists(paste0(arg[1], "/BWAS_", arg[2], ".moa"))) {
bwasResult=read.table(paste0(arg[1], "/BWAS_", arg[2], ".moa"), header=T, stringsAsFactors = F)  } else {
bwasResult=read.table(paste0(arg[1], "/BWAS_", arg[2], ".moment"), header=T, stringsAsFactors = F)
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
res=clusterIdentification(path=arg[1], bwas = bwas,  moda = moda, signifThreshold = signifThreshold, KNearestNeighbours = 10, phenotypeLabel = arg[2])
colnames(res)=paste0(colnames(res), "_", moda)
resTot=c(resTot, res)
}

# Write final output - with all iterations
write.table(resTot, paste0(arg[3], "/Results_clusterFWER_", arg[2] ,  ".txt" ))


}




atlasPath="path/to/atlas/folder"


#####################
# Functions
formatBWAScortical<-function(BWASsumstat, hemi, mod){
BWASsumstat$ProbeID<-BWASsumstat$Probe


atlas<-read.table(paste0(atlasPath, "/fsaverage_vertex_labels_names_order2.txt"), header=T, stringsAsFactors = F)
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


atlas<-read.table(paste0(atlasPath, "/fsaverage_subcortical_vertices_label.txt"), header=T, stringsAsFactors = F)
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

BWASPlotSimpleSimulatedPheno<-function(path , iter, NNN, yMax, mlmaOrLinear, ntotvar, typeAssociatedVertices){
library(png)
library(plyr)
library(qqman)

maxx<-yMax
BWASsignif=NULL
chi2All=NULL
pAll=NULL
lambda=NULL
signifLog10= -log( (0.05/652283)/ntotvar, base=10 )

# Open list of true associations
TP=read.table(paste0(path, "/../02_phenotypesSimul/NVertices", NNN, "_", typeAssociatedVertices, "_Niter", iter, ".par"), header = T, stringsAsFactors = F)
TP$vertexNum=substr(TP$Causal, 5, nchar(TP$Causal))
TP$hemi=substr(TP$Causal, 1,2)

png(paste0(path, "/Manhathan_", typeAssociatedVertices, "_N",NNN, "_iter", iter, "_simple.png"), width = 50, height = 16, units = "cm", res = 400, type="cairo")
plot.new()
layout(matrix(c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,6,7,8), nrow =1 , ncol = 20, byrow = F))

 BWASsumstatfilePath = paste0(path ,"/BWAS_NVertices",NNN,"_", typeAssociatedVertices, "_Niter",iter,".", mlmaOrLinear)
    BWASsumstatfile<-  read.table(BWASsumstatfilePath, header=T, stringsAsFactors = F)
    BWASsumstatfile$P<-BWASsumstatfile$p
    BWASsumstatfile$log10p<- -log10(BWASsumstatfile$P)

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

sum(table(BWASsumstat$subcvName))


  # Get all chi2 and pvalues
BWASsumstat$CHI2=(BWASsumstat$b/BWASsumstat$se)**2
chi2All=c(chi2All, BWASsumstat$CHI2)
pAll=c(pAll, BWASsumstat$p)

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
plot(BWASsumstat$xax, BWASsumstat$log10p, col=BWASsumstat$col, pch=20, cex=1,  xlab="", xaxt="n", cex.lab=1, ylim=c(0, maxx+1), ylab="-log10(p-value)", bty="n", main="", cex.main=0.7, xaxs="i", yaxs="i", xlim=c(-10, max(BWASsumstat$xax)+10))
points(BWASsumstat$xax[which(BWASsumstat$Probe %in% TP$Causal)], BWASsumstat$log10p[which(BWASsumstat$Probe %in% TP$Causal)], pch=17, cex=1.5, col="darkgrey")
title("Cortical thickness \n left", line = -1, cex.main=0.7)
title(paste0("N", NNN, "_iter", iter), line = 1, cex.main=0.9, adj=0)
} else{

  if(mod=="area"){
      moda<- "Cortical area"
} else if (mod=="thickness"){
     moda<- "Cortical thickness"
}else if (mod=="LogJacs"){
     moda<- "Subcortical area"
}else if (mod=="thick"){
     moda<- "Subcortical thickness" }
  titl<-paste( moda , "\n",ifelse(hemi=="lh", "Left", "Right"))
  par(mar=c(8,0.5,2,0.5))
  plot(BWASsumstat$xax, BWASsumstat$log10p, col=BWASsumstat$col, pch=20, cex=1,  xlab="", xaxt="n",  yaxt="n", cex.lab=0.9, ylim=c(0, maxx+1), ylab="", bty="n", main="", cex.main=0.7, xaxs="i", yaxs="i", xlim=c(-10, max(BWASsumstat$xax)+10) )
  points(BWASsumstat$xax[which(BWASsumstat$Probe %in% TP$Causal)], BWASsumstat$log10p[which(BWASsumstat$Probe %in% TP$Causal)], pch=17, cex=1.5, col="darkgrey")
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

write.csv(BWASsignif, paste0(path, "_BWAS_fullSummary_N",NNN, "_iter", iter, ".csv" ))
write.csv(BWASsignif[which(BWASsignif$log10p > signifLog10),], paste0(path, "_BWAS_signif_N",NNN, "_iter", iter, ".csv" ))


}


#####################
summaryBWAS<-function(path ,  mlmaOrLinear, ntotvar, iter, NNN, typeAssociatedVertices){

BWASsignif=NULL
signifLog10= -log( (0.05/652283)/ntotvar, base=10 )
chi2All=NULL
res=NULL
pAll=NULL
lambda=NULL

# Open files and format them

    BWASsumstatfilePath = paste0(path ,"/BWAS_NVertices",NNN,"_", typeAssociatedVertices, "_Niter",iter,".", mlmaOrLinear)
    BWASsumstatfile<-  read.table(BWASsumstatfilePath, header=T, stringsAsFactors = F)
    BWASsumstatfile$P<-BWASsumstatfile$p
    BWASsumstatfile$log10p<- -log10(BWASsumstatfile$P)
    associatedVertices=read.table(paste0(path ,"/../01_vertexSelection/NVertices",NNN,"_", typeAssociatedVertices, "_Niter",iter))

    for (mod in c("area", "thickness", "thick", "LogJacs")){

BWASsumstat=NULL
  if (mod %in% c("area", "thickness")){

    for (hemi in c("lh", "rh")){
BWASsumstatHemi<-formatBWAScortical(BWASsumstat = BWASsumstatfile,hemi = hemi, mod=mod)
BWASsumstat=rbind(BWASsumstat,BWASsumstatHemi)
    }
  }
   if (mod %in% c("LogJacs", "thick")){
BWASsumstat<-formatBWASsubcortical(BWASsumstat = BWASsumstatfile)
BWASsumstat=BWASsumstat[grep(BWASsumstat$Probe, pattern = mod),]
   }

# Calculate values of interest
BWASsumstat$CHI2=(BWASsumstat$b/BWASsumstat$se)**2
BWASsumstatTP=BWASsumstat[which(BWASsumstat$Probe %in% associatedVertices$V1),]
BWASsumstatNA=BWASsumstat[which(!BWASsumstat$Probe %in% associatedVertices$V1),]


res=cbind(res, cbind(signif(mean(BWASsumstat$CHI2),3 ), signif( median(BWASsumstat$CHI2)/qchisq(0.5,1),3), ifelse(length(which(BWASsumstatNA$log10p>signifLog10))>0, 1, 0 ), length(which(BWASsumstatNA$log10p>signifLog10))/length(BWASsumstatNA$log10p), length(which(BWASsumstatNA$P<0.05))/length(BWASsumstatNA$log10p) , length(which(BWASsumstatTP$log10p>signifLog10))/length(BWASsumstatTP$log10p) ))

# get all chi2 and pvalues
chi2All=rbind(chi2All, cbind(BWASsumstat$CHI2, BWASsumstat$log10p,BWASsumstat$P,  as.character(BWASsumstat$Probe), mod) )
  }

# Get global metrics
chi2All=as.data.frame(chi2All)
colnames(chi2All)=c("CHI2", "log10p","P", "Probe", "moda")
chi2All$CHI2=as.numeric(levels(chi2All$CHI2))[chi2All$CHI2]
chi2All$log10p=as.numeric(levels(chi2All$log10p))[chi2All$log10p]
chi2All$P=as.numeric(levels(chi2All$P))[chi2All$P]

chi2All$isAssociated=ifelse(chi2All$Probe %in% associatedVertices$V1, 1, 0)
chi2All$isSignificant=ifelse(chi2All$log10p > signifLog10, 1, 0)

chi2AllTP=chi2All[which(chi2All$Probe %in% associatedVertices$V1),]
chi2AllNA=chi2All[-which(chi2All$Probe %in% associatedVertices$V1),]

chi2OtherModa=chi2All[-which(chi2All$moda %in% typeAssociatedVertices),]

res=cbind(res, cbind(signif(mean(chi2AllNA$CHI2),3 ), signif( median(chi2AllNA$CHI2)/qchisq(0.5,1),3), ifelse(length(which(chi2AllNA$log10p>signifLog10))>0, 1, 0 ), length(which(chi2AllNA$log10p>signifLog10)), length(which(chi2AllNA$P<0.05))/length(chi2AllNA$log10p) , length(which(chi2AllTP$log10p>signifLog10))/length(chi2AllTP$log10p), auc(chi2All$isAssociated, chi2All$isSignificant), length(which(chi2OtherModa$log10p>signifLog10)), signif( median(chi2OtherModa$CHI2)/qchisq(0.5,1),3), length(which(chi2OtherModa$P<0.05))/length(chi2OtherModa$P) ))

colnames(res)=c(colnames(res), paste0(c("meanX2H0", "lambdaH0", "FWER", "nFP", "FPR_0.05_H0", "TPR"), rep(c("_area", "_thickness", "_thick", "_LogJacs", "_total"), each=6) ), "AUC_total", "nFP_OtherModa" , "lambdaOtherModa","FPR_0.05_OtherModa" )

return(res)
}

#####################################
# All moda
summaryBWAS_allModa<-function(path , iter, NNN){

BWASsignif=NULL
signifLog10= -log( (0.05/652283), base=10 )
chi2All=NULL
res=NULL
pAll=NULL
lambda=NULL

# Open files and format them
if (file.exists(paste0(path ,"/BWAS_NVertices",NNN,"_Niter",iter, ".linear"))){
BWASsumstatfile=read.table(paste0(path ,"/BWAS_NVertices",NNN,"_Niter",iter, ".linear"), header=T, stringsAsFactors = F) } else if (file.exists(paste0(path ,"/BWAS_NVertices",NNN,"_Niter",iter, ".moa"))){
BWASsumstatfile=read.table(paste0(path ,"/BWAS_NVertices",NNN,"_Niter",iter, ".moa"), header=T, stringsAsFactors = F)  } else {
BWASsumstatfile=read.table(paste0(path ,"/BWAS_NVertices",NNN,"_Niter",iter, ".moment"), header=T, stringsAsFactors = F)
}
    BWASsumstatfile$P<-BWASsumstatfile$p
    BWASsumstatfile$log10p<- -log10(BWASsumstatfile$P)
    associatedVertices=read.table(paste0(path ,"/../01_vertexSelection/NVertices",NNN,"_Niter",iter))

    for (mod in c("area", "thickness", "thick", "LogJacs")){

BWASsumstat=NULL
  if (mod %in% c("area", "thickness")){

    for (hemi in c("lh", "rh")){
BWASsumstatHemi<-formatBWAScortical(BWASsumstat = BWASsumstatfile,hemi = hemi, mod=mod)
BWASsumstat=rbind(BWASsumstat,BWASsumstatHemi)
    }
  }
   if (mod %in% c("LogJacs", "thick")){
BWASsumstat<-formatBWASsubcortical(BWASsumstat = BWASsumstatfile)
BWASsumstat=BWASsumstat[grep(BWASsumstat$Probe, pattern = mod),]
   }

# Calculate values of interest
BWASsumstat$CHI2=(BWASsumstat$b/BWASsumstat$se)**2
BWASsumstatTP=BWASsumstat[which(BWASsumstat$Probe %in% associatedVertices$V1),]
BWASsumstatNA=BWASsumstat[which(!BWASsumstat$Probe %in% associatedVertices$V1),]


res=cbind(res, cbind(signif(mean(BWASsumstat$CHI2),3 ), signif( median(BWASsumstat$CHI2)/qchisq(0.5,1),3), ifelse(length(which(BWASsumstatNA$log10p>signifLog10))>0, 1, 0 ), length(which(BWASsumstatNA$log10p>signifLog10))/length(BWASsumstatNA$log10p), length(which(BWASsumstatNA$P<0.05))/length(BWASsumstatNA$log10p) , length(which(BWASsumstatTP$log10p>signifLog10))/length(BWASsumstatTP$log10p) ))

# get all chi2 and pvalues
chi2All=rbind(chi2All, cbind(BWASsumstat$CHI2, BWASsumstat$log10p,BWASsumstat$P,  as.character(BWASsumstat$Probe), mod) )
  }

# Get global metrics
chi2All=as.data.frame(chi2All)
colnames(chi2All)=c("CHI2", "log10p","P", "Probe", "moda")
chi2All$CHI2=as.numeric(levels(chi2All$CHI2))[chi2All$CHI2]
chi2All$log10p=as.numeric(levels(chi2All$log10p))[chi2All$log10p]
chi2All$P=as.numeric(levels(chi2All$P))[chi2All$P]

chi2All$isAssociated=ifelse(chi2All$Probe %in% associatedVertices$V1, 1, 0)
chi2All$isSignificant=ifelse(chi2All$log10p > signifLog10, 1, 0)

chi2AllTP=chi2All[which(chi2All$Probe %in% associatedVertices$V1),]
chi2AllNA=chi2All[-which(chi2All$Probe %in% associatedVertices$V1),]

chi2OtherModa=chi2All

res=cbind(res, cbind(signif(mean(chi2AllNA$CHI2),3 ), signif( median(chi2AllNA$CHI2)/qchisq(0.5,1),3), ifelse(length(which(chi2AllNA$log10p>signifLog10))>0, 1, 0 ), length(which(chi2AllNA$log10p>signifLog10)), length(which(chi2AllNA$P<0.05))/length(chi2AllNA$log10p) , length(which(chi2AllTP$log10p>signifLog10))/length(chi2AllTP$log10p), auc(chi2All$isAssociated, chi2All$isSignificant), length(which(chi2OtherModa$log10p>signifLog10)), signif( median(chi2OtherModa$CHI2)/qchisq(0.5,1),3), length(which(chi2OtherModa$P<0.05))/length(chi2OtherModa$P),  mean( roc(chi2All$isAssociated, chi2All$isSignificant)$sensitivities[2], roc(chi2All$isAssociated, chi2All$isSignificant)$specificities[2] ) ))

colnames(res)=c(colnames(res), paste0(c("meanX2H0", "lambdaH0", "FWER", "nFP", "FPR_0.05_H0", "TPR"), rep(c("_area", "_thickness", "_thick", "_LogJacs", "_total"), each=6) ), "AUC_total", "nFP_OtherModa" , "lambdaOtherModa","FPR_0.05_OtherModa", "BalAccuracy" )

return(res)
}

#################################
# FWHM20
# All moda
summaryBWAS_allModa_fwhm20<-function(path , iter, NNN){

BWASsignif=NULL
signifLog10= -log( (0.05/652283), base=10 )
chi2All=NULL
res=NULL
pAll=NULL
lambda=NULL

# Open files and format them
if (file.exists(paste0(path ,"/BWAS_NVertices",NNN,"_Niter",iter, "_fwhm20.linear"))){
BWASsumstatfile=read.table(paste0(path ,"/BWAS_NVertices",NNN,"_Niter",iter, "_fwhm20.linear"), header=T, stringsAsFactors = F) } else if (file.exists(paste0(path ,"/BWAS_NVertices",NNN,"_Niter",iter, "_fwhm20.moa"))){
BWASsumstatfile=read.table(paste0(path ,"/BWAS_NVertices",NNN,"_Niter",iter, "_fwhm20.moa"), header=T, stringsAsFactors = F)  } else {
BWASsumstatfile=read.table(paste0(path ,"/BWAS_NVertices",NNN,"_Niter",iter, "_fwhm20.moment"), header=T, stringsAsFactors = F)
}
    BWASsumstatfile$P<-BWASsumstatfile$p
    BWASsumstatfile$log10p<- -log10(BWASsumstatfile$P)
    associatedVertices=read.table(paste0(path ,"/../01_vertexSelection/NVertices",NNN,"_Niter",iter))

    for (mod in c("area", "thickness", "thick", "LogJacs")){

BWASsumstat=NULL
  if (mod %in% c("area", "thickness")){

    for (hemi in c("lh", "rh")){
BWASsumstatHemi<-formatBWAScortical(BWASsumstat = BWASsumstatfile,hemi = hemi, mod=mod)
BWASsumstat=rbind(BWASsumstat,BWASsumstatHemi)
    }
  }
   if (mod %in% c("LogJacs", "thick")){
BWASsumstat<-formatBWASsubcortical(BWASsumstat = BWASsumstatfile)
BWASsumstat=BWASsumstat[grep(BWASsumstat$Probe, pattern = mod),]
   }

# Calculate values of interest
BWASsumstat$CHI2=(BWASsumstat$b/BWASsumstat$se)**2
BWASsumstatTP=BWASsumstat[which(BWASsumstat$Probe %in% associatedVertices$V1),]
BWASsumstatNA=BWASsumstat[which(!BWASsumstat$Probe %in% associatedVertices$V1),]


res=cbind(res, cbind(signif(mean(BWASsumstat$CHI2),3 ), signif( median(BWASsumstat$CHI2)/qchisq(0.5,1),3), ifelse(length(which(BWASsumstatNA$log10p>signifLog10))>0, 1, 0 ), length(which(BWASsumstatNA$log10p>signifLog10))/length(BWASsumstatNA$log10p), length(which(BWASsumstatNA$P<0.05))/length(BWASsumstatNA$log10p) , length(which(BWASsumstatTP$log10p>signifLog10))/length(BWASsumstatTP$log10p) ))

# get all chi2 and pvalues
chi2All=rbind(chi2All, cbind(BWASsumstat$CHI2, BWASsumstat$log10p,BWASsumstat$P,  as.character(BWASsumstat$Probe), mod) )
  }

# Get global metrics
chi2All=as.data.frame(chi2All)
colnames(chi2All)=c("CHI2", "log10p","P", "Probe", "moda")
chi2All$CHI2=as.numeric(levels(chi2All$CHI2))[chi2All$CHI2]
chi2All$log10p=as.numeric(levels(chi2All$log10p))[chi2All$log10p]
chi2All$P=as.numeric(levels(chi2All$P))[chi2All$P]

chi2All$isAssociated=ifelse(chi2All$Probe %in% associatedVertices$V1, 1, 0)
chi2All$isSignificant=ifelse(chi2All$log10p > signifLog10, 1, 0)

chi2AllTP=chi2All[which(chi2All$Probe %in% associatedVertices$V1),]
chi2AllNA=chi2All[-which(chi2All$Probe %in% associatedVertices$V1),]

chi2OtherModa=chi2All

res=cbind(res, cbind(signif(mean(chi2AllNA$CHI2),3 ), signif( median(chi2AllNA$CHI2)/qchisq(0.5,1),3), ifelse(length(which(chi2AllNA$log10p>signifLog10))>0, 1, 0 ), length(which(chi2AllNA$log10p>signifLog10)), length(which(chi2AllNA$P<0.05))/length(chi2AllNA$log10p) , length(which(chi2AllTP$log10p>signifLog10))/length(chi2AllTP$log10p), auc(chi2All$isAssociated, chi2All$isSignificant), length(which(chi2OtherModa$log10p>signifLog10)), signif( median(chi2OtherModa$CHI2)/qchisq(0.5,1),3), length(which(chi2OtherModa$P<0.05))/length(chi2OtherModa$P) ))

colnames(res)=c(colnames(res), paste0(c("meanX2H0", "lambdaH0", "FWER", "nFP", "FPR_0.05_H0", "TPR"), rep(c("_area", "_thickness", "_thick", "_LogJacs", "_total"), each=6) ), "AUC_total", "nFP_OtherModa" , "lambdaOtherModa","FPR_0.05_OtherModa" )

return(res)
}

###################################
# Null phenotypes
# All moda
  #path="/gpfs1/scratch/30days/uqbcouvy/56_BWAS/03_BWAS_uncorrected_allModa"
  #iter=1
  #NNN=0

summaryBWAS_allModa_NULL<-function(path , iter, NNN){

BWASsignif=NULL
signifLog10= -log( (0.05/652283), base=10 )
chi2All=NULL
res=NULL
pAll=NULL
lambda=NULL

# Open files and format them
if (file.exists(paste0(path ,"/BWAS_NVertices",NNN,"_Niter",iter, ".linear"))){
BWASsumstatfile=read.table(paste0(path ,"/BWAS_NVertices",NNN,"_Niter",iter, ".linear"), header=T, stringsAsFactors = F) } else if (file.exists(paste0(path ,"/BWAS_NVertices",NNN,"_Niter",iter, ".moa"))){
BWASsumstatfile=read.table(paste0(path ,"/BWAS_NVertices",NNN,"_Niter",iter, ".moa"), header=T, stringsAsFactors = F)  } else {
BWASsumstatfile=read.table(paste0(path ,"/BWAS_NVertices",NNN,"_Niter",iter, ".moment"), header=T, stringsAsFactors = F)
}
    BWASsumstatfile$P<-BWASsumstatfile$p
    BWASsumstatfile$log10p<- -log10(BWASsumstatfile$P)

    for (mod in c("area", "thickness", "thick", "LogJacs")){

BWASsumstat=NULL
  if (mod %in% c("area", "thickness")){

    for (hemi in c("lh", "rh")){
BWASsumstatHemi<-formatBWAScortical(BWASsumstat = BWASsumstatfile,hemi = hemi, mod=mod)
BWASsumstat=rbind(BWASsumstat,BWASsumstatHemi)
    }
  }
   if (mod %in% c("LogJacs", "thick")){
BWASsumstat<-formatBWASsubcortical(BWASsumstat = BWASsumstatfile)
BWASsumstat=BWASsumstat[grep(BWASsumstat$Probe, pattern = mod),]
   }

# Calculate values of interest
BWASsumstat$CHI2=(BWASsumstat$b/BWASsumstat$se)**2
BWASsumstatNA=BWASsumstat

res=cbind(res, cbind(signif(mean(BWASsumstat$CHI2),3 ), signif( median(BWASsumstat$CHI2)/qchisq(0.5,1),3), ifelse(length(which(BWASsumstatNA$log10p>signifLog10))>0, 1, 0 ), length(which(BWASsumstatNA$log10p>signifLog10))/length(BWASsumstatNA$log10p), length(which(BWASsumstatNA$P<0.05))/length(BWASsumstatNA$log10p) , NA ))

# get all chi2 and pvalues
chi2All=rbind(chi2All, cbind(BWASsumstat$CHI2, BWASsumstat$log10p,BWASsumstat$P,  as.character(BWASsumstat$Probe), mod) )
  }

# Get global metrics
chi2All=as.data.frame(chi2All)
colnames(chi2All)=c("CHI2", "log10p","P", "Probe", "moda")
chi2All$CHI2=as.numeric(levels(chi2All$CHI2))[chi2All$CHI2]
chi2All$log10p=as.numeric(levels(chi2All$log10p))[chi2All$log10p]
chi2All$P=as.numeric(levels(chi2All$P))[chi2All$P]

chi2All$isAssociated=0
chi2All$isSignificant=ifelse(chi2All$log10p > signifLog10, 1, 0)
chi2AllNA=chi2All
chi2OtherModa=chi2All

res=cbind(res, cbind(signif(mean(chi2AllNA$CHI2),3 ), signif( median(chi2AllNA$CHI2)/qchisq(0.5,1),3), ifelse(length(which(chi2AllNA$log10p>signifLog10))>0, 1, 0 ), length(which(chi2AllNA$log10p>signifLog10)), length(which(chi2AllNA$P<0.05))/length(chi2AllNA$log10p) , NA, NA, length(which(chi2OtherModa$log10p>signifLog10)), signif( median(chi2OtherModa$CHI2)/qchisq(0.5,1),3), length(which(chi2OtherModa$P<0.05))/length(chi2OtherModa$P) ))

colnames(res)=c(colnames(res), paste0(c("meanX2H0", "lambdaH0", "FWER", "nFP", "FPR_0.05_H0", "TPR"), rep(c("_area", "_thickness", "_thick", "_LogJacs", "_total"), each=6) ), "AUC_total", "nFP_OtherModa" , "lambdaOtherModa","FPR_0.05_OtherModa" )

return(res)
}


##################################
# Add coordinates to formated BWAS file

bwasAddCoordinatesCortical=function(bwasFormatted, hemi, moda){
# Get vertices coordinates
coor=read.table(paste0(atlasPath ,"/", hemi, ".vertexCoordinates.asc"))
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

significantClusterAroundVertex=function(bwas, signifThreshold, truePositiveLabel, TPorFP, nbNearestNeighbours ){
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

bwas[,paste0("cluster_", TPorFP, "_", truePositiveLabel)]=ifelse(bwas$voxelStatus %in% c("belongInCluster", "truePositive"), 1 , 0 )
print(paste("Found ", sum(bwas[,paste0("cluster_", TPorFP, "_", truePositiveLabel)]),"significant vertices in ", TPorFP, " cluster including ", truePositiveLabel, "located in", bwas$LABELNAME[which(bwas$Probe==truePositiveLabel)] ))
bwas$voxelStatus=NULL

}

return(bwas)
}


############################################
# Extract cluster FWER metrics results
clusterFWER=function(path, bwas, tpList, NNN, moda, iter, signifThreshold, KNearestNeighbours){

library(plyr)
bwas11=NULL
twoTPinCLuster=fwerProblemIncreaseK=0

bwas$AssociatedVertices=ifelse(bwas$Probe %in% tpList$V1, 1 ,0)
bwas$signifVoxel=ifelse(bwas$p < signifThreshold, 1 , 0)
bwas$inTPcluster=0

# Depending on modality - loop on hemisphere or subcortical structures
if (moda %in% c("area", "thickness")){ ROIlist=c("lh", "rh")}
if (moda %in% c("LogJacs", "thick")){ ROIlist=c(10, 11, 12, 13, 17, 18, 26, 49, 50, 51, 52, 53, 54, 58)}

# Loop on the different structures and concatenate results
  for (ROI in ROIlist){

  print(ROI)
  if (moda %in% c("LogJacs", "thick")){ bwas1=bwas[which(bwas$LABEL==ROI),] }
  if (moda %in% c("area", "thickness")){ bwas1=bwas[which(bwas$hemi==ROI),] }

  # Find True Positive clusters
  for (iii in tpList$V1[which(tpList$V1 %in% bwas1$Probe)]){
   #  check TP vertex has not been attributed to a cluster already
    if (bwas1$inTPcluster[which(bwas1$Probe==iii)]==0){
  bwas1=significantClusterAroundVertex(bwas=bwas1, signifThreshold=signifThreshold, truePositiveLabel = iii, TPorFP = "TP", nbNearestNeighbours = KNearestNeighbours)
  bwas1$inTPcluster=rowSums(bwas1[,grep(x = colnames(bwas1), pattern = "cluster_TP"), drop=FALSE])
    }
  }

  # Find False Positive clusters
  bwas1$inFPcluster=0
  while (length(which(bwas1$inTPcluster==0 & bwas1$signifVoxel==1 & bwas1$inFPcluster==0))>0){
  FPvertex=bwas1$Probe[which(bwas1$inTPcluster==0 & bwas1$signifVoxel==1  & bwas1$inFPcluster==0)][1]
  bwas1=significantClusterAroundVertex(bwas=bwas1, signifThreshold=signifThreshold, truePositiveLabel = FPvertex, TPorFP = "FP", nbNearestNeighbours = KNearestNeighbours)
  bwas1$inFPcluster=rowSums(bwas1[,grep(x = colnames(bwas1), pattern = "cluster_FP"), drop=FALSE])
  bwas1$inFPandTPcluster=rowSums( cbind(ifelse(bwas1$inFPcluster>0, 1, 0), ifelse(bwas1$inTPcluster>0, 1, 0)))

  if (length(table(bwas1$inFPcluster))>2){
 print("WARNING: some vertices have been attributed to several FP clusters. Check the results! To solve problem, you may need to modify the nearest neighbour search option k, in mcNNindex() (used in function significantClusterAroundVertex )")  }

    if (length(table(bwas1$inFPandTPcluster))>2){
 print("WARNING: false positive cluster overlaps true positive one ! The 2 clusters are merged into the true positive to avoid inflating FWER. Flag fwerProblemIncreaseK is set to 1. To avoid this from happening too often, increase the number of nearest neigbours")

    # Merge clusters and resets variable
    TPvertex=bwas1$Probe[which(bwas1$inFPandTPcluster==2 & bwas1$AssociatedVertices==1)][1]
    bwas1[which(bwas1[,paste0("cluster_FP_", FPvertex)]==1),paste0("cluster_TP_", TPvertex)]=1
    bwas1[,paste0("cluster_FP_", FPvertex)]=NULL
    fwerProblemIncreaseK=1
    bwas1$inTPcluster=rowSums(bwas1[,grep(x = colnames(bwas1), pattern = "cluster_TP"), drop=FALSE])
    bwas1$inFPcluster=rowSums(bwas1[,grep(x = colnames(bwas1), pattern = "cluster_FP"), drop=FALSE])
    bwas1$inFPandTPcluster=rowSums( cbind(ifelse(bwas1$inFPcluster>0, 1, 0), ifelse(bwas1$inTPcluster>0, 1, 0)))
      }
  }
    # Combine all ROI - hemisphere of subcortical volumes
    bwas11=rbind.fill(bwas11, bwas1)
  }
    # Extract statistics and information on the full data
    FWERcluster=ifelse( length(grep(x = colnames(bwas11), pattern = "cluster_FP"))>0, 1 , 0)
    nbFPclusters=length(grep(x = colnames(bwas11), pattern = "cluster_FP"))
    nbTPclusters=length(grep(x = colnames(bwas11), pattern = "cluster_TP"))
    maxSizeFPclusters=max(colSums(bwas11[,grep(x = colnames(bwas11), pattern = "cluster_FP"), drop=FALSE], na.rm = T), na.rm = T)
    medianSizeFPclusters=median(colSums(bwas11[,grep(x = colnames(bwas11), pattern = "cluster_FP"), drop=FALSE], na.rm = T), na.rm = T)
    minSizeFPclusters=min(colSums(bwas11[,grep(x = colnames(bwas11), pattern = "cluster_FP"), drop=FALSE], na.rm = T), na.rm = T)
    maxSizeTPcluster=max(colSums(bwas11[,grep(x = colnames(bwas11), pattern = "cluster_TP"), drop=FALSE], na.rm = T), na.rm = T)
    medianSizeTPcluster=median(colSums(bwas11[,grep(x = colnames(bwas11), pattern = "cluster_TP"), drop=FALSE], na.rm = T), na.rm = T)
    minSizeTPcluster=min(colSums(bwas11[,grep(x = colnames(bwas11), pattern = "cluster_TP"), drop=FALSE], na.rm = T), na.rm = T)
    isFPclusterLargerTPcluster=ifelse(maxSizeFPclusters>minSizeTPcluster, 1 ,0)

    FWERcluster_ROI=0
    for (iii in grep(x = colnames(bwas11), pattern = "cluster_TP")){
      if( length(table(bwas11$LABELNAME[which(bwas11[,iii]==1)]))>1 ){ FWERcluster_ROI=1   }   }

     # Format results
    res=as.data.frame(t(c(FWERcluster, nbFPclusters, maxSizeFPclusters,medianSizeFPclusters, minSizeFPclusters, isFPclusterLargerTPcluster,nbTPclusters, maxSizeTPcluster,medianSizeTPcluster,minSizeTPcluster, FWERcluster_ROI, twoTPinCLuster, fwerProblemIncreaseK)))
    colnames(res)=c("FWERcluster", "NumberFalsePositiveClusters", "maxSizeFalsePositiveCluster", "medianSizeFalsePositiveCluster", "minSizeFalsePositiveCluster" ,"IsThereAFPCLusterLargerThanTPCluster","nbTPclusters","maxSizeTruePositiveCluster", "MedianSizeTruePositiveCluster","minSizeTPcluster", "FWERcluster_ROI", "twoTPinCLuster", "fwerProblemIncreaseK")
    # Write outputs and bwas table with all the informations (for plots etc..)

    write.table(bwas11, file = paste0(path, "/BWAS_NVertices", NNN, "_", moda, "_Niter", iter, "_clustersAndCoordinates"), quote = F )
    return(res)
}



######
clusterFWER_fwhm20=function(path, bwas, tpList, NNN, moda, iter, signifThreshold, KNearestNeighbours){

library(plyr)
bwas11=NULL
twoTPinCLuster=fwerProblemIncreaseK=0

bwas$AssociatedVertices=ifelse(bwas$Probe %in% tpList$V1, 1 ,0)
bwas$signifVoxel=ifelse(bwas$p < signifThreshold, 1 , 0)
bwas$inTPcluster=0

# Depending on modality - loop on hemisphere or subcortical structures
if (moda %in% c("area", "thickness")){ ROIlist=c("lh", "rh")}
if (moda %in% c("LogJacs", "thick")){ ROIlist=c(10, 11, 12, 13, 17, 18, 26, 49, 50, 51, 52, 53, 54, 58)}

# Loop on the different structures and concatenate results
  for (ROI in ROIlist){

  print(ROI)
  if (moda %in% c("LogJacs", "thick")){ bwas1=bwas[which(bwas$LABEL==ROI),] }
  if (moda %in% c("area", "thickness")){ bwas1=bwas[which(bwas$hemi==ROI),] }

  # Find True Positive clusters
  for (iii in tpList$V1[which(tpList$V1 %in% bwas1$Probe)]){
   #  check TP vertex has not been attributed to a cluster already
    if (bwas1$inTPcluster[which(bwas1$Probe==iii)]==0){
  bwas1=significantClusterAroundVertex(bwas=bwas1, signifThreshold=signifThreshold, truePositiveLabel = iii, TPorFP = "TP", nbNearestNeighbours = KNearestNeighbours)
  bwas1$inTPcluster=rowSums(bwas1[,grep(x = colnames(bwas1), pattern = "cluster_TP"), drop=FALSE])
    }
  }

  # Find False Positive clusters
  bwas1$inFPcluster=0
  while (length(which(bwas1$inTPcluster==0 & bwas1$signifVoxel==1 & bwas1$inFPcluster==0))>0){
  FPvertex=bwas1$Probe[which(bwas1$inTPcluster==0 & bwas1$signifVoxel==1  & bwas1$inFPcluster==0)][1]
  bwas1=significantClusterAroundVertex(bwas=bwas1, signifThreshold=signifThreshold, truePositiveLabel = FPvertex, TPorFP = "FP", nbNearestNeighbours = KNearestNeighbours)
  bwas1$inFPcluster=rowSums(bwas1[,grep(x = colnames(bwas1), pattern = "cluster_FP"), drop=FALSE])
  bwas1$inFPandTPcluster=rowSums( cbind(ifelse(bwas1$inFPcluster>0, 1, 0), ifelse(bwas1$inTPcluster>0, 1, 0)))

  if (length(table(bwas1$inFPcluster))>2){
 print("WARNING: some vertices have been attributed to several FP clusters. Check the results! To solve problem, you may need to modify the nearest neighbour search option k, in mcNNindex() (used in function significantClusterAroundVertex )")  }

    if (length(table(bwas1$inFPandTPcluster))>2){
 print("WARNING: false positive cluster overlaps true positive one ! The 2 clusters are merged into the true positive to avoid inflating FWER. Flag fwerProblemIncreaseK is set to 1. To avoid this from happening too often, increase the number of nearest neigbours")

    # Merge clusters and resets variable
    TPvertex=bwas1$Probe[which(bwas1$inFPandTPcluster==2 & bwas1$AssociatedVertices==1)][1]
    bwas1[which(bwas1[,paste0("cluster_FP_", FPvertex)]==1),paste0("cluster_TP_", TPvertex)]=1
    bwas1[,paste0("cluster_FP_", FPvertex)]=NULL
    fwerProblemIncreaseK=1
    bwas1$inTPcluster=rowSums(bwas1[,grep(x = colnames(bwas1), pattern = "cluster_TP"), drop=FALSE])
    bwas1$inFPcluster=rowSums(bwas1[,grep(x = colnames(bwas1), pattern = "cluster_FP"), drop=FALSE])
    bwas1$inFPandTPcluster=rowSums( cbind(ifelse(bwas1$inFPcluster>0, 1, 0), ifelse(bwas1$inTPcluster>0, 1, 0)))
      }
  }
    # Combine all ROI - hemisphere of subcortical volumes
    bwas11=rbind.fill(bwas11, bwas1)
  }
    # Extract statistics and information on the full data
    FWERcluster=ifelse( length(grep(x = colnames(bwas11), pattern = "cluster_FP"))>0, 1 , 0)
    nbFPclusters=length(grep(x = colnames(bwas11), pattern = "cluster_FP"))
    nbTPclusters=length(grep(x = colnames(bwas11), pattern = "cluster_TP"))
    maxSizeFPclusters=max(colSums(bwas11[,grep(x = colnames(bwas11), pattern = "cluster_FP"), drop=FALSE], na.rm = T), na.rm = T)
    medianSizeFPclusters=median(colSums(bwas11[,grep(x = colnames(bwas11), pattern = "cluster_FP"), drop=FALSE], na.rm = T), na.rm = T)
    minSizeFPclusters=min(colSums(bwas11[,grep(x = colnames(bwas11), pattern = "cluster_FP"), drop=FALSE], na.rm = T), na.rm = T)
    maxSizeTPcluster=max(colSums(bwas11[,grep(x = colnames(bwas11), pattern = "cluster_TP"), drop=FALSE], na.rm = T), na.rm = T)
    medianSizeTPcluster=median(colSums(bwas11[,grep(x = colnames(bwas11), pattern = "cluster_TP"), drop=FALSE], na.rm = T), na.rm = T)
    minSizeTPcluster=min(colSums(bwas11[,grep(x = colnames(bwas11), pattern = "cluster_TP"), drop=FALSE], na.rm = T), na.rm = T)
    isFPclusterLargerTPcluster=ifelse(maxSizeFPclusters>minSizeTPcluster, 1 ,0)

    FWERcluster_ROI=0
    for (iii in grep(x = colnames(bwas11), pattern = "cluster_TP")){
      if( length(table(bwas11$LABELNAME[which(bwas11[,iii]==1)]))>1 ){ FWERcluster_ROI=1   }   }

     # Format results
    res=as.data.frame(t(c(FWERcluster, nbFPclusters, maxSizeFPclusters,medianSizeFPclusters, minSizeFPclusters, isFPclusterLargerTPcluster,nbTPclusters, maxSizeTPcluster,medianSizeTPcluster,minSizeTPcluster, FWERcluster_ROI, twoTPinCLuster, fwerProblemIncreaseK)))
    colnames(res)=c("FWERcluster", "NumberFalsePositiveClusters", "maxSizeFalsePositiveCluster", "medianSizeFalsePositiveCluster", "minSizeFalsePositiveCluster" ,"IsThereAFPCLusterLargerThanTPCluster","nbTPclusters","maxSizeTruePositiveCluster", "MedianSizeTruePositiveCluster","minSizeTPcluster", "FWERcluster_ROI", "twoTPinCLuster", "fwerProblemIncreaseK")
    # Write outputs and bwas table with all the informations (for plots etc..)

    write.table(bwas11, file = paste0(path, "/BWAS_NVertices", NNN, "_", moda, "_Niter", iter, "_clustersAndCoordinates_fwhm20"), quote = F )
    return(res)
}




#####################################
# PLOT cluster analysis

PlotClustersTP_FP_cortical=function(path, NNN, hemi, moda, iter){
library(rgl)
library(viridis)
cols=viridis(5, option = "C")

bwasPlot=read.table(paste0(path, "/BWAS_NVertices", NNN, "_", moda, "_Niter", iter, "_clustersAndCoordinates"), header=T)
bwasPlot=bwasPlot[which(bwasPlot$hemi==hemi),]

bwasPlot$color="lightgrey"
bwasPlot$color[which(bwasPlot$inTPcluster>0)]=cols[5]
bwasPlot$color[which(bwasPlot$inFPcluster>0)]=cols[3]
bwasPlot$color[which(bwasPlot$inTPcluster>1)]=cols[2]
bwasPlot$color[which(bwasPlot$inFPandTPcluster>1)]=cols[4]
bwasPlot$color[which(bwasPlot$AssociatedVertices==1)]=cols[1]

bwasPlot$radius=0.5
bwasPlot$radius[which(bwasPlot$inTPcluster>0)]=0.8
bwasPlot$radius[which(bwasPlot$inFPcluster>0)]=0.8
bwasPlot$radius[which(bwasPlot$AssociatedVertices==1)]=1.5
bwasPlot$radius[which(bwasPlot$AssociatedVertices==1 & bwasPlot$signifVoxel==0)]=0.9
#bwasPlot$X=bwasPlot$X*(-1)

# Downsample
bwasPlot2=bwasPlot[sample(1:length(bwasPlot$ProbeID), size = 10000),]
bwasPlot2=rbind(bwasPlot2, bwasPlot[which(bwasPlot$AssociatedVertices==1 | bwasPlot$inTPcluster>0 | bwasPlot$inFPcluster>0),])

# Draw plot
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot2[,c( "Y","X", "Z")]), col=bwasPlot2$color, radius = bwasPlot2$radius, alpha=bwasPlot2$alpha)
rgl.snapshot(paste0(path, "/BWAS_NVertices", NNN, "_", hemi, "_", moda, "_Niter", iter, "_clustersAndCoordinates_inside.linear.png"))
rgl.close()

bwasPlot2$X=bwasPlot2$X*(-1)
bwasPlot2$Y=bwasPlot2$Y*(-1)
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot2[,c( "Y","X", "Z")]), col=bwasPlot2$color, radius = bwasPlot2$radius, alpha=bwasPlot$alpha)
rgl.snapshot(paste0(path, "/BWAS_NVertices", NNN, "_", hemi, "_", moda, "_Niter", iter, "_clustersAndCoordinates_outside.linear.png"))
rgl.close()

}




######################
# ADD subcortical coordinates

bwasAddCoordinatesSubcortical=function(bwasFormatted, hemi, moda){
# Get vertices coordinates
coor=read.table(paste0(atlasPath,"/Atlas_coordinates_subortical.txt"), stringsAsFactors = F, header=T)
# Add coordinates
bwasFormatted=merge(bwasFormatted, coor, by="ProbeID")
return(bwasFormatted)
}



########################
# Plot subcortical

PlotClustersTP_FP_subcortical=function(path, NNN, hemi, moda, iter){

library(rgl)
library(viridis)
cols=viridis(5, option = "C")

bwasPlot=read.table(paste0(path, "/BWAS_NVertices", NNN, "_", moda, "_Niter", iter, "_clustersAndCoordinates"), header=T)
bwasPlot=bwasPlot[which(bwasPlot$hemi==hemi),]

bwasPlot$color="lightgrey"
bwasPlot$color[which(bwasPlot$inTPcluster>0)]=cols[5]
bwasPlot$color[which(bwasPlot$inFPcluster>0)]=cols[3]
bwasPlot$color[which(bwasPlot$inTPcluster>1)]=cols[2]
bwasPlot$color[which(bwasPlot$inFPandTPcluster>1)]=cols[4]
bwasPlot$color[which(bwasPlot$AssociatedVertices==1)]=cols[1]

bwasPlot$radius=0.5
bwasPlot$radius[which(bwasPlot$inTPcluster>0)]=0.6
bwasPlot$radius[which(bwasPlot$inFPcluster>0)]=0.6
bwasPlot$radius[which(bwasPlot$AssociatedVertices==1)]=1.5
bwasPlot$radius[which(bwasPlot$AssociatedVertices==1 & bwasPlot$signifVoxel==0)]=0.9
bwasPlot$X=bwasPlot$X*(-1)
bwasPlot$Y=bwasPlot$Y*(-1)

# Downsample
bwasPlot2=bwasPlot[sample(1:length(bwasPlot$ProbeID), size = 7000),]
bwasPlot2=rbind(bwasPlot2, bwasPlot[which(bwasPlot$AssociatedVertices==1 | bwasPlot$inTPcluster>0 | bwasPlot$inFPcluster>0),])


# Draw plot
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot2[,c( "Z",  "X", "Y")]), col=bwasPlot2$color, radius = bwasPlot2$radius)
rgl.snapshot(paste0(path, "/BWAS_NVertices", NNN, "_", hemi, "_", moda, "_Niter", iter, "_clustersAndCoordinates_outside.linear.png"))
rgl.close()

bwasPlot2$X=bwasPlot2$X*(-1)
bwasPlot2$Z=bwasPlot2$Z*(-1)
par3d(windowRect = c(0, 0, 800, 800)*1.5, zoom=0.8)
spheres3d(as.matrix(bwasPlot2[,c( "Z",  "X", "Y")]), col=bwasPlot2$color, radius = bwasPlot2$radius)
rgl.snapshot(paste0(path, "/BWAS_NVertices", NNN, "_", hemi, "_", moda, "_Niter", iter, "_clustersAndCoordinates_inside.linear.png"))
rgl.close()

}


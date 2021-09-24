# Functions that opens BWAS sum-stats, to create BRS and evaluate them

makeOSCAscorefile=function(scenario, iter, NNN){

vert=vertSig=NULL
for (moda in c("area", "thickness", "LogJacs", "thick")){

# Open sum stats with cluster identified
bwas=read.table(paste0(scenario, "_clusterFWER/BWAS_NVertices", NNN, "_", moda,"_Niter", iter, "_clustersAndCoordinates"), stringsAsFactors = F)

# Top vertex per cluster
clusVars=NULL
clusVars=colnames(bwas)[grep("cluster_", x = colnames(bwas))]
if (length(clusVars)>0){
for (iii in 1:length(clusVars)){
cluster=bwas[which(bwas[,clusVars[iii]]==1),]
cluster=cluster[order(cluster$p),]
vert=rbind(vert, c(cluster$Probe[1], cluster$b[1]) )
}}

# All signif vertices
vertSig=rbind(vertSig, bwas[which(bwas$signifVoxel==1),c("Probe", "b")] )
}

# Export tables for OSCA score calculation
write.table(vert, paste0( scenario, "_clusterFWER/BWAS_topVerticesCluster_NVertices", NNN, "_","Niter", iter,  ".txt"), col.names = F, row.names = F , quote=F)
write.table(vertSig, paste0( scenario, "_clusterFWER/BWAS_allSignifVertices_NVertices", NNN, "_", "Niter", iter,  ".txt"), col.names = F, row.names = F , quote=F)

}


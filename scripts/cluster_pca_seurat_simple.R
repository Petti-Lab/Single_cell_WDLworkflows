#!/usr/local/bin/Rscript

.libPaths( c("/storage1/fs1/allegra.petti/Active/R_libs_scratch/RLibs_4.0.3",.libPaths()) )
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(Matrix)
library(RColorBrewer)
library(ggthemes)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 5) {
  args <- c("--help")
}

Seurat_rds <- as.character(args[1])
sample_name <- as.character(args[2])
cluster.res =  as.numeric(args[3])
nPC= as.integer(args[4])
CellCycleTirosh =  as.character(args[5])

scrna_GEX <- readRDS(Seurat_rds)

genome <- "GRCh38";
date = gsub("2022-","22",Sys.Date(),perl=TRUE);
date = gsub("-","",date);

control = "Cycling"

print("Performing cell cycle analysis...")
cell.cycle.tirosh <- read.table(CellCycleTirosh, sep='\t', header=FALSE);
s.genes = cell.cycle.tirosh$V2[which(cell.cycle.tirosh$V1 == "G1/S")];
g2m.genes = cell.cycle.tirosh$V2[which(cell.cycle.tirosh$V1 == "G2/M")];
scrna_GEX <- CellCycleScoring(object=scrna_GEX, s.features=s.genes, g2m.features=g2m.genes, set.ident=FALSE)
##
##
### Normalize the data
#print("Normalizing data...")
scrna_GEX <- NormalizeData(object = scrna_GEX, normalization.method = "LogNormalize", scale.factor = 1e6); # 1e6 is new as of 1/8/20
#
## Detection of variable genes across the single cells for downstream analysis see HVFInfo(object=scrna_GEX[["RNA"]])
#
scrna_GEX <- FindVariableFeatures(object = scrna_GEX, selection.method = 'vst', mean.cutoff = c(0.1,8), dispersion.cutoff = c(1, Inf))
print(paste("Number of Variable Features: ",length(x = VariableFeatures(object = scrna_GEX))));
#
VG.file = sprintf("%s.%s.variableGenes.%s.pdf", sample_name, control, date);
pdf(VG.file, useDingbats=FALSE)
vg <- VariableFeaturePlot(scrna_GEX)
print(vg);
dev.off()
###
#
if (control == "Cycling") { # This removes all signal associated with the cell cycle
  scrna_GEX <- ScaleData(object = scrna_GEX, features = rownames(x = scrna_GEX), vars.to.regress = c("S.Score","G2M.Score"), display.progress=FALSE);
} else if (control == "CyclingRB") {
  scrna_GEX <- ScaleData(object = scrna_GEX, features = rownames(x = scrna_GEX), vars.to.regress = c("S.Score","G2M.Score","percent.ribo"), display.progress=FALSE);
} else if (control == "CyclingDiff") {
  scrna_GEX <- ScaleData(object = scrna_GEX, features = rownames(x = scrna_GEX), vars.to.regress = c("CC.Difference"), display.progress=FALSE);
} else {
  scrna_GEX <- ScaleData(object = scrna_GEX, features = rownames(x = scrna_GEX), display.progress=FALSE);
}
#
##
print("Saving scrna_GEX object after scTransform...");
saveRDS(scrna_GEX, file = sprintf("%s.%s.VST.%s.rds", sample_name,control,date))
print("Object saved.")

# get_significant_pcs <- function(scrna_GEX) {
#   control='Cycling'
# #   if(length(unique(scrna_GEX$orig.ident))==1){
# #     nPC=20 
# #   }else{
# #     nPC=50
# #   }
#   # print('Does the error happen inside get_significant_pcs?')
#   scrna_GEX <- RunPCA(scrna_GEX, npcs = nPC, verbose = FALSE)
#   scrna_GEX <- JackStraw(object = scrna_GEX, num.replicate = 100, dims=nPC)
#   scrna_GEX <- ScoreJackStraw(object = scrna_GEX, dims = 1:nPC)
#   jpeg(sprintf("PCA.jackstraw.%s.%s.jpg", control, date), width = 10, height = 6, units="in", res=300);
#   js <- JackStrawPlot(object = scrna_GEX, dims = 1:nPC)
#   print(js);
#   dev.off();
#   pc.pval <- scrna_GEX@reductions$pca@jackstraw@overall.p.values
#   print(pc.pval);
#   nPC=length( pc.pval[,'Score'][pc.pval[,'Score'] <= 0.001]) #using more stringent criteria of p-value <= 0.001
#   # print('No the error does not happen inside get_significant_pcs!')
#   #redefine nPCs based on number of significant prinicipal components in jackstraw plot
#   return(nPC)
# } #pipeline might break downstream as it expects same no of pcs and cluster resolution downstream (which needs to be changed to be more flexible)

# nPC = get_significant_pcs(scrna_GEX = scrna_GEX)

scrna_GEX <- RunPCA(object = scrna_GEX, npcs = nPC, verbose = FALSE);
print("Saving scrna_GEX object...");
saveRDS(scrna_GEX, file = sprintf("%s.%s.PCA.%s.rds", sample_name,control,date));
print("Object saved.");
# print ("Elbow Running...");
# elbow <- ElbowPlot(object = scrna_GEX)
# jpeg(sprintf("%s.PCA30.elbow.%s.%s.jpg", sample_name, control, date), width = 6, height = 8, units="in", res=300);
# print(elbow);
# dev.off();
scrna_GEX <- RunUMAP(object = scrna_GEX, reduction = "pca", dims = 1:nPC)
scrna_GEX <- RunTSNE(object = scrna_GEX, reduction = "pca", dims = 1:nPC)
##
print("Saving scrna_GEX object...");
saveRDS(scrna_GEX, file = sprintf("%s.%s.SCT.PCA.UMAP.TSNE.%s.rds", sample_name,control,date))
print("Object saved.")
#
print(sprintf("Plot default UMAP & tSNE with %d  PCs...",nPC))
jpeg(sprintf("%s.UMAP.%s.%s.jpg", sample_name, control, date), width = 10, height = 8, units="in", res=300);
#p1 <- DimPlot(object = scrna_GEX, reduction = "tsne", group.by = "Batch", pt.size=0.1)
p2 <- DimPlot(object = scrna_GEX, reduction = "umap", group.by = "Sample", pt.size=0.1)
#print(plot_grid(p1, p2));
print(p2);
dev.off();
##
print ("VizDimLoadings Running...");
jpeg(sprintf("%s.VizDimLoadings.%s.%s.jpg", sample_name, control, date), width = 8, height = 30, units="in", res=300);
vdl <- VizDimLoadings(object = scrna_GEX, dims = 1:3)
print(vdl);
dev.off();
##
### ProjectDim scores each gene in the dataset (including genes not included in the PCA) based on their correlation
### with the calculated components. Though we don't use this further here, it can be used to identify markers that
### are strongly correlated with cellular heterogeneity, but may not have passed through variable gene selection.
### The results of the projected PCA can be explored by setting use.full=T in the functions above
print ("ProjectDim Running...");
scrna_GEX <- ProjectDim(object = scrna_GEX)
#
#### In particular `DimHeatmap` allows for easy exploration of the primary sources of heterogeneity in a dataset, and can be useful when trying to decide which PCs to include for further downstream analyses. Both cells and genes are ordered according to their PCA scores. Setting cells.use to a number plots the 'extreme' cells on both ends of the spectrum, which dramatically speeds plotting for large datasets. Though clearly a supervised analysis, we find this to be a valuable tool for exploring correlated gene sets.
#
print ("DimHeatmap Running...");
jpeg(sprintf("%s.PCA.heatmap.top.%s.%s.jpg", sample_name, control, date), width = 8.5, height = 11, units="in", res=300);
hm <- DimHeatmap(object = scrna_GEX, dims = 1, cells = 500, balanced = TRUE);
print(hm);
dev.off();

jpeg(sprintf("%s.PCA.heatmap.multi.%s.%s.jpg", sample_name, control, date), width = 8.5, height = 24, units="in", res=300);
hm.multi <- DimHeatmap(object = scrna_GEX, dims = 1:10, cells = 500, balanced = TRUE);
print(hm.multi);
dev.off();
##
##### A more ad hoc method for determining which PCs to use is to look at a plot of the standard deviations of the principle components and draw your cutoff where there is a clear elbow in the graph. 
##
print ("Elbow Running...");
elbow <- ElbowPlot(object = scrna_GEX)
jpeg(sprintf("%s.PCA.elbow.%s.%s.jpg", sample_name, control, date), width = 6, height = 8, units="in", res=300);
print(elbow);
dev.off();
#
print("Performing JackStraw analysis...")
scrna_GEX <- JackStraw(object = scrna_GEX, num.replicate = 100, dims=nPC)
scrna_GEX <- ScoreJackStraw(object = scrna_GEX, dims = 1:nPC)
jpeg(sprintf("%s.PCA.jackstraw.%s.%s.jpg", sample_name, control, date), width = 10, height = 6, units="in", res=300);
js <- JackStrawPlot(object = scrna_GEX, dims = 1:nPC)
print(js);
dev.off();
##
## print overall pvalues for each PC:
pc.pval <- scrna_GEX@reductions$pca@jackstraw@overall.p.values
print(pc.pval);
write.table(pc.pval, file=sprintf("%s.PCA.jackstraw.scores.%s.%s.xls", sample_name, control, date), quote=FALSE, sep='\t', col.names=TRUE);
#
print("Saving scrna_GEX object...");
# save everything up to here:
saveRDS(scrna_GEX, file = sprintf("%s.%s.SCT.PCA.UMAP.TSNE.%s.rds", sample_name, control, date))
print("Object saved.")
#
scrna_GEX <- FindNeighbors(object=scrna_GEX, dims=1:nPC);
scrna_GEX <- FindClusters(object=scrna_GEX, resolution=cluster.res);
scrna_GEX <- StashIdent(object = scrna_GEX, save.name = sprintf("ClusterNames_%.1f_%dPC", cluster.res, nPC))
DEGs <- FindAllMarkers(object=scrna_GEX); # output is a matrix!
write.table(DEGs, file=sprintf("%s.DEGs.Wilcox.cluster.res.%.1f.PCA.%d.%s.%s.xls", sample_name,cluster.res,nPC, control, date), quote=FALSE, sep="\t", row.names=FALSE) # must save cluster-specific marker genes
#
print("Saving scrna_GEX object...");
saveRDS(scrna_GEX, file = sprintf("%s.%s.SCT.PCA.UMAP.TSNE.CLUST.%s.rds", sample_name, control, date))
#print("Object saved.") 
#scrna_GEX <- readRDS('/scratch1/fs1/allegra.petti/khan.saad/MouseTIL_analysis/MouseTIL_Seurat/Cycling.SCT.PCA.UMAP.TSNE.CLUST.20200527.rds')
# scrna_GEX <- readRDS('/scratch1/fs1/allegra.petti/khan.saad/MouseTIL_analysis/MouseTIL_Seurat_mc0.05/Cycling.SCT.PCA.UMAP.TSNE.CLUST.20200529.rds')
print("Plot graph-based clusters on UMAP...");
control <- 'Cycling'
n.graph = length(unique(scrna_GEX[[sprintf("ClusterNames_%.1f_%dPC",cluster.res, nPC)]][,1]));
rainbow.colors = rainbow(n.graph, s=0.6, v=0.9);

jpeg(sprintf("%s.UMAP.clusters.%d.%.1f.%s.%s.jpg", sample_name, nPC, cluster.res, control, date), width = 10, height = 8, units="in", res=300);
p2 <- DimPlot(object = scrna_GEX, reduction = "umap", group.by = sprintf("ClusterNames_%.1f_%dPC",cluster.res, nPC), cols = rainbow.colors, pt.size=0.1) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
print(p2);
dev.off();

jpeg(sprintf("%s.UMAP.clusters.labeled.%d.%.1f.%s.%s.jpg", sample_name, nPC, cluster.res, control, date), width = 10, height = 8, units="in", res=300);
p2 <- DimPlot(object = scrna_GEX, reduction = "umap", group.by = sprintf("ClusterNames_%.1f_%dPC",cluster.res, nPC), cols = rainbow.colors, pt.size=0.1, label=TRUE,label.size = 5) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
print(p2);
dev.off();

mito.genes <- grep(pattern = "^MT-", x = rownames(x = scrna_GEX), value = TRUE,ignore.case = TRUE);
percent.mito <- Matrix::colSums(x = GetAssayData(object = scrna_GEX, slot = 'counts')[mito.genes, ]) / Matrix::colSums(x = GetAssayData(object = scrna_GEX, slot = 'counts'));
scrna_GEX[['percent.mito']] <- percent.mito;

# ribosomal genes
ribo.genes <- grep(pattern = "^RP[SL][[:digit:]]", x = rownames(x = scrna_GEX), value = TRUE,ignore.case = TRUE);
percent.ribo <- Matrix::colSums(x = GetAssayData(object = scrna_GEX, slot = 'counts')[ribo.genes, ]) / Matrix::colSums(x = GetAssayData(object = scrna_GEX, slot = 'counts'));

scrna_GEX[['percent.ribo']] <- percent.ribo;

feature.pal = rev(colorRampPalette(brewer.pal(11,"Spectral"))(50));

print("Making additional UMAP plots");
# color UMAP plots by parameters of interest:
print ("color by UMI");
jpeg(sprintf("%s.umap.%d.%.1f.colorby.UMI.%s.%s.jpg", sample_name, nPC, cluster.res, control, date), width = 10, height = 8,  units="in", res=300);
fp2 <- FeaturePlot(object = scrna_GEX, features = c("nCount_RNA"), cols = feature.pal, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
print(fp2);
dev.off();

print ("color by % mito");
jpeg(sprintf("%s.umap.%d.%.1f.colorby.MC.%s.%s.jpg", sample_name, nPC, cluster.res, control, date), width = 10, height = 8, units="in", res=300);
fp2 <- FeaturePlot(object = scrna_GEX, features = c("percent.mito"), cols = feature.pal, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
print(fp2);
dev.off();

print ("color by % RB");
jpeg(sprintf("%s.umap.%d.%.1f.colorby.RB.%s.%s.jpg", sample_name, nPC, cluster.res, control, date), width = 10, height = 8, units="in", res=300);
fp2 <- FeaturePlot(object = scrna_GEX, features = c("percent.ribo"), cols = feature.pal, pt.size=0.1, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
print(fp2);
dev.off();

print ("color by nGene");
jpeg(sprintf("%s.umap.%d.%.1f.colorby.nGene.%s.%s.jpg", sample_name, nPC, cluster.res, control, date), width = 10, height = 8, units="in", res=300);
fp2 <- FeaturePlot(object = scrna_GEX, features = c("nFeature_RNA"), cols = feature.pal, pt.size=0.1, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
print(fp2);
dev.off();

print ("color by Phase");
jpeg(sprintf("%s.umap.%d.%.1f.colorby.Phase.%s.%s.jpg", sample_name, nPC, cluster.res, control, date), width = 10, height = 8, units="in", res=300);
phase.colors = ptol_pal()(3)
umapplot <- DimPlot(object = scrna_GEX, cols=phase.colors, group.by="Phase", pt.size=0.1, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank())
print(umapplot);
dev.off();

print ("color UMAP by Principal Components");
jpeg(sprintf("%s.UMAP.%d.%.1f.colorby.PCs.%s.%s.jpg", sample_name, nPC, cluster.res, control, date), width = 12, height = 6, units="in", res=100);
redblue=c("blue","gray","red");
fp1 <- FeaturePlot(object = scrna_GEX, features = 'PC_1', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
fp2 <- FeaturePlot(object = scrna_GEX, features = 'PC_2', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
fp3 <- FeaturePlot(object = scrna_GEX, features = 'PC_3', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
fp4 <- FeaturePlot(object = scrna_GEX, features = 'PC_4', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
fp5 <- FeaturePlot(object = scrna_GEX, features = 'PC_5', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
fp6 <- FeaturePlot(object = scrna_GEX, features = 'PC_6', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
fp7 <- FeaturePlot(object = scrna_GEX, features = 'PC_7', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
fp8 <- FeaturePlot(object = scrna_GEX, features = 'PC_8', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
fp9 <- FeaturePlot(object = scrna_GEX, features = 'PC_9', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
fp10 <- FeaturePlot(object = scrna_GEX, features = 'PC_10', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
print(plot_grid(fp1, fp2, fp3, fp4, fp5, fp6, fp7, fp8, fp9, fp10));
# print(plot_grid(fp1, fp2, fp3, fp4, fp5, fp6, fp7, fp8, fp9));
dev.off();

print("Saving...");
# save everything up to here
saveRDS(scrna_GEX, file = sprintf("%s.%s.SCT.PCA.UMAP.TSNE.CLUST.%s.rds", sample_name, control, date))


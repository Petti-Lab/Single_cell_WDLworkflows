#!/usr/local/bin/Rscript

.libPaths( c("/storage1/fs1/allegra.petti/Active/R_libs_scratch/RLibs_4.0.3",.libPaths()) )
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(Matrix)
#source('./scripts/Plot_QC_scrnaseq.R')

plot_qc_metrics <- function(seurat_obj,filename,mt_pat="^MT-") {
seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
seurat_obj$mitoRatio <- PercentageFeatureSet(object = seurat_obj, pattern = mt_pat)
seurat_obj$mitoRatio <- seurat_obj@meta.data$mitoRatio / 100
metadata <- seurat_obj@meta.data
metadata$cells <- rownames(metadata)

metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

Ncells <- metadata %>% 
  ggplot(aes(x=Sample, fill=Sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

cell_dens <- metadata %>% 
  ggplot(aes(color=Sample, x=nUMI, fill= Sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)+ggtitle("Cell density")

genes_per_cell <- metadata %>% 
  ggplot(aes(color=Sample, x=nGene, fill= Sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300) + ggtitle("Genes per cell")

box_genespercell <- metadata %>% 
  ggplot(aes(x=Sample, y=log10(nGene), fill=Sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

umibygenesmito <- metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~Sample) + ggtitle("UMI vs Genes colored by mito ratio")


mitoratio_dens <- metadata %>% 
  ggplot(aes(color=Sample, x=mitoRatio, fill=Sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2) + ggtitle("mitoratio density")

genes_per_umi <- metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = Sample, fill=Sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8) + ggtitle("Genes per UMI")

pdf(filename)
print(Ncells)
print(cell_dens)
print(genes_per_cell)
print(box_genespercell)
print(umibygenesmito)
print(mitoratio_dens)
print(genes_per_umi)
dev.off()
}

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2) {
  args <- c("--help")
}

Seurat_10x_directory <- as.character(args[1])
sample_name <- as.character(args[2])

genome <- "GRCh38";
date = gsub("2023-","23",Sys.Date(),perl=TRUE);
date = gsub("-","",date);

# data.10x = Read10X_h5(Seurat_h5file, use.names = TRUE, unique.features = TRUE)
data.10x  <- Read10X(data.dir = Seurat_10x_directory);
scrna_GEX = CreateSeuratObject(counts = data.10x, min.cells=10, min.features=100, project=sample_name);
scrna_GEX$Sample = sample_name

# ##################################################################################
# # This bit is only for testing purposes for the rscript within wdl would be disabled
# # in the main workflow
set.seed(100)
random_sample_of_cells = sample(Cells(scrna_GEX),2500)
#  randomly select 2500 cells for testing the pipeline.
scrna_GEX <- subset(scrna_GEX,cells=random_sample_of_cells)
###################################################################################


mito.genes <- grep(pattern = "^MT-", x = rownames(x = scrna_GEX), value = TRUE,ignore.case = TRUE);
percent.mito <- Matrix::colSums(x = GetAssayData(object = scrna_GEX, slot = 'counts')[mito.genes, ]) / Matrix::colSums(x = GetAssayData(object = scrna_GEX, slot = 'counts'));
scrna_GEX[['percent.mito']] <- percent.mito;

ribo.genes <- grep(pattern = "^RP[SL][[:digit:]]", x = rownames(x = scrna_GEX), value = TRUE,ignore.case = TRUE);
percent.ribo <- Matrix::colSums(x = GetAssayData(object = scrna_GEX, slot = 'counts')[ribo.genes, ]) / Matrix::colSums(x = GetAssayData(object = scrna_GEX, slot = 'counts'));
scrna_GEX[['percent.ribo']] <- percent.ribo;

print("Making violin plots...");
#
png(sprintf("%s.VlnPlot.%s.png", sample_name, date), width = 13, height = 6, units="in", res=300);
vln <- VlnPlot(object = scrna_GEX, features = c("percent.mito", "percent.ribo"), pt.size=0, ncol = 2, group.by="Sample");
print(vln);
dev.off();

png(sprintf("%s.VlnPlot.nCount.25Kmax.%s.png", sample_name,date), width = 10, height = 10, units="in", res=300)
vln <- VlnPlot(object = scrna_GEX, features = "nCount_RNA", pt.size=0, group.by="Sample", y.max=25000)
print(vln)
dev.off();

png(sprintf("%s.VlnPlot.nFeature.%s.png", sample_name,date), width = 10, height = 10, units="in", res=300)
vln <- VlnPlot(object = scrna_GEX, features = "nFeature_RNA", pt.size=0, group.by="Sample")
print(vln)
dev.off()

print("Making scatter plots...")
pdf(sprintf("%s.Scatter1.%s.pdf", sample_name, date), width = 8, height = 6, useDingbats=FALSE);
scatter <- FeatureScatter(object = scrna_GEX, feature1 = "nCount_RNA", feature2 = "percent.mito", pt.size=0.1,group.by = 'Sample')
print(scatter);
dev.off();
pdf(sprintf("%s.Scatter2.%s.pdf", sample_name, date), width = 8, height = 6, useDingbats=FALSE);
scatter <- FeatureScatter(object = scrna_GEX, feature1 = "nCount_RNA", feature2 = "percent.ribo", pt.size=0.1,group.by = 'Sample')
print(scatter);
dev.off();
pdf(sprintf("%s.Scatter3.%s.pdf", sample_name, date), width = 8, height = 6, useDingbats=FALSE);
scatter <- FeatureScatter(object = scrna_GEX, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.1,group.by = 'Sample')
print(scatter);
dev.off();

saveRDS(scrna_GEX, file = sprintf("%s.unfilteredSeuratObject_GEX.%s.rds",sample_name , date))

plot_qc_metrics(seurat_obj=scrna_GEX,filename=sprintf("%s.Pre_filtering_QC.%s.pdf",sample_name , date))

mc.hi = 0.05;
m <- median(scrna_GEX@meta.data$nFeature_RNA)
s <- sd(scrna_GEX@meta.data$nFeature_RNA)
Feature95 <- quantile(scrna_GEX@meta.data$nFeature_RNA, 0.95);
Feature05 <- quantile(scrna_GEX@meta.data$nFeature_RNA, 0.05);
m1 <- mean(scrna_GEX@meta.data$nCount_RNA)
s1 <- sd(scrna_GEX@meta.data$nCount_RNA)
Count93 <- quantile(scrna_GEX@meta.data$nCount_RNA, 0.93);
scrna_GEX <- subset(x = scrna_GEX, subset = nFeature_RNA > 700 & nCount_RNA < Count93 & percent.mito < mc.hi)

plot_qc_metrics(seurat_obj=scrna_GEX,filename=sprintf("%s.Post_filtering_QC.%s.pdf",sample_name , date))

png(sprintf("%s.VlnPlot.Filtered.nCount.25Kmax.control.%s.png", sample_name, date), width = 10, height = 10, units="in", res=300)
vln <- VlnPlot(object = scrna_GEX, features = "nCount_RNA", pt.size=0, group.by="Sample", y.max=25000)
print(vln)
dev.off();
#
png(sprintf("%s.VlnPlot.Filtered.nFeature.control.%s.png", sample_name,date), width = 10, height = 10, units="in", res=300)
vln <- VlnPlot(object = scrna_GEX, features = "nFeature_RNA", pt.size=0, group.by="Sample")
print(vln)
dev.off()

# get_significant_pcs <- function(scrna_GEX) {
#   control='Cycling'
#   if(length(unique(scrna_GEX$orig.ident))==1){
#     nPC=20 
#   }else{
#     nPC=50
#   }
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
#   nPC=length( pc.pval[,'Score'][pc.pval[,'Score'] <= 0.05]) 
#   # print('No the error does not happen inside get_significant_pcs!')
#   #redefine nPCs based on number of significant prinicipal components in jackstraw plot
#   return(nPC)
# } ##this is used in another script

# the commented out code is used in another script so not using it here

# scrna_GEX <- NormalizeData(object = scrna_GEX, normalization.method = "LogNormalize", scale.factor = 1e6); # 1e6 is new as of 1/8/20

# scrna_GEX <- FindVariableFeatures(object = scrna_GEX, selection.method = 'vst', mean.cutoff = c(0.1,8), dispersion.cutoff = c(1, Inf))
# print(paste("Number of Variable Features: ",length(x = VariableFeatures(object = scrna_GEX))));
# #
# VG.file = sprintf("%s.variableGenes.%s.pdf",control, date);
# pdf(VG.file, useDingbats=FALSE)
# vg <- VariableFeaturePlot(scrna_GEX)
# print(vg);
# dev.off()

# nPC <- get_significant_pcs(scrna_GEX)
# scrna_GEX <- RunPCA(scrna_GEX, npcs = nPC, verbose = FALSE)

# scrna_GEX <- RunUMAP(object = scrna_GEX, reduction = "pca", dims = 1:nPC)
# scrna_GEX <- RunTSNE(object = scrna_GEX, reduction = "pca", dims = 1:nPC)
# scrna_GEX <- JackStraw(object = scrna_GEX, num.replicate = 100, dims=nPC)
# scrna_GEX <- ScoreJackStraw(object = scrna_GEX, dims = 1:nPC)
# jpeg(sprintf("PCA.jackstraw.%s.%s.%s.jpg",sample,control, date), width = 10, height = 6, units="in", res=300);
# js <- JackStrawPlot(object = scrna_GEX, dims = 1:nPC)
# print(js);
# dev.off();

# jpeg(sprintf("UMAP.%s.%s.%s.jpg",sample,control, date), width = 10, height = 8, units="in", res=300);
# p2 <- DimPlot(object = scrna_GEX, reduction = "umap", group.by = "Sample", pt.size=0.1)
# print(p2);
# dev.off();

# print ("VizDimLoadings Running...");
# jpeg(sprintf("VizDimLoadings.%s.%s.%s.jpg",sample,control, date), width = 8, height = 30, units="in", res=300);
# vdl <- VizDimLoadings(object = scrna_GEX, dims = 1:3)
# print(vdl);
# dev.off();

# print ("ProjectDim Running...");
# scrna_GEX <- ProjectDim(object = scrna_GEX)

# # saveRDS(scrna_GEX, file = sprintf("%s.SCT.PCA.UMAP.TSNE.%s.rds",control, date))

# print ("DimHeatmap Running...");
# jpeg(sprintf("PCA.heatmap.top.%s.%s.%s.jpg",sample,control, date), width = 8.5, height = 11, units="in", res=300);
# hm <- DimHeatmap(object = scrna_GEX, dims = 1, cells = 500, balanced = TRUE);
# print(hm);
# dev.off();

# jpeg(sprintf("PCA.heatmap.multi.%s.%s.%s.jpg",sample,control, date), width = 8.5, height = 24, units="in", res=300);
# hm.multi <- DimHeatmap(object = scrna_GEX, dims = 1:10, cells = 500, balanced = TRUE);
# print(hm.multi);
# dev.off();

# scrna_GEX <- FindNeighbors(scrna_GEX, reduction = "pca", dims = 1:nPC)
# scrna_GEX <- FindClusters(scrna_GEX, resolution = 0.5)
# scrna_GEX[[sprintf("ClusterNames_%.1f_%dPC",0.5, nPC)]] <- Idents(scrna_GEX)
# scrna_GEX <- FindClusters(scrna_GEX, resolution = 0.7)
# scrna_GEX[[sprintf("ClusterNames_%.1f_%dPC",0.7, nPC)]] <- Idents(scrna_GEX)
# scrna_GEX <- FindClusters(scrna_GEX, resolution = 0.9)
# scrna_GEX[[sprintf("ClusterNames_%.1f_%dPC",0.9, nPC)]] <- Idents(scrna_GEX)

# scrna_GEX <- FindClusters(scrna_GEX, resolution = 1.2)
# scrna_GEX[[sprintf("ClusterNames_%.1f_%dPC",1.2, nPC)]] <- Idents(scrna_GEX)
# Idents(scrna_GEX) <- sprintf("ClusterNames_%.1f_%dPC",1.2, nPC)


# n.graph = length(unique(Idents(scrna_GEX)))
# print(n.graph)
# rainbow.colors = rainbow(n.graph, s=0.6, v=0.9);
# names(rainbow.colors) <- sort(unique(Idents(scrna_GEX)))
# print(rainbow.colors)
# print(length(rainbow.colors))

# jpeg(sprintf("UMAP.clusters.%s.%d.%.1f.%s.jpg",sample,nPC, 1.2, date), width = 10, height = 8, units="in", res=300);
# p2 <- DimPlot(object = scrna_GEX, reduction = "umap", group.by = sprintf("ClusterNames_%.1f_%dPC",0.7, nPC), cols = rainbow.colors, pt.size=0.1) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
# print(p2);
# dev.off();

# jpeg(sprintf("UMAP.clusters.labeled.%s.%d.%.1f.%s.jpg",sample,nPC, 1.2, date), width = 10, height = 8, units="in", res=300);
# p2 <- DimPlot(object = scrna_GEX, reduction = "umap", group.by = sprintf("ClusterNames_%.1f_%dPC",0.7, nPC), cols = rainbow.colors, pt.size=0.1, label=TRUE,label.size = 5) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
# print(p2);
# dev.off();

# print ("color UMAP by Principal Components");
# jpeg(sprintf("UMAP.%d.%.1f.colorby.PCs.%s.%s.jpg",sample,nPC, 1.2, date), width = 12, height = 6, units="in", res=100);
# redblue=c("blue","gray","red");
# fp1 <- FeaturePlot(object = scrna_GEX, features = 'PC_1', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
# fp2 <- FeaturePlot(object = scrna_GEX, features = 'PC_2', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
# fp3 <- FeaturePlot(object = scrna_GEX, features = 'PC_3', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
# fp4 <- FeaturePlot(object = scrna_GEX, features = 'PC_4', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
# fp5 <- FeaturePlot(object = scrna_GEX, features = 'PC_5', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
# fp6 <- FeaturePlot(object = scrna_GEX, features = 'PC_6', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
# fp7 <- FeaturePlot(object = scrna_GEX, features = 'PC_7', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
# fp8 <- FeaturePlot(object = scrna_GEX, features = 'PC_8', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
# fp9 <- FeaturePlot(object = scrna_GEX, features = 'PC_9', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
# fp10 <- FeaturePlot(object = scrna_GEX, features = 'PC_10', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
# print(plot_grid(fp1, fp2, fp3, fp4, fp5, fp6, fp7, fp8, fp9, fp10));
# print(plot_grid(fp1, fp2, fp3, fp4, fp5, fp6, fp7, fp8, fp9));
# dev.off();

#only for testing purposes
# set.seed(100)
# select 2500 or all cells from single cell object after filtering (whichever is minimum) of cells from filtered subset randomly for testing the script.
# random_sample_of_cells = sample(Cells(scrna_GEX),min(4000,length(Cells(scrna_GEX))))
# scrna_GEX <- subset(scrna_GEX,cells=random_sample_of_cells)


saveRDS(scrna_GEX, file = sprintf("%s.MergedFilteredSeuratObject.%s.rds", sample_name, date))

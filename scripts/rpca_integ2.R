.libPaths( c("/storage1/fs1/allegra.petti/Active/R_libs_scratch/RLibs_4.0.3",.libPaths()) )
library(SingleR)
library(SingleCellExperiment)
library(scuttle)
library(scran)
library(scater)
library(scRNAseq)
library(Seurat)
library(BiocParallel)
library(yaml)
library(Ckmeans.1d.dp)
library(stringr)
library(dplyr)
library(sctransform)
library(ggplot2)
library(sctransform)
library(future)

library(RColorBrewer)


options(future.globals.maxSize = 16000 * 1024^3)

args = commandArgs(trailingOnly=TRUE);
print (length(args))

if(length(args) != 5)
{
  stop("\n adequate number of arguments not given \n")
}

rds_file <- args[1]
prefix <- args[2]
subset_col <- args[3]
subset_sel <- args[4]
ident_name <- args[5]


seurat_obj <- readRDS(rds_file)

date = gsub("2021-","21",Sys.Date(),perl=TRUE);
date = gsub("-","",date);

# #debugging mode
# print('#####rds_file####')
# print(rds_file)

# print('#####prefix####')
# print(prefix)

# print('#####subset_col####')
# print(subset_col)

# print('#####subset_sel####')
# print(subset_sel)

print('#####BEFORE SUBSETTING SEURAT OBJECT####')
print(seurat_obj)


DefaultAssay(seurat_obj) <- 'RNA'

if(subset_col!='NA'){
Idents(seurat_obj) <- subset_col
seurat_obj <- subset(seurat_obj,idents=subset_sel)
}


print('#####AFTER SUBSETTING SEURAT OBJECT####')
print(seurat_obj)

print('#####Default assay####')
print(DefaultAssay(seurat_obj))

# I had a similar issue when I switched from seurat 3 to 4, for me 
# the issue was the introduction of the k.weight parameter in IntegrateData, 
# which by default is set at 100, but should be less than the smallest number of cells 
# you're integrating with. So, for example, if organoid.list contains objects with 100,
 # 40, 200, and 300 cells, I would run
# organoid.combined <- IntegrateData(anchorset = organoid.anchors, dims = 1:20,
 # k.weight = 39)
# I hope this helps!


rpca_integration <- function(seurat_obj,prefix) {
cluster <- SplitObject(seurat_obj, split.by = "orig.ident")

def_k_wt = min(table(seurat_obj$orig.ident))-1

 for (j in 1:length(cluster)) {
    
    if (length(colnames(cluster[[j]]@assays$RNA@counts)) < 40) {
      # cluster[[j]] <- NULL
      
      print(paste0("Excluding ", prefix, " ", names(cluster[j])))
      exclude <- c(exclude, names(cluster[j]))
      next
      
    } else {
      
      print(paste0("Normalizing ", prefix, " ", names(cluster[j])))
      cluster[[j]] <- NormalizeData(cluster[[j]], verbose = FALSE)
      print(paste0("Finding variable features of ", prefix, " ", names(cluster[j])))
      cluster[[j]] <- FindVariableFeatures(cluster[[j]], selection.method = "vst", nfeatures = 3000, verbose = FALSE)
      
    }
  }

  features <- SelectIntegrationFeatures(object.list = cluster)

for (j in 1:length(cluster)) {
    print(paste0("Scaling ",prefix, " ", names(cluster[j])))
    cluster[[j]] <- ScaleData(cluster[[j]], features = features, verbose = FALSE)
    cluster[[j]] <- RunPCA(cluster[[j]], features = features, verbose = FALSE, npcs = 30)
  }

  anchors <- FindIntegrationAnchors(object.list = cluster, anchor.features = features, reduction = "rpca", verbose = FALSE)  

integ.features <- rownames(seurat_obj@assays$RNA@counts)

  print(paste0("Integrating ", prefix))
  combined <- IntegrateData(anchorset = anchors, k.weight = def_k_wt,features.to.integrate = integ.features,verbose = FALSE)

DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = c(0.1,0.3,0.5, 0.7, 0.9))

# merge.tech <- DimPlot(combined, label = T, group.by = "technique", raster = F)
# merge.samp <- DimPlot(combined, label = T, group.by = "orig.ident", raster = F)

# jpeg(filename = paste0("umap_technique_merged_rpca-integ",prefix,"_cells.jpg"),height = 5,width = 5,res = 300,units = "in")
# print(merge.tech)
# dev.off()

# jpeg(filename = paste0("umap_samples_merged_rpca-integ",prefix,"_cells.jpg"),height = 5, width = 5,res = 300,units = "in")
# print(merge.samp)
# dev.off()

return(combined)
}

integ_obj <- rpca_integration(seurat_obj,prefix)

DefaultAssay(integ_obj) <- "RNA"

Idents(integ_obj) <- ident_name

DEGs <- FindAllMarkers(object=integ_obj,only.pos = T);

write.table(DEGs, file=sprintf("DEGs.Wilcox.%s.%s_onlypos.%s.txt",ident_name, prefix, date), quote=FALSE, sep="\t", row.names=FALSE) # must save cluster-specific marker genes

qual_cols = brewer.pal.info[brewer.pal.info$category == "qual", ]
qual_cols <- qual_cols[qual_cols$colorblind==T,]
qual_cols
col_vector = unlist(mapply(brewer.pal, qual_cols$maxcolors, rownames(qual_cols)))

sel_sig_colors <- sample(col_vector,length(sort(unique(Idents(integ_obj)))))
names(sel_sig_colors) <- sort(unique(Idents(integ_obj)))

jpeg(sprintf("%s.%s.labeled.%s.jpg",ident_name, prefix, date), width = 15, height = 15, units="cm", res=300);
p2 <- DimPlot(object = integ_obj, reduction = "umap", group.by = ident_name, cols = sel_sig_colors, pt.size=0.1,label=T,label.size = 6)+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
print(p2);
dev.off()

jpeg(sprintf("%s.%s.%s.jpg",ident_name,prefix, date), width = 15, height = 15, units="cm", res=300);
p2 <- DimPlot(object = integ_obj, reduction = "umap", group.by = ident_name, cols = sel_sig_colors, pt.size=0.1)+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
print(p2);
dev.off()


Idents(integ_obj) <- 'orig.ident'

sel_sig_colors <- sample(col_vector,length(sort(unique(Idents(integ_obj)))))
names(sel_sig_colors) <- sort(unique(Idents(integ_obj)))

jpeg(sprintf("%s.orig.ident.%s.labeled.jpg", prefix, date), width = 15, height = 15, units="cm", res=300);
p2 <- DimPlot(object = integ_obj, reduction = "umap", group.by = 'orig.ident', cols = sel_sig_colors, pt.size=0.1,label=T,label.size = 6)+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
print(p2);
dev.off()

jpeg(sprintf("%s.orig.ident.%s.jpg", prefix, date), width = 15, height = 15, units="cm", res=300);
p2 <- DimPlot(object = integ_obj, reduction = "umap", group.by = 'orig.ident', cols = sel_sig_colors, pt.size=0.1)+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
print(p2);
dev.off()

saveRDS(integ_obj,file=sprintf("%s_rpca_%s_ident.%s.rds", prefix,ident_name,date))

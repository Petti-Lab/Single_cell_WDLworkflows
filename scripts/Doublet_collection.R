.libPaths( c("/storage1/fs1/allegra.petti/Active/R_libs_scratch/RLibs_4.0.3",.libPaths()) )
library(DoubletCollection)
library(DoubletFinder)
library(scds)
library(scDblFinder)
library(scran)
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(Matrix)
library(data.table)

# koetjen/rstudio:seurat
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2) {
  args <- c("--help")
}

Seurat_10x_directory <- as.character(args[1])
# raw_barcodes <- as.character(args[2])
sample_name <- as.character(args[2])
# /storage1/fs1/allegra.petti/Active/spatial_snRNAseq_gbm/scrnaseq/SAMPLES/B183_1/outs/raw_feature_bc_matrix/barcodes.tsv.gz
# Seurat_10x_directory='/storage1/fs1/allegra.petti/Active/spatial_snRNAseq_gbm/scrnaseq/SAMPLES/B183_1/outs/filtered_feature_bc_matrix'
# dt = fread("/storage1/fs1/allegra.petti/Active/spatial_snRNAseq_gbm/scrnaseq/SAMPLES/B183_1/outs/raw_feature_bc_matrix/barcodes.tsv.gz")
print('input Seurat_10x_directory')
print(Seurat_10x_directory)
# dt = fread(raw_barcodes,header=F)
# raw_cell_count <- length(dt$V1)

# B178-01 (fresh)
# B183-01 (fresh)
# B185-03 (unsorted)
# B186
# B189
# WU-1200 Core CD45-
# WU1200_Core_CD45neg
#   11:32
# Also, Neftel 10x and Neftel Smartseq
# Error: package or namespace load failed for ‘DoubletCollection’ in loadNamespace(i, c(lib.loc, .libPaths()), versionChec
# k = vI[[i]]):
#/storage1/fs1/allegra.petti/Active/spatial_snRNAseq_gbm/scrnaseq/SAMPLES/B186/outs/filtered_feature_bc_matrix
# /raw_feature_bc_matrix/barcodes.tsv.gz
#data.10x <- Read10X(data.dir = Seurat_10x_directory);

data.10x <- gcs_get_object(object_name = Seurat_10x_directory, parseFunction = Read10X())


methods <- c('doubletCells','cxds','bcds','hybrid','scDblFinder','DoubletFinder')
score.list <- FindScores(data.10x, methods)

#Not meant for this kind of doublet
# get_doublet_rate_10x <- function(no_of_cells,raw_cell_count=750000){
# fraction_no_cells <- exp(-no_of_cells/raw_cell_count)/1
# fraction_1_cell <- (no_of_cells/raw_cell_count)*exp(-no_of_cells/raw_cell_count)/1
# fraction_doublet <- 1 - fraction_no_cells - fraction_1_cell
# return(fraction_doublet)
# }

#10x doublet rates 

doublet_results <- FindDoublets(score.list=score.list,rate= 0.08)

myres_list <- list(doublet_scores=score.list,doublet_results=doublet_results,cellnames=colnames(data.10x))

saveRDS(object = myres_list,file = sprintf("%s_Doublet_collection_results.rds",sample_name))

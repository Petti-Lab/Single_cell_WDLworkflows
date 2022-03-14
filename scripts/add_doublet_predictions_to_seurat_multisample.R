.libPaths( c("/storage1/fs1/allegra.petti/Active/R_libs_scratch/RLibs_4.0.3",.libPaths()) )
library(DoubletFinder)
library(scds)
library(scDblFinder)
library(scran)
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(Matrix)
library(DoubletCollection)
library(Seurat)
library(BiocParallel)
library(yaml)
library(Ckmeans.1d.dp)
library(stringr)

Seurat_file <- as.character(args[1])
input_tsv_file <- as.character(args[2])
doublet_files <- as.character(args[3])


add_doublet_predictions_to_seurat_singlesample <- function(seurat_object,doublet_object){
  doublet_object[['doublet_results']] <- FindDoublets(score.list=doublet_object$doublet_scores,rate=0.08)
  for(i in names(doublet_object[['doublet_results']])){
    sel_index <- doublet_object[['doublet_results']][[i]]
    doublet_cells <- doublet_object$cellnames[sel_index]
    non_doublet_cells <- setdiff(doublet_object$cellnames,doublet_cells)
    doublet_cells_pred <- rep('yes',length(doublet_cells))
    names(doublet_cells_pred) <- doublet_cells
    non_doublet_cells_pred <- rep('no',length(non_doublet_cells))
    names(non_doublet_cells_pred) <- non_doublet_cells
    doublet_preds <- c(doublet_cells_pred,non_doublet_cells_pred)
    sel_cells <- intersect(Cells(seurat_object),names(doublet_preds))
    sel_doublet_preds <- doublet_preds[sel_cells]
    index <- match(Cells(seurat_object),names(sel_doublet_preds))
    seurat_object[[sprintf("doublet_results_%s",i)]] <- sel_doublet_preds[index]
  }
  for(i in names(doublet_object[['doublet_scores']])){
    doubscores <- doublet_object[['doublet_scores']][[i]]
    names(doubscores) <- doublet_object$cellnames
    sel_cells <- intersect(Cells(seurat_object),names(doubscores))
    sel_doubscores <- doubscores[sel_cells]
    index <- match(Cells(seurat_object),names(sel_doubscores))
    seurat_object[[sprintf("doublet_scores_%s",i)]] <- sel_doubscores[index]
  }
  return(seurat_object)
}

# df$count <- rowSums(df[c(1,3)] == "Yes")

input_df <- read.table(input_tsv_file,sep="\t",header=FALSE)
colnames(input_df) <- c('Sample','cellranger_10x_directory')

doublet_sample_list <- list()


for(i in 1:nrow(input_df)){
doublet_sample_list[[input_df$Sample[i]]] <- doublet_files[grep(input_df$Sample[i],doublet_files)]
}

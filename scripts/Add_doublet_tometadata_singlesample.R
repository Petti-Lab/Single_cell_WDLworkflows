# .libPaths( c("/storage1/fs1/allegra.petti/Active/R_libs_scratch/RLibs_4.0.3",.libPaths()) )
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
library(RColorBrewer)
library(gridExtra)


args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2) {
  args <- c("--help")
}

Seurat_file <- as.character(args[1])
doublet_file <- as.character(args[2])
sample_name <- as.character(args[3])


print(Seurat_file)
print(doublet_file)

date = gsub("2023-","23",Sys.Date(),perl=TRUE);
date = gsub("-","",date);

seurat_object <- readRDS(Seurat_file)
doublet_object <- readRDS(doublet_file)

add_doublet_predictions_to_seurat_singlesample <- function(seurat_object,doublet_object){
  # doublet_object[['doublet_results']] <- FindDoublets(score.list=doublet_object$doublet_scores,rate=0.08)
  for(i in names(doublet_object$doublet_results)){
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
  for(i in names(doublet_object$doublet_scores)){
    doubscores <- doublet_object[['doublet_scores']][[i]]
    names(doubscores) <- doublet_object$cellnames
    sel_cells <- intersect(Cells(seurat_object),names(doubscores))
    sel_doubscores <- doubscores[sel_cells]
    index <- match(Cells(seurat_object),names(sel_doubscores))
    seurat_object[[sprintf("doublet_scores_%s",i)]] <- sel_doubscores[index]
  }
  return(seurat_object)
}

seurat_object <- add_doublet_predictions_to_seurat_singlesample(seurat_object=seurat_object,doublet_object=doublet_object)


All_doublet_preds <- seurat_object@meta.data[grep('doublet_results_',names(seurat_object@meta.data))]

All_doublet_scores <- seurat_object@meta.data[grep('doublet_scores_',names(seurat_object@meta.data))]
# All_doublet_preds <- do.call(rbind,unname(doublet_prediction_df))

ncol_doubs <- ncol(All_doublet_preds)
maj <- round((ncol_doubs/2)+1, digits = 0)

All_doublet_preds$prediction_aggrement <- rowSums(All_doublet_preds=='yes')
#doublet predictor aggrement for each cell

All_doublet_preds$majority_doublet_predictions <- 'no'
All_doublet_preds$majority_doublet_predictions[All_doublet_preds$prediction_aggrement >= maj]  <- 'yes'
#majority doublets : predictions for which majority of doublet predictors agree

All_doublet_preds$prediction_aggrement <- NULL

index <- match(Cells(seurat_object),rownames(All_doublet_preds))

seurat_object$majority_doublet_predictions <- All_doublet_preds$majority_doublet_predictions[index]

# majority_doublet_predictions

#Any doublet prediction

sel_cells_index = unique(unname(unlist(doublet_object$doublet_results)))
any_doublet_pred_meth = doublet_object$cellnames[sel_cells_index]
non_doublet_cells <- setdiff(doublet_object$cellnames,any_doublet_pred_meth)
doublet_cells_pred <- rep('yes',length(any_doublet_pred_meth))
names(doublet_cells_pred) <- any_doublet_pred_meth
non_doublet_cells_pred <- rep('no',length(non_doublet_cells))
names(non_doublet_cells_pred) <- non_doublet_cells
doublet_preds <- c(doublet_cells_pred,non_doublet_cells_pred)
sel_cells <- intersect(Cells(seurat_object),names(doublet_preds))
sel_doublet_preds <- doublet_preds[sel_cells]
index <- match(Cells(seurat_object),names(sel_doublet_preds))
seurat_object[[sprintf("doublet_results_%s",'any_method')]] <- sel_doublet_preds[index]

All_doublet_preds <- seurat_object@meta.data[grep('doublet_results_',names(seurat_object@meta.data))]

feature.pal = rev(colorRampPalette(brewer.pal(11,"Spectral"))(20))

gplotlist = list()
for(i in colnames(All_doublet_scores)){
  dm <-FeaturePlot(seurat_object,features =c(i),cols = feature.pal,label = T,pt.size = 0.8,label.size = 15)
  dm <- dm + theme(text = element_text(size = 22)) + theme(legend.title=element_text(color="black",size=20))+ theme(legend.text=element_text(size=20))+guides(fill = guide_legend(override.aes = list(size=15)),colour = guide_colourbar(barwidth =10,barheight=20))
  dm <- dm +theme(axis.text.y = element_text(color="black",size=22))+theme(axis.text.x = element_text(color="black",size=22))+theme(axis.title.x = element_text(color="black",size=22))+theme(axis.title.y = element_text(color="black",size=22))
  dm <- dm+ggtitle(i) + theme(plot.title = element_text(hjust = 0.5))
  gplotlist[[i]] = dm
  # ggsave(sprintf("featureplot_%s_%s_lab_%s.png",i,sample_name,date),plot = dm, width = 30, height = 30, units = "in",dpi = 300,device = "png",scale = 1)
}

ml <- marrangeGrob(gplotlist, nrow=2, ncol=2,top=quote(sample_name))
ggsave(sprintf("featureplot_%s_%s_lab_%s.pdf","doublet_scores",sample_name,date),ml,width=20, height=20)
 # /storage1/fs1/allegra.petti/Active/Users/khan.saad/WDL_pipelines/scatter_doublet/b4473306-6ff6-4a81-86b0-b44e275c90b7/call-add_doublets_metadata_tomultisample_seurat/execution/

doublet_cols <- c("darkred","grey88")
names(doublet_cols) <- c("yes","no")

gplotlist = list()
for(i in c(colnames(All_doublet_preds),'majority_doublet_predictions')){
  Idents(seurat_object) <- i
  mydmplt <-DimPlot(seurat_object,cols= doublet_cols,pt.size = 0.8) + ggtitle(i) + theme(plot.title = element_text(hjust = 0.5))
  mydmplt <- mydmplt + theme(text = element_text(size = 22)) + theme(legend.title=element_text(color="black",size=20))+ theme(legend.text=element_text(size=20))+guides(fill = guide_legend(override.aes = list(size=15)))
  mydmplt <- mydmplt +theme(axis.text.y = element_text(color="black",size=22))+theme(axis.text.x = element_text(color="black",size=22))+theme(axis.title.x = element_text(color="black",size=22))+theme(axis.title.y = element_text(color="black",size=22))
  gplotlist[[i]] = mydmplt
  }
ml <- marrangeGrob(gplotlist, nrow=2, ncol=2,top=quote(sample_name))
ggsave(sprintf("Dimplot_%s_%s_%s.pdf","doublet_results",sample_name,date),ml,width=20, height=20)
# ggsave(sprintf("Dimplot_%s_%s_%s.png",i,sample_name,date),plot = mydmplt, width = 30, height = 30, units = "in",dpi = 300,device = "png",scale = 1)

saveRDS(seurat_object, file = paste0('Cycling.SCT.PCA.UMAP.TSNE.CLUST.',sample_name,".doublet_calls.",date,".rds"))
#saveRDS(seurat_object, file = outfile_seurat)

# seurat_object <- add_doublet_predictions_to_seurat_singlesample(seurat_object=seurat_object,doublet_object=doublet_object)
# saveRDS(object = seurat_object,file=outfile_seurat)

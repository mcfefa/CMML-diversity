# Import mitoclone object, plot in palantir embeddings

library(mitoClone)
library(deepSNV)
library(GenomicRanges)
library(pheatmap)
library(Seurat)
library(ggplot2)
library(gridExtra)
library(dplyr)
library("reticulate")
setwd('~')
Sys.setenv(RETICULATE_PYTHON = "./anaconda3/envs/mitoclone/bin/python", force = T)
knitr::opts_chunk$set(echo = TRUE)
py_config()

rm(list = ls())

# Set dir to save to and from which we pull countTables
dir <- '/Users/Brian/scCode/mitoclone/'
saveDir <- paste0(dir, 'mito_in_palantir_params_from_02022022_run_05082022/')
mito_object_dir <- paste0(dir, "mitoclone_objects_post_gurobi_02022022/")

# Import seurat stuff and organize it (now including cell type from singleR)
#source(paste0(dir, 'seurat_import_cleaning_02022022.R'))
#saveRDS(umap_embeddings, '~/scCode/mitoclone/umap_embeddings_pre_mito_02022022.rds')
umap_embeddings <- readRDS('~/scCode/mitoclone/umap_embeddings_pre_mito_02022022.rds')

#Run each sample separately (no blacklist)
ptm <- proc.time()

mito_objects <- list.files(mito_object_dir, full.names = T)
plot_list <- list()

for (i in c(1:39)){
  
  # Get current patient "name"
  currentSample <- paste0("CMML_", i)
  
  # 1, 8, and 29 had no sites (only one clone, plot accordingly)
  if (i %in% c(1, 8, 29)) {
    
    # Get only embeddings for current patient
    embeddings_cut <- umap_embeddings[grep(paste0(currentSample, "_"), row.names(umap_embeddings)),]
    
    p1 <- ggplot(embeddings_cut) + geom_jitter(aes(x = palantir_embedding_X, y = palantir_embedding_Y, color = "#F8766D"), size = .6) +
      theme_bw()+ theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(face = "bold", size = 15),
                        axis.ticks = element_blank(), axis.text = element_blank(),
                        axis.title = element_blank())+ggtitle(paste0(currentSample))+
      guides(color="none")
  } else {
    if (paste0(mito_object_dir, '/cmml_', i, '_clustered_mitoclone_object_02022022.rds') %in% mito_objects) {
      CMML <- readRDS(paste0(mito_object_dir, 'cmml_', i, '_clustered_mitoclone_object_02022022.rds'))
    } else {
      next
    }
    
    #Get clonal info from mutaCluster output
    cmml_clonalInfo <- data.frame(CMML@cell2clone)
    clone <- data.frame("CB" = row.names(cmml_clonalInfo), "Clone" = NA)
    for (j in 1:length(cmml_clonalInfo[,1])){
      clone[j, "Clone"] <- which.max(cmml_clonalInfo[j,])[[1]]
    }
    
    barcodesFromMito <- clone[,1]
    barcodesOnlyFromMito <- sub("CB_", "", barcodesFromMito)
    barcodesOnlyMito <- gsub("-1.*","",barcodesOnlyFromMito)
    clone$CB <- barcodesOnlyMito
    clone$Sample <- substr(clone$CB, 1, 7)
    
    
    #Match seurat and mitoclone barcodes
    combine.df <- data.frame("CB" = NA, "X" = NA, "Y" = NA, "Clone" = NA, "Clus" = NA, "Sample" = NA)
    barcodes <- row.names(umap_embeddings)
    for (j in 1:length(clone$CB)){
      if (clone$CB[j] %in% barcodes){
        combine.df[j,"Clone"] <- clone$Clone[j]
        combine.df[j,"CB"] <- clone$CB[j]
        index <- which(row.names(umap_embeddings) %in% clone$CB[j])
        combine.df[j,"X"] <- umap_embeddings[index, "X"]
        combine.df[j,"Y"] <- umap_embeddings[index,"Y"]
        combine.df[j,"Palantir_X"] <- umap_embeddings[index, "palantir_embedding_X"]
        combine.df[j,"Palantir_Y"] <- umap_embeddings[index,"palantir_embedding_Y"]
        combine.df[j, "Clus"] <- as.numeric(as.character(umap_embeddings[index, "Clus"]))
        combine.df[j, "cellType"] <- umap_embeddings[index, "cellType"]
        combine.df[j, "Sample"] <- clone$Sample[j]
      }
    }
    
    # Add binary mutational status to combine.df
    binary_mut_status <- mutate_all(as.data.frame(CMML@ternary)[1:dim(combine.df)[1],], function(x) as.integer(x))
    row.names(binary_mut_status) <- barcodesOnlyMito[1:dim(combine.df)[1]]
    combine_muts <- cbind(combine.df, binary_mut_status)
    
    # Remove na's (cells that were filtered by our QC)
    combine_muts <- combine_muts[!is.na(combine_muts$CB),]
    
    #pdf(paste0(saveDir, 'cmml', i, '_palantir_embeddings_mitoclone.pdf'), width = 6, height = 4.5)
    p1 <- ggplot(combine_muts) + geom_jitter(aes(x = Palantir_X, y = Palantir_Y, color = as.factor(Clone)), size = .6) +
      theme_bw()+ theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(face = "bold", size = 15),
                        axis.ticks = element_blank(), axis.text = element_blank(),
                        axis.title = element_blank())+ggtitle(paste0(currentSample))+
      guides(color="none")
  }

  plot_list[[i]] <- p1
  rm(CMML)
}
proc.time() - ptm


plot_list <- plot_list[lengths(plot_list) != 0]

pdf("~/scCode/plot_all_mito_in_palantir_v2.pdf", height = 16, width = 11)
do.call("grid.arrange", c(plot_list, ncol=5))
dev.off()

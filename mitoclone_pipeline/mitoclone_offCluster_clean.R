library(mitoClone)
library(deepSNV)
library(GenomicRanges)
library(devtools)
library(pheatmap)
library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(RColorBrewer)
library(gplots)


rm(list=ls())

#Change plotTree function so that it actually works (period gave errors, gsub to _)
plotTree <- function(mutcalls, file = "mytree.ps") {
  gv <- toGraphviz(mutcalls@tree)
  gv <- gsub("\\.", "_", gv)
  tmp <- tempfile()
  writeLines(gv, con = tmp)
  system(sprintf("dot -Tps %s > %s", tmp, file),wait = T)
  file.remove(tmp)
}

#Clean mitoclone analysis (once off the cluster)
dir <- '/Users/4472241/scCode/mitoclone/mitoclone_RDS_defaultParams_MINREADS3_MINCELL10_MINFRAC0.05_MINRelative.Other0.99/'
outdir <- '/Users/4472241/scCode/mitoclonePNG_04-23-2021/'
n_samples <- 20
doneRunning <- c(1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21:39)
for(i in 20:20){
  
  #Read in the data from the cluster
  P_sample <- readRDS(paste0(dir, 'P', i, '_results_mitoclone.rds'))
  
  #Cluster metaclones
  #CMML_sample_tree <- quick_cluster(P_sample)
  CMML_sample_tree <- clusterMetaclones(P_sample, min.lik = 1)
  
  #Plot the phylogenetic tree and save
  plotTree(CMML_sample_tree, file = paste0(outdir, "P", i, "_phylogenetic.ps"))
  
  #Plot the clones
  png(paste0(outdir, "p", i, "plotClones.png"))    
  print(plotClones(CMML_sample_tree))
  dev.off()
  
  #Configure into palantir gene expression FDL
  paldir <- '/Users/4472241/scCode/runPalantir/outputNew/'
  palantir_embeddings <- read.csv(paste0(paldir, 'cmml', i,'fdl.csv'))
  #cmml_cellNames <- read.csv(paste0(dir, 'cmml', i, 'CellNames.csv'))
  barcodes <- palantir_embeddings$X
  barcodes <- sub(".*_", "", barcodes)
  palantir_embeddings$X <- barcodes
  colnames(palantir_embeddings) <- c("CB", "X", "Y")
  
  #Get clonal info from mutaCluster output
  cmml_clonalInfo <- data.frame(CMML_sample_tree@cell2clone)
  clone <- data.frame("CB"=row.names(cmml_clonalInfo), "Clone" = NA)
  for (j in 1:length(cmml_clonalInfo[,1])){
    clone[j, "Clone"] <- which.max(cmml_clonalInfo[j,])[[1]]
  }
  
  barcodesFromMito <- clone[,1]
  barcodesOnlyFromMito <- sub(".*CB_", "", barcodesFromMito)
  barcodesOnlyMito <- gsub("-1.*","",barcodesOnlyFromMito)
  clone$CB <- barcodesOnlyMito
  
  #Match palantir and mitoclone barcodes
  combine.df <- data.frame("CB" = NA, "X" = NA, "Y" = NA, "Clone" = NA)
  for (j in 1:length(clone$CB)){
    if (clone$CB[j] %in% barcodes){
      combine.df[j,"Clone"] <- clone$Clone[j]
      combine.df[j,"CB"] <- clone$CB[j]
      index <- which(palantir_embeddings$CB %in% clone$CB[j])
      combine.df[j,"X"] <- palantir_embeddings[index, "X"]
      combine.df[j,"Y"] <- palantir_embeddings[index,"Y"]
    }
  }
  
  #Remove nas, convert clone to factor in combine.df and plot
  combine.df <- combine.df[!is.na(combine.df$CB),]
  combine.df$Clone <- as.factor(combine.df$Clone)
  
  png(paste0(outdir, 'palantirFDL_mitocloneColored_CMML_', i, '.png'))
  print(ggplot(combine.df, aes(x=X, y=Y, color = Clone))+ geom_point(size = .5) +
          theme(panel.background = element_rect(fill = "white"))+
          scale_color_manual(values=c("blue", "green", "red", "pink", 'orange', 'black', 'yellow',
                                      'purple', 'gray', 'cyan', 'turquoise', "grey1", "grey2",
                                      "grey3", "grey4", "grey5", "grey6", "grey7",
                                      "grey8", "grey8", "grey9", "grey10")))
  dev.off()
  
  ## Now get the palantir embeddings from the normal reference
  palNormdir <- '/Users/4472241/scCode/palantirEmbeddings/mapToRep1/'
  cmml_embeddings <- read.csv(paste0(palNormdir, 'cmml', i,'PalantirEmbeddings.csv'))
  cmml_cellNames <- read.csv(paste0(palNormdir, 'cmml', i, 'CellNames.csv'))
  barcodes <- cmml_cellNames$x
  barcodes <- sub(".*_", "", barcodes)
  cmml_embeddings$X <- barcodes
  colnames(cmml_embeddings) <- c("CB", "X", "Y")
  
  #Get clonal info from mutaCluster output
  cmml_clonalInfo <- data.frame(CMML_sample_tree@cell2clone)
  clone <- data.frame("CB"=row.names(cmml_clonalInfo), "Clone" = NA)
  for (j in 1:length(cmml_clonalInfo[,1])){
    clone[j, "Clone"] <- which.max(cmml_clonalInfo[j,])[[1]]
  }
  
  barcodesFromMito <- clone[,1]
  barcodesOnlyFromMito <- sub(".*CB_", "", barcodesFromMito)
  barcodesOnlyMito <- gsub("-1.*","",barcodesOnlyFromMito)
  clone$CB <- barcodesOnlyMito
  
  #Match palantir and mitoclone barcodes
  combine.df <- data.frame("CB" = NA, "X" = NA, "Y" = NA, "Clone" = NA)
  for (j in 1:length(clone$CB)){
    if (clone$CB[j] %in% barcodes){
      combine.df[j,"Clone"] <- clone$Clone[j]
      combine.df[j,"CB"] <- clone$CB[j]
      index <- which(cmml_embeddings$CB %in% clone$CB[j])
      combine.df[j,"X"] <- cmml_embeddings[index, "X"]
      combine.df[j,"Y"] <- cmml_embeddings[index,"Y"]
    }
  }
  
  #Remove nas, convert clone to factor in combine.df and plot
  combine.df <- combine.df[!is.na(combine.df$CB),]
  combine.df$Clone <- as.factor(combine.df$Clone)
  
  png(paste0(outdir, 'palantirNormalEmbedding_mitocloneColored_CMML_', i, '.png'))
  print(ggplot(combine.df, aes(x=X, y=Y, color = Clone))+ geom_point(size = .5) +
          theme(panel.background = element_rect(fill = "white"))+
          scale_color_manual(values=c("blue", "green", "red", "pink", 'orange', 'black', 'yellow',
                                      'purple', 'gray', 'cyan', 'turquoise', "grey1", "grey2",
                                      "grey3", "grey4", "grey5", "grey6", "grey7",
                                      "grey8", "grey8", "grey9", "grey10")))
  dev.off()
}


dir <- '/Users/4472241/scCode/mitoclone/mitoclone_RDS_defaultParams_MINREADS3_MINCELL10_MINFRAC0.05_MINRelative.Other0.99/'
outdir <- '/Users/4472241/scCode/mitoclonePNG_04-23-2021/'
n_samples <- 20
doneRunning <- c(1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,
                 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39)

allSeurat <- readRDS('/Users/4472241/scCode/setNewClustering+keepOldUMAP/NewClusteringOldUMAP_04-30-2021_all_46_samples_umap_layered_noHarmony.rds')

#Configure into UMAP space
umap_embeddings <- data.frame(allSeurat@reductions[['umap']]@cell.embeddings)
#cmml_cellNames <- read.csv(paste0(dir, 'cmml', i, 'CellNames.csv'))
barcodes <- row.names(umap_embeddings)
barcodes <- sub(".*_", "", barcodes)
umap_embeddings$CB <- barcodes
colnames(umap_embeddings) <- c( "X", "Y", "CB")

for(i in 20:20){
  i <- 1
  #Read in the data from the cluster
  P_sample <- readRDS(paste0(dir, 'P', i, '_results_mitoclone.rds'))
  
  #Cluster metaclones
  #CMML_sample_tree <- quick_cluster(P_sample)
  CMML_sample_tree <- clusterMetaclones(P_sample, min.lik = 1)
  
  #Plot the phylogenetic tree and save
  #plotTree(CMML_sample_tree, file = paste0(outdir, "P", i, "_phylogenetic.ps"))
  
  #Plot the clones
  #png(paste0(outdir, "p", i, "plotClones.png"))    
  #print(plotClones(CMML_sample_tree))
  #dev.off()
  
  #Get clonal info from mutaCluster output
  cmml_clonalInfo <- data.frame(CMML_sample_tree@cell2clone)
  clone <- data.frame("CB"=row.names(cmml_clonalInfo), "Clone" = NA)
  for (j in 1:length(cmml_clonalInfo[,1])){
    clone[j, "Clone"] <- which.max(cmml_clonalInfo[j,])[[1]]
  }
  
  barcodesFromMito <- clone[,1]
  barcodesOnlyFromMito <- sub(".*CB_", "", barcodesFromMito)
  barcodesOnlyMito <- gsub("-1.*","",barcodesOnlyFromMito)
  clone$CB <- barcodesOnlyMito
  
  #Match palantir and mitoclone barcodes
  combine.df <- data.frame("CB" = NA, "X" = NA, "Y" = NA, "Clone" = NA)
  for (j in 14:length(clone$CB)){
    if (clone$CB[j] %in% barcodes){
      combine.df[j,"Clone"] <- clone$Clone[j]
      combine.df[j,"CB"] <- clone$CB[j]
      index <- which(umap_embeddings$CB %in% clone$CB[j])
      combine.df[j,"X"] <- umap_embeddings[index, "X"]
      combine.df[j,"Y"] <- umap_embeddings[index,"Y"]
    }
  }
  
  #Remove nas, convert clone to factor in combine.df and plot
  combine.df <- combine.df[!is.na(combine.df$CB),]
  combine.df$Clone <- as.factor(combine.df$Clone)
  
  png(paste0(outdir, 'palantirFDL_mitocloneColored_CMML_', i, '.png'))
  print(ggplot(combine.df, aes(x=X, y=Y, color = Clone))+ geom_point(size = .5) +
          theme(panel.background = element_rect(fill = "white"))+
          scale_color_manual(values=c("blue", "green", "red", "pink", 'orange', 'black', 'yellow',
                                      'purple', 'gray', 'cyan', 'turquoise', "grey1", "grey2",
                                      "grey3", "grey4", "grey5", "grey6", "grey7",
                                      "grey8", "grey8", "grey9", "grey10")))
  dev.off()
  
  ## Now get the palantir embeddings from the normal reference
  palNormdir <- '/Users/4472241/scCode/palantirEmbeddings/mapToRep1/'
  cmml_embeddings <- read.csv(paste0(palNormdir, 'cmml', i,'PalantirEmbeddings.csv'))
  cmml_cellNames <- read.csv(paste0(palNormdir, 'cmml', i, 'CellNames.csv'))
  barcodes <- cmml_cellNames$x
  barcodes <- sub(".*_", "", barcodes)
  cmml_embeddings$X <- barcodes
  colnames(cmml_embeddings) <- c("CB", "X", "Y")
  
  #Get clonal info from mutaCluster output
  cmml_clonalInfo <- data.frame(CMML_sample_tree@cell2clone)
  clone <- data.frame("CB"=row.names(cmml_clonalInfo), "Clone" = NA)
  for (j in 1:length(cmml_clonalInfo[,1])){
    clone[j, "Clone"] <- which.max(cmml_clonalInfo[j,])[[1]]
  }
  
  barcodesFromMito <- clone[,1]
  barcodesOnlyFromMito <- sub(".*CB_", "", barcodesFromMito)
  barcodesOnlyMito <- gsub("-1.*","",barcodesOnlyFromMito)
  clone$CB <- barcodesOnlyMito
  
  #Match palantir and mitoclone barcodes
  combine.df <- data.frame("CB" = NA, "X" = NA, "Y" = NA, "Clone" = NA)
  for (j in 1:length(clone$CB)){
    if (clone$CB[j] %in% barcodes){
      combine.df[j,"Clone"] <- clone$Clone[j]
      combine.df[j,"CB"] <- clone$CB[j]
      index <- which(cmml_embeddings$CB %in% clone$CB[j])
      combine.df[j,"X"] <- cmml_embeddings[index, "X"]
      combine.df[j,"Y"] <- cmml_embeddings[index,"Y"]
    }
  }
  
  #Remove nas, convert clone to factor in combine.df and plot
  combine.df <- combine.df[!is.na(combine.df$CB),]
  combine.df$Clone <- as.factor(combine.df$Clone)
  
  png(paste0(outdir, 'palantirNormalEmbedding_mitocloneColored_CMML_', i, '.png'))
  print(ggplot(combine.df, aes(x=X, y=Y, color = Clone))+ geom_point(size = .5) +
          theme(panel.background = element_rect(fill = "white"))+
          scale_color_manual(values=c("blue", "green", "red", "pink", 'orange', 'black', 'yellow',
                                      'purple', 'gray', 'cyan', 'turquoise', "grey1", "grey2",
                                      "grey3", "grey4", "grey5", "grey6", "grey7",
                                      "grey8", "grey8", "grey9", "grey10")))
  dev.off()
}

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
#directory for mitoclone results:
dir <- '/Users/4472241/scCode/mitoclone/mitoclone_RDS_defaultParams_MINREADS3_MINCELL10_MINFRAC0.05_MINRelative.Other0.99/'
#directory to save to:
outdir <- '/Users/4472241/scCode/mitoclonePNG_04-23-2021/'
#directory for csv files containing FDL embeddings from palantir:
paldir <- '/Users/4472241/scCode/runPalantir/outputNew/'
#directory for csv files containing 
palNormdir <- '/Users/4472241/scCode/palantirEmbeddings/mapToRep1/'

## NOTE: Only samples with more than one clone will plot, others will throw error...just skip those samples
n_samples <- 39
for(i in 1:n_samples){
  
  #Read in the data from the mitoclone analysis
  P_sample <- readRDS(paste0(dir, 'P', i, '_results_mitoclone.rds'))
  
  #Cluster metaclones
  #CMML_sample_tree <- quick_cluster(P_sample)
  CMML_sample_tree <- clusterMetaclones(P_sample, min.lik = 1)
  
  #Plot the phylogenetic tree and save
  plotTree(CMML_sample_tree, file = paste0(outdir, "P", i, "_phylogenetic.ps"))
  dev.off()
  
  #Plot the clones
  png(paste0(outdir, "p", i, "plotClones.png"))    
  print(plotClones(CMML_sample_tree))
  dev.off()
  
  #Configure into palantir gene expression FDL
  palantir_embeddings <- read.csv(paste0(paldir, 'cmml', i,'fdl.csv'))
  #cmml_cellNames <- read.csv(paste0(dir, 'cmml', i, 'CellNames.csv'))
  #Get barcodes from embeddings and keep only the actual barcode
  barcodes <- palantir_embeddings$X
  barcodes <- sub(".*_", "", barcodes)
  palantir_embeddings$X <- barcodes
  colnames(palantir_embeddings) <- c("CB", "X", "Y")
  
  #Get clonal info and barcodes from mutaCluster output
  cmml_clonalInfo <- data.frame(CMML_sample_tree@mainClone)
  clone <- data.frame("CB"=row.names(cmml_clonalInfo), "Clone" = NA)
  for (j in 1:length(cmml_clonalInfo[,1])){
    clone[j, "Clone"] <- which.max(cmml_clonalInfo[j,])[[1]]
  }
  
  #Get barcodes from mitoclone and keep only the actual barcode
  barcodesFromMito <- clone[,1]
  barcodesOnlyFromMito <- sub(".*CB_", "", barcodesFromMito)
  barcodesOnlyMito <- gsub("-1.*","",barcodesOnlyFromMito)
  clone$CB <- barcodesOnlyMito
  
  #Match palantir and mitoclone barcodes, assign matching clonal info and embedding info to combine.df
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
  cmml_embeddings <- read.csv(paste0(palNormdir, 'cmml', i,'PalantirEmbeddings.csv'))
  cmml_cellNames <- read.csv(paste0(palNormdir, 'cmml', i, 'CellNames.csv'))
  barcodes <- cmml_cellNames$x
  barcodes <- sub(".*_", "", barcodes)
  cmml_embeddings$X <- barcodes
  colnames(cmml_embeddings) <- c("CB", "X", "Y")
  
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

library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
library(patchwork)
library(BiocNeighbors)
library(tidyverse)
library(RColorBrewer)
library(gplots)
library(Rmagic)
library(gdata)
library(ggpubr)

rm(list = ls())

#Import the seurat object with all of our samples assigned to clusters
dir <- '/Users/4472241/scCode/redoEntirePipeline_05-05-2021/standard_pipeline_withHarmony/'
allSeurat <- readRDS(paste0(dir, 'seuratObj_05-11-2021_postStandardPipeline_withHarmony_allSamples_39+8.rds'))

#Import the Wu signature
sigDir <- '/Users/4472241/scCode/reRunScores_seurat/'
wuSig <- read.csv(paste0(sigDir, 'Wu_Signatures.csv'), skip = 2)

#Convert each specific column to element of list
wuSig.list <- list()
for (i in 1:length(wuSig[1,])){
  wuSig.list[[i]] <- wuSig[2:101,i]
}
names(wuSig.list) <- colnames(wuSig)

allSeurat <- AddModuleScore(allSeurat, wuSig.list, search = T)
allSeurat$wu_HSC <- allSeurat$Cluster1
allSeurat$wu_MEP <- allSeurat$Cluster2
allSeurat$wu_GMP <- allSeurat$Cluster3
allSeurat$wu_ProB <- allSeurat$Cluster4
allSeurat$wu_ETP <- allSeurat$Cluster5

FeaturePlot(allSeurat, features = "wu_GMP")


#Import the VanGalen signature
VGSig <- read.xls(paste0(sigDir, 'VanGalen_mmc3_GeneListClassifier.xlsx'), skip = 1)
#Cut columns that aren't signatures
VGSig <- VGSig[1:50,c(2:9, 11:13)]

#Convert each specific column to element of list
VGSig.list <- list()
for (i in 1:length(VGSig[1,])){
  VGSig.list[[i]] <- VGSig[,i]
}
names(VGSig.list) <- colnames(VGSig)

allSeurat <- AddModuleScore(allSeurat, VGSig.list, search = T)
allSeurat$VG_HSC.Prog <- allSeurat$Cluster1
allSeurat$VG_GMP <- allSeurat$Cluster2
allSeurat$VG_Myeloid <- allSeurat$Cluster3
allSeurat$VG_HSC.Prog.like <- allSeurat$Cluster4
allSeurat$VG_GMP.like <- allSeurat$Cluster5
allSeurat$VG_Myeloid.like <- allSeurat$Cluster6
allSeurat$VG_HSC.like <- allSeurat$Cluster7
allSeurat$VG_Prog.like <- allSeurat$Cluster8
allSeurat$VG_ProMono.like <- allSeurat$Cluster9
allSeurat$VG_Mono.like <- allSeurat$Cluster10
allSeurat$VG_cDC.like <- allSeurat$Cluster11

FeaturePlot(allSeurat, features = "VG_GMP.like")

#Import the Karamitros signature
karamitros_sig <- read.csv(paste0(sigDir, 'karamitrosSignatures.csv'))

#Convert each specific column to element of list
karamitros_sig.list <- list()
for (i in 1:length(karamitros_sig[1,])){
  karamitros_sig.list[[i]] <- karamitros_sig[,i]
}
names(karamitros_sig.list) <- colnames(karamitros_sig)

allSeurat <- AddModuleScore(allSeurat, karamitros_sig.list, search = T)
allSeurat$karamitros_GMP <- allSeurat$Cluster1
allSeurat$karamitros_MLP <- allSeurat$Cluster2
allSeurat$karamitros_LMPP <- allSeurat$Cluster3

FeaturePlot(allSeurat, features = "karamitros_LMPP")

#Delete temp metadata that has been renamed
allSeurat$Cluster1 <- NULL
allSeurat$Cluster2 <- NULL
allSeurat$Cluster3 <- NULL
allSeurat$Cluster4 <- NULL
allSeurat$Cluster5 <- NULL
allSeurat$Cluster6 <- NULL
allSeurat$Cluster7 <- NULL
allSeurat$Cluster8 <- NULL
allSeurat$Cluster9 <- NULL
allSeurat$Cluster10 <- NULL
allSeurat$Cluster11 <- NULL


#Import the Eppert gene signature
eppert <- read.csv(paste0(sigDir, 'eppertGeneList.csv'))

allSeurat <- AddModuleScore(allSeurat, list(eppert$x), search = T)
allSeurat$EppertHSC <- allSeurat$Cluster1

FeaturePlot(allSeurat, features = "EppertHSC")

saveRDS(allSeurat, paste0(sigDir, 'allSeurat_39+8_postStandardPipeline_withHarmony_withSeuratGeneScores_05-20-2021.rds'))

#Map values to a separate cluster 2 metadata entry
allSeurat$clus2 <- plyr::mapvalues(
  x = allSeurat$RNA_snn_res.0.05, 
  from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"), 
  to = c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
)

for(i in 18:37){
  scoreName <- names(allSeurat@meta.data)[i]
  combined_df <- data.frame("Cluster2" = allSeurat$clus2, "Score" = allSeurat@meta.data[[i]])
  png(paste0(sigDir, scoreName, '_vlnPlot_with_p_value_mean_sd.png'))
  print(ggplot(combined_df, aes(x=Cluster2, y=Score, group = Cluster2, fill = Cluster2)) + 
    geom_violin() + ggtitle(paste0(scoreName, "Cluster 2")) + stat_compare_means(method = "wilcox.test") + 
    stat_summary(fun.data=mean_sdl, geom="pointrange", color="black"))
  dev.off()
}



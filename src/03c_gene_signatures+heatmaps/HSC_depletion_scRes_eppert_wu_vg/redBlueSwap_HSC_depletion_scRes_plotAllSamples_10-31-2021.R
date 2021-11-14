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

rm(list = ls())

#Save all heatmaps in a separate folder
filedirSave <- '/Users/Brian/scCode/Paper_Figures/'

#Import the Wu Gene List
WuGeneList <- read.csv('/Users/Brian/scCode/geneListsProgenitorCellState_Wu2020BloodAdvances.csv', header = T, skip = 2)[2:101,]
row.names(WuGeneList) <- 1:100
genes_hscWu <- WuGeneList$HSC[1:50]

#Import the van galen gene list
VanGalenGeneList <- read.csv('/Users/Brian/scCode/VanGalen_mmc3_GeneListClassifier.csv', header = T, skip = 1)
genes_hscVG <- VanGalenGeneList$HSC.Prog[1:50]

#Try the nature paper hsc list 
EppertGeneList <- read.xls('/Users/Brian/scCode/eppert/natureEppert_HSC+LSC_signatures.xls', header = T, skip = 1, sheet = 8)
genes_hscEp <- EppertGeneList$Gene.Symbol[c(1:53, 63:71)]
#Fix a few genes which have synonyms
genes_hscEp[2] <- 'ABCB4'
genes_hscEp[19] <- 'CFH'
genes_hscEp[48] <- 'GNL1'

#Combine vectors
genes_hscCombined <- c(genes_hscWu, genes_hscVG, genes_hscEp)
allGenes <- genes_hscCombined

#get the last gene in the VanGalen and Wu signatures so we know where to draw the line
lastWu <- genes_hscWu[50]
lastVG <- c(genes_hscWu, genes_hscVG)[length(c(genes_hscWu, genes_hscVG))]


#Import the seurat object with all of our samples assigned to clusters. Scale the data
allSeurat <- readRDS('/Users/Brian/Downloads/seuratObj_05-11-2021_postStandardPipeline_withHarmony_allSamples_39+8.rds')
allSeurat <- ScaleData(allSeurat, features = allGenes)

#Make vector of all the names of each sample
orig.ident.vec <- unique(allSeurat$orig.ident)

#Take all normals as the normal example
normals <- subset(allSeurat, orig.ident %in% c("HuaPt1", "HuaPt2", "HuaPt3", "HuaPt4", "SettyPt1",
                                               "SettyPt2", "SettyPt3", "CD34"), features = allGenes)

for (j in 1:length(orig.ident.vec)){

  cmml <- subset(allSeurat, orig.ident == orig.ident.vec[j], features = allGenes)
  
  #Make into dataframes (separate by gene list then combine)
  cmmlWu.df <- data.frame(cmml@assays[["RNA"]]@data[genes_hscWu,])[1:50,]
  genesKeep_VG <- genes_hscVG[genes_hscVG %in% cmml@assays[['RNA']]@data@Dimnames[[1]]]
  cmmlVG.df <- data.frame(cmml@assays[["RNA"]]@data[genesKeep_VG,])
  genesKeep_Ep <- genes_hscEp[genes_hscEp %in% cmml@assays[['RNA']]@data@Dimnames[[1]]]
  cmmlEp.df <- data.frame(cmml@assays[["RNA"]]@data[genesKeep_Ep,])[1:50,]
  
  normalsWu.df <- data.frame(normals@assays[['RNA']]@data[genes_hscWu,])[1:50,]
  normalsVG.df <- data.frame(normals@assays[['RNA']]@data[genesKeep_VG,])
  normalsEp.df <- data.frame(normals@assays[['RNA']]@data[genesKeep_Ep,])[1:50,]
  
  cmml.df <- do.call("rbind", list(cmmlWu.df, cmmlVG.df, cmmlEp.df))
  normals.df <- do.call("rbind", list(normalsWu.df, normalsVG.df, normalsEp.df))
  
  #Downsample the number of cells for equal comparison/easier plotting
  n_samplesCMML <- 858
  n_samplesNormal <- 1000
  cmmlDS <- cmml.df[,sample(ncol(cmml.df), n_samplesCMML)]
  normalsDS <- normals.df[,sample(ncol(normals.df), n_samplesNormal)]
  dfAll <- cbind(cmmlDS, normalsDS)
  
  #Scale by z-score of the included samples, set max at +/- 3 for color scale
  minV <- -3
  maxV <- 3
  for (i in 1:dim(dfAll)[[1]]){
    dfAll[i,] <- (dfAll[i,]-rowMeans(dfAll)[i])/sd(dfAll[i,])
    dfAll[i,] <- sapply(dfAll[i,], function(y) min(max(y,minV),maxV))
  }
  
  #convert to matrix
  matrixNormal <- as.matrix(dfAll)
  
  #Find where to draw lines based on last gene in Van Galen and Wu
  genesKeep <- allGenes[allGenes %in% row.names(cmml.df)]
  colBreak1 <- which(allGenes %in% lastWu)
  colBreak2 <- which(allGenes %in% lastVG)
  
  #Define colormap 
  coul <- colorRampPalette(c('blue','white', 'red'))
  
  #Plot and save
  png(paste0('/Users/Brian/scCode/Paper_Figures/redBlueSwap_hsc_depletion_allSamples_relativeToAllNormals/hscSignatureHeatmap_', orig.ident.vec[j], '_normalAll_withEppertSignature_allowDuplicates.png'))#, width = 8, height = 6)
  heatmap.2(t(matrixNormal), rowsep = n_samplesCMML,
            colsep = c(colBreak1, colBreak2),
            sepwidth=c(.3,3),
            sepcolor="black", trace="none",
            Rowv=F,Colv=F, scale="none", dendrogram="none",key=F, 
            lhei = c(0.01,.5),margins=c(1,8), col = coul)
  dev.off()
}

library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
library(patchwork)
library(BiocNeighbors)
library(uwot)
library(data.table)
library(matrixStats)
library(FactoClass)
library(gplots)

#Perform a pseudo bulk analysis of the samples, then run UMAP
#and hierarchical clustering

#Clear workspace
rm(list=ls())

#Set file directory to pull from and open seurat object 
## with all the samples (already log-normalized)
filedirOpen <- '/Users/Brian/Downloads/'
allSeurat <- readRDS(paste0(filedirOpen, 'seuratObj_05-11-2021_postStandardPipeline_withHarmony_allSamples_39+8.rds'))
allSeurat@assays[['RNA']]@data['CD38',]

#Import gene list and keep first 60 genes from each type
WuGeneList <- read.csv('/Users/Brian/scCode/geneListsProgenitorCellState_Wu2020BloodAdvances.csv', header = T, skip = 2)[2:101,]
row.names(WuGeneList) <- 1:100
genes <- c(WuGeneList$HSC[1:60], WuGeneList$GMP[1:60], WuGeneList$MEP[1:60])

allSeurat <- ScaleData(allSeurat, features=genes)

#Split by sample ("Identity")
all.list <- SplitObject(allSeurat, split.by = "orig.ident")

#Pre-allocate matrix: rows are samples, columns are genes (avg. gene expression)
bulkMatrix <- matrix(nrow = 47, ncol = 180)

#Fill bulkMatrix using scaled average expression
count <- 1
for (object in all.list){
  #Get scaled data only for variable features
  data <- object@assays[["RNA"]]@scale.data[genes,]
  avgGeneExp <- rowMeans(data)
  avgGeneT <- t(avgGeneExp)
  #Take the avg of every row (gene)
  #if (count < 9){
  #  bulkMatrix[count+39,] <- avgGeneT
  #}else{
  #  bulkMatrix[count-8,] <- avgGeneT
  #}
  
  bulkMatrix[count,] <- avgGeneT
  count <- count + 1
}

#row.names(bulkMatrix) <- names(all.list)
#row.names(bulkMatrix) <- c(40:47, 1:22, 25, 23:24, 26:39)
bulkMatrix <- bulkMatrix[c(9:30, 32:33, 31, 34:47, 1:8),]
names(all.list)

colnames(bulkMatrix) <- genes
#Save files in this directory
filedir <- "/Users/Brian/scCode/Paper_Figures/"
coul <- colorRampPalette(c('blue', 'white', 'red'))

group1 <- c(9,14,15,18,23,25,27,31,34)
group3 <- c(1,5,10,11,19,24,28,29,30,36,38)
group4 <- c(2, 3,4,6, 7,8, 12,13,16,17, 20, 21,22,26,32,33,35,37,39)
normals <- c(40,41,42,43,44,45,46, 47)
orderList <- c(group3, group1, group4, normals)

bulkMatrixOrdered <- bulkMatrix[orderList,]

GMP.sums <- rowSums(bulkMatrixOrdered[,1:60])
MEP.sums <- rowSums(bulkMatrixOrdered[,61:120])
HSC.sums <- rowSums(bulkMatrixOrdered[,121:180])

length(c(group1, group3, group4, normals))
pdf('/Users/Brian/scCode/Paper_Figures/redBlueSwap_WuGeneHeatmap_UMAPGroups_noDendrogram_60Genes_10-31-2021.pdf', width = 12, height = 6)
heatmap.2(bulkMatrixOrdered, rowsep=c(11,20,39), colsep=c(60,120),
          #sepwidth=c(10,.1),
          sepcolor="black", trace="none",
          distfun = function(x) dist(x, method = "euclidean"),
          hclustfun = function(x) hclust(x, method="ward.D2"),
          srtCol = 70,
          Rowv=F,Colv=F, scale="none", dendrogram="none",key=F, 
          lhei = c(0.01,.5),margins=c(5,8), col = coul)
dev.off()

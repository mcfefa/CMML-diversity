library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
library(patchwork)
library(BiocNeighbors)
library(uwot)
library(data.table)
library(matrixStats)

#Perform a pseudo bulk analysis of the samples, then run UMAP
#and hierarchical clustering

#Clear workspace
rm(list=ls())

#Set file directory to pull from and open seurat object 
## with all the samples (already log-normalized)
filedirOpen <- '/Users/4472241/scCode/correctedFullIntegration/'
allSeurat <- readRDS(paste0(filedirOpen, 'integratedObject_Louvainclustering_11162020.rds'))

#Import gene list and keep first 15 genes from each type
WuGeneList <- read.csv('/Users/4472241/scCode/geneListsProgenitorCellState_Wu2020BloodAdvances.csv', header = T, skip = 2)[2:101,]
row.names(WuGeneList) <- 1:100
genes <- c(WuGeneList$GMP[1:30], WuGeneList$MEP[1:30], WuGeneList$HSC[1:30])

#Split by sample ("Identity")
all.list <- SplitObject(allSeurat, split.by = "Identity")

#Pre-allocate matrix: rows are samples, columns are genes (avg. gene expression)
bulkMatrix <- matrix(nrow = 46, ncol = 90)

#Fill bulkMatrix using scaled average expression
count <- 1
for (object in all.list){
  #Get scaled data only for variable features
  data <- object@assays[["RNA"]]@data[genes,]
  avgGeneExp <- rowMeans(data)
  avgGeneT <- t(avgGeneExp)
  #Take the avg of every row (gene)
  if (count < 8){
    bulkMatrix[count+39,] <- avgGeneT
  }else{
    bulkMatrix[count-7,] <- avgGeneT
  }
  count <- count + 1
}

#Define color based on treatment, condition
hmaList <- data.matrix(c(5, 6, 7, 17, 32, 33, 35, 37, 39)) + 7
RuxList <- data.matrix(c(34, 38)) +7
color <- matrix(nrow = 46, ncol = 1)
count1 <- 1
count2 <- 1
for (i in 1:46){
  
  if (i > 7){
    color[i-7,] <- 'CMML'
  } else{
    color[i+39,] <- 'Normal'
  }
  
  #Check for hma and rux
  if (count1 < 10){
    if (i == round(hmaList[count1])){
      color[i-7,] <- 'HMA'
      count1 <- count1 + 1
    }
  }
  if (count2 < 3){
    if (i == round(RuxList[count2])){
      color[i-7,] <- 'Rux'
      count2 <- count2 + 1
    }
  }
  if ((i-7)==36){
    color[i-7,] <- "Chemo"
  }
  if ((i-7)==5){
    color[i-7,] <- "AML+HMA"
  }
}

colnames(bulkMatrix) <- genes
#Save files in this directory
filedir <- "/Users/4472241/scCode/Paper_Figures/"
coul <- colorRampPalette(c('red', 'white', 'blue'))
group1 <- c(9,14,15,18,23,25,27,31,34)
group3 <- c(1,5,10,11,19,24,28,29,30,36,38)
group4 <- c(2,3,4,6,7,8,12,13,16,17,20,21,22,26,32,33,35,37,39)
normals <- c(40,41,42,43,44,45,46)
length(c(group1, group3, group4, normals))
heatmap.2(t(bulkMatrix[c(group1, group3, group4, normals),]), colsep=c(9,20,37), rowsep=c(30,60),
          #sepwidth=c(10,.1),
          sepcolor="black", trace="none",
          Rowv=F,Colv=F, scale="row", dendrogram="none",key=F, 
          lhei = c(0.01,.5),margins=c(1,8), col = coul)

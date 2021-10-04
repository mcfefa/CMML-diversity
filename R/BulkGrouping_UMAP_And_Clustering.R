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
filedirOpen <- '/Users/Brian/Downloads/'
allSeurat <- readRDS(paste0(filedirOpen, 'seuratObj_05-11-2021_postStandardPipeline_withHarmony_allSamples_39+8.rds'))

#Find var. features and scale (2000 is default)
allSeurat <- FindVariableFeatures(allSeurat)
allSeurat <- ScaleData(allSeurat, features = VariableFeatures(allSeurat))

#Split by sample ("Identity")
all.list <- SplitObject(allSeurat, split.by = "orig.ident")

#Pre-allocate matrix: rows are samples, columns are genes (avg. gene expression)
bulkMatrix <- matrix(nrow = 47, ncol = 2000)

#Fill bulkMatrix using scaled average expression
count <- 1
for (object in all.list){
  #Get scaled data only for variable features
  data <- object@assays[["RNA"]]@scale.data
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

row.names(bulkMatrix) <- names(all.list)
row.names(bulkMatrix) <- c(40:47, 1:22, 25, 23:24, 26:39)

#Export so COMET can use it
#counts <- data.frame(t(bulkMatrix))
#row.names(counts) <- VariableFeatures(allSeurat)
#write.table(counts, file = "/Users/4472241/scCode/forCometPseudoBulk/pseudoBulkMarkers.txt", sep = "\t",
#            row.names = TRUE, col.names = NA)

#Define color based on treatment, condition
#hmaList <- data.matrix(c(5, 6, 7, 17, 32, 33, 35, 37, 39)) + 8
#RuxList <- data.matrix(c(34, 38)) +8
#color <- matrix(nrow = 47, ncol = 1)
#count1 <- 1
#count2 <- 1
#for (i in 1:47){
#  
#  if (i > 8){
#    color[i-8,] <- 'CMML'
#  } else{
#    color[i+39,] <- 'Normal'
#  }
#  
##  #Check for hma and rux
#  if (count1 < 10){
#    if (i == round(hmaList[count1])){
#      color[i-8,] <- 'HMA'
#      count1 <- count1 + 1
#    }
#  }
#  if (count2 < 3){
#    if (i == round(RuxList[count2])){
#      color[i-8,] <- 'Rux'
##      count2 <- count2 + 1
#    }
#  }
#  if ((i-8)==36){
#    color[i-8,] <- "Chemo"
#  }
#  if ((i-8)==5){
#    color[i-8,] <- "AML+HMA"
#  }
#}

#Save files in this directory
filedir <- "/Users/Brian/scCode/Paper_Figures/"

#Run UMAP (comments show previous parameters for umap that were saved)
#params for umap, n_neighbors = 6, dist = euclidean, min_dist = 0.4, n_epochs = 500
#params for umap alternate n_neighbors = 5, dist = cosine, min_dist = 0.4 n_epochs = 500
#params for umapBest n_neighbors = 5, dist = cosine, min_dist = 0.5, n_epochs = 500
#params for UMAPFINAL n_neighbors = 12, metric = euclidean, min_dist = 0.2, n_epochs = 500
#params for UMAPGlobal n_neighbors = 39 metric = euclidean min_dist = 0.03 n_epochs = 500
#params for UMAPLastTime n_neighbors = 39, metric = euclidean, min_dist = 0.05, n_epochs = 500

#Run UMAP
umapIt <- umap(bulkMatrix, n_neighbors = 39, metric = 'euclidean', min_dist = 0.05, n_epochs = 500, scale = 'none', ret_model = TRUE)
color <- c(rep("normal", 8), rep("CMML", 39))

#Make embedding into dataframe
umapEmbedding <- data.frame("X" = umapIt$embedding[,1], "Y" = umapIt$embedding[,2], 'Color' = color)
#umapEmbedding <- data.frame("X" = umapIt$embedding[,1], "Y" = umapIt$embedding[,2])
row.names(umapEmbedding) <- row.names(bulkMatrix)
#Or load embeddings that we have saved
write.csv(umapEmbedding, "/Users/Brian/scCode/Paper_Figures/pseudoBulkUMAPLastTime.csv")
#umapEmbedding <- read.csv("/Users/Brian/scCode/Paper_Figures/pseudoBulkUMAPLastTime.csv")

#plot all (optionally save)
pdf('/Users/Brian/scCode/Paper_Figures/pseudoBulkUMAPLastTimeNoTx.pdf')
ggplot(umapEmbedding, aes(x=X, y = Y, color = Color)) + geom_point(size = 7) + geom_text(label=rownames(umapEmbedding),
                                                                                         nudge_x = 0.0, nudge_y = 0.0, color = 'black', check_overlap = F)+
  theme(panel.background = element_rect(fill = "white")) +
  scale_color_manual(values=c("light blue", "gray"))
dev.off()

#subset dataframe so we just have normal and cmml
embeddingNoTx <- subset(umapEmbedding, umapEmbedding$Color %in% c('CMML', 'Normal'))

#Plot subsetted data (normal and No Tx only) and optionally save
#pdf('/Users/Brian/scCode/Paper_Figures/pseudoBulkUMAPLastTimeNoTx.pdf')
#ggplot(embeddingNoTx, aes(x=X, y = Y, color = Color)) + geom_point(size = 7) + geom_text(label=rownames(embeddingNoTx),
#                                                                                         nudge_x = 0.0, nudge_y = 0.0, color = 'black', check_overlap = F)+
#  theme(panel.background = element_rect(fill = "white"))+
#  scale_color_manual(values=c( "light blue",  "gray"))
#dev.off()

#Optionally Write umap embeddings to csv 
#write.csv(umapEmbedding, "/Users/4472241/scCode/Paper_Figures/dataMain/pseudoBulkUMAPLastTime.csv")


#*********Cluster using Ward Hierarchical Clustering**************
d <- dist(bulkMatrix, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward.D2")

#Plot and save
pdf(paste0(filedir, 'hierarchicalClustering.pdf'), width = 12, height = 6)
plot(fit) # display dendogram
groups <- cutree(fit, k=4) # cut tree into 4 clusters
# draw dendogram with red borders around the 4 clusters
rect.hclust(fit, k=4, border="red")
dev.off()

#Save groupings as csv
write.csv(groups, file = paste0(filedir, "hierarchicalClusterGroupings.csv"))

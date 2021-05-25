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

#Perform a pseudo bulk analysis of the samples, then run UMAP
#and hierarchical clustering

#Clear workspace
rm(list=ls())

#Load the processed seurat object
dirOpen <- '/Users/4472241/scCode/redoEntirePipeline_05-05-2021/standard_pipeline_withHarmony/'
allSeurat <- readRDS(paste0(dirOpen, 'seuratObj_05-11-2021_postStandardPipeline_withHarmony_allSamples_39+8.rds'))


#Import gene list and keep first 60 genes from each type
WuGeneList <- read.csv('/Users/4472241/scCode/reRunScores_seurat/Wu_Signatures.csv', header = T, skip = 2)[2:101,]
row.names(WuGeneList) <- 1:100
genes <- c(WuGeneList$HSC[1:60], WuGeneList$GMP[1:60], 'CD38', WuGeneList$MEP[1:60])

allSeurat <- ScaleData(allSeurat, features=genes)

#Split by sample ("Identity")
all.list <- SplitObject(allSeurat, split.by = "orig.ident")
names(all.list)[4] <- "Zheng"

#Pre-allocate matrix: rows are samples, columns are genes (avg. gene expression)
bulkMatrix <- matrix(nrow = 47, ncol = 181)

#Fill bulkMatrix using scaled average expression
count <- 1
for (object in all.list){
  #Get scaled data only for features in Wu signatures
  data <- object@assays[["RNA"]]@scale.data[genes,]
  avgGeneExp <- rowMeans(data)
  avgGeneT <- t(avgGeneExp)
  
  #Take the avg of every row (gene)
  bulkMatrix[count,] <- avgGeneT
  count <- count + 1
}
row.names(bulkMatrix) <- names(all.list)
colnames(bulkMatrix) <- genes

#Save files in this directory
filedir <- "/Users/4472241/scCode/Paper_Figures/"
coul <- colorRampPalette(c('red', 'white', 'blue'))

group1 <- c("SF100109106293", "SF12062800475", "SF14072200012", "2V001", "SF100109101914",
            "6AD001", "SF14031800065", "SF16112300029", "SF14111400033")
group3 <- c("LTB3966", "4K001", "SF100109111451", "SF100109110236", "SF14101000049", 
            "SF12042500035", "SF14050700419", "SF16026800045", "SF16072200003", "SF14061300036", 
            "SF15010200008")
group4 <- c("LTB4121", "LTB6169", "4J003", "4Q001", "5E001", "5H001", "SF14040100158",
            "SF14060200025", "SF13061200056", "4S001", "SF16112900158", "6AE001",    
            "6AC001", "SF12092600014", "SF13032800016", "SF14110400108", "SF14092500135", 
            "SF14080400065", "SF13070900171")
normals <- c("SettyPt1", "SettyPt2", "SettyPt3", "Zheng",     "HuaPt1",   "HuaPt2",
             "HuaPt3",   "HuaPt4" )
orderList <- c(normals, group4, group1, group3)

bulkMatrixOrdered <- bulkMatrix[orderList,]

GMP.sums <- rowSums(bulkMatrixOrdered[,1:45])
MEP.sums <- rowSums(bulkMatrixOrdered[,46:90])
HSC.sums <- rowSums(bulkMatrixOrdered[,91:135])

length(c(group1, group3, group4, normals))
pdf(paste0(filedir, 'WuGeneHeatmap_updated_05-25-2021_UMAPGroups_noDendrogram_60GenesEach.pdf'), width = 9, height = 8)
heatmap.2(t(bulkMatrixOrdered), colsep=c(8,27,36), rowsep=c(60,121),
          #sepwidth=c(10,.1),
          sepcolor="black", trace="none",
          distfun = function(x) dist(x, method = "euclidean"),
          hclustfun = function(x) hclust(x, method="ward.D2"),
          labRow = F,
          srtCol = 80, offsetCol = -.5,
          Rowv=F,Colv=F, scale="none", dendrogram="none",key=F, 
          lhei = c(0.0005,.025),margins=c(8,3), col = coul, 
          colCol = c(rep("black", 8), rep("#D43030", length(group4)),
                     rep("#33CEFF", length(group1)), rep("#7FC93C", length(group3))))
      
 dev.off()

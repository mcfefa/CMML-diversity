library(Seurat)
library(harmony)
library(Matrix)
library(stringr)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(RColorBrewer)
library(gplots)
library(clustree)

rm(list =ls())

#Set dir
dir <- '/Users/4472241/scCode/redoEntirePipeline_05-05-2021/'
outdir <- '/Users/4472241/scCode/redoEntirePipeline_05-05-2021/standard_pipeline_withHarmony/'
date <- '_05-06-2021'

#Read in seurat object with all samples
allSeurat <- readRDS(paste0(dir, 'postQC_CMML39+healthy8_05-06-2021.rds'))

#Log normalize (using default params, but include anyway for clarity)
allSeurat <- NormalizeData(allSeurat, normalization.method = "LogNormalize", scale.factor = 10000)

#Find variable features (using default params, but include for clarity)
allSeurat <- FindVariableFeatures(allSeurat, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes and plot
top10 <- head(VariableFeatures(allSeurat), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(allSeurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Scale data (only use HVG), regressing out effects of nCountRNA and percent.mito
allSeurat <- ScaleData(allSeurat, features = VariableFeatures(allSeurat), vars.to.regress = c("nCount_RNA","percent.mito"))

#Run, plot and save PCA
allSeurat <- RunPCA(allSeurat, features = VariableFeatures(object = allSeurat))
pdf(paste0(outdir, 'PCA_allSeurat_standard_seurat_pipeline_preHarmony', date, '.pdf'),
    width = 10, height = 2)
DimPlot(allSeurat, reduction = "pca", split.by = "tech", pt.size = 0.0001)
dev.off()

#Determine dimensionality of dataset
ElbowPlot(allSeurat, ndims = 50)

allSeurat <- RunHarmony(allSeurat, group.by.vars = "tech")

#Cluster at various resolutions to make clustree plot
allSeurat <- FindNeighbors(allSeurat, dims = 1:50, reduction = "harmony")
resVec <- c(0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2)
for (res in resVec){
  allSeurat <- FindClusters(allSeurat, resolution = res)
}

pdf(paste0(outdir, 'clustree_standard_pipeline_withHarmony.pdf'), width = 9, height = 7)
clustree(allSeurat, prefix = "RNA_snn_res.")
dev.off()

#Run, plot, and save UMAP
resUse <- 0.1

#Fix ordering of resolution we use (default is that clusters 10, 11, come before 2)
allSeurat$clusterResolution_0.1 <- as.factor(as.numeric(as.character(allSeurat$RNA_snn_res.0.1)))

allSeurat <- RunUMAP(allSeurat, dims = 1:50, reduction = "harmony")
pdf(paste0(outdir, 'UMAP_res=', resUse, '_standard_pipeline_withHarmony.pdf'))
DimPlot(allSeurat, reduction = "umap", group.by = paste0("clusterResolution_", resUse))
dev.off()

pdf(paste0(outdir, 'UMAP_res=', resUse, '_standard_pipeline_withHarmony_splitByTech.pdf'),
    width = 10, height = 4)
DimPlot(allSeurat, reduction = "umap", split.by = "tech", group.by = paste0("clusterResolution_", resUse))
dev.off()

############   Export the cluster composition to csv
date <- '05-07-2021'

# create temporary files for storing names and numbers
file1 <- "./outputDataNames_20210506.csv"
file2 <- "./outputData_20210506.csv"

# save the cluster identity for each individual cell
# NOTE: you need to change the $RNA_snn_res.0.05 part if you identify a better resoultion for clustering your data 
cat(allSeurat@meta.data$clusterResolution_0.1, file=file2, sep=",\n")
# save the orig.ident per cell 
listNames <- allSeurat@meta.data$orig.ident
write.table(data.frame(listNames),
            row.names=FALSE,
            col.names = FALSE, 
            file = file1,
            sep=",")
# merging the two files together in grouping to ultimately write out the number of cells per identity per cluster 
#(so in this case, we'd end up with A being a table that has 5 samples x number of cluster idenfied number of rows and two columns )
mydat1 <- read.csv(file2)
mydat2 <- read.csv(file1)
fulldat <- cbind(mydat2[1],mydat1[1])
fulltab <- as.data.table(fulldat)
# name table columns
names(fulltab)[1] <- paste("UMI")
names(fulltab)[2] <- paste("cluster")
# group data based on clusters
group_by(fulltab, cluster)
# create a table counting unqiue UMIs/cells per cluster
tabPerClus <- fulltab %>% group_by(cluster) %>% count()
type <- sub("\\_.*","",fulltab$UMI)
fulltab <- cbind(fulltab, type)
A <- fulltab %>% group_by(cluster) %>% count(type)

#Order A and sum same components
A <- A[order(A$type),]
count <- 0
for (i in 1:length(A$type)){
  if ( i == 1){
    A[i,'PartialSum'] <- A[i,'n']
    count <- count+1
  }else if (A[i, 'type'] %in% A[i-1, 'type']){
    A[i, 'PartialSum'] <- A[i-1, 'PartialSum'] + A[i, 'n']
    count <- count + 1
    if (i == length(A$type)){
      A[c((i-count):(i)), 'Sum'] <- A[i, 'PartialSum']
      break
    }
  }else{
    A[i, 'PartialSum'] <- A[i, 'n']
    if (i-1-count > 0){
      A[c((i-1-count):(i-1)), 'Sum'] <- A[i-1, 'PartialSum']
    }else{
      A[c((i-count):(i-1)), 'Sum'] <- A[i-1, 'PartialSum']
    }
    count <- 0
  }
}

A$Fraction <- A$n/A$Sum


# saving matrix/table A as a CSV file that will later be read into Julia for diversity calculations
divout <- paste(outdir,"CellBreakdown_PerClusterPerType_standard_pipeline_withHarmony_res-", resUse, "_",date,".csv",sep="") 
write.csv(A, file=divout)




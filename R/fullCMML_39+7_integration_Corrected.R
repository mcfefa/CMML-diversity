#Here we just check the projection of query pca in reference pca space
library(Seurat)
library(dplyr)
library(SeuratData)
library(harmony)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggplot2)
library(harmony)
library(cowplot)
library(BiocNeighbors)
library(data.table)
library(leiden)
library(reticulate)
library(igraph)
library(matrixStats)

#Clear workspace var.
rm(list = ls())

#Set dir to save to
filedir <- '/Users/4472241/scCode/correctedFullIntegration/'

#Start the clock
ptm <- proc.time()

#Define function to scale query seurat objects and run pca in reference space
ScaleAndProjectPCA <- function(newObject, refObject, scale_max = 10){
  
  # inspired by https://www.r-bloggers.com/a-faster-scale-function/
  refMat <- data.matrix(refObject@assays[['RNA']]@data[VariableFeatures(refObject),])
  mat <- data.matrix(newObject@assays[['RNA']]@data[VariableFeatures(refObject),])
  
  #Find reference mean
  rmRef <- rowMeans2(x = refMat, na.rm = TRUE)
  
  #Find Reference SD
  rsdRef <- rowSds(refMat, center = rmRef)
  
  #Center
  mat <- mat - rmRef
  
  #Scale
  mat <- mat / rsdRef
  
  #Set Scaled Max to adjusted data
  if (scale_max != Inf) {
    mat[mat > scale_max] <- scale_max
  }
  
  #Transpose
  B <- t(mat)
  
  #Find basis vector
  basisVec <- data.matrix(refObject@reductions[["pca"]]@feature.loadings)
  
  #T.matrix is your embedding result
  T.matrix <- B %*% basisVec
  
  #add row (cell) names
  row.names(T.matrix) <- newObject@assays[["RNA"]]@counts@Dimnames[[2]]
  
  #Add reduction to seurat object
  newObject[["pca"]] <- CreateDimReducObject(embeddings = T.matrix, key = "PCA_", assay = DefaultAssay(newObject))
  
  #return Seurat object in reductions->pca
  return(newObject)
}

#Load seven normal patient data
#sevenHealthy <- readRDS("/Users/4472241/scCode/normalData/Normal_CD34only_Zheng+Setty+Hua_QCprocessedCohort_SeuratObj_20200321.rds")

#Specify the seven Healthy object as reference in metadata
#sevenHealthy[["Integration"]] <- "Reference"

#merge, normalize, runPCA on reference cohort (7 normal)
#sevenHealthy <- NormalizeData(sevenHealthy)
#sevenHealthy <- FindVariableFeatures(sevenHealthy, selection.method = "vst", nfeatures = 2000)
#all.genes <- rownames(sevenHealthy)
#sevenHealthy <- ScaleData(sevenHealthy, features = all.genes)
#sevenHealthy <- RunPCA(sevenHealthy, features = VariableFeatures(object = sevenHealthy))

#RunHarmony and UMAP and plot before and after
#DimPlot(sevenHealthy, reduction = 'pca', group.by = "Identity")
#sevenHealthy <- RunHarmony(sevenHealthy, group.by.vars = c("Identity", "Tech"), theta = c(2,0),
                       #reference_values = "Setty1", 
                       #max.iter.harmony = 20, plot_convergence = TRUE)
#DimPlot(sevenHealthy, reduction = 'harmony', group.by = "Identity")
#sevenHealthy <- RunUMAP(sevenHealthy, reduction = 'harmony', dims = 1:10)
#DimPlot(sevenHealthy, reduction = "umap", group.by = "Identity")

#Load/save rds file with all the integrated info from the 7 normal samples
#saveRDS(sevenHealthy, "/Users/4472241/scCode/normalData/sevenHealthyHarmonizedForIntegration.rds")
sevenHealthy <- readRDS('/Users/4472241/scCode/normalData/sevenHealthyHarmonizedForIntegration.rds')

all.genes <- rownames(sevenHealthy)

#*****Normal 7 have been harmonized above, now layer on CMML******
integratedObject <- sevenHealthy
integratedObject@reductions[["postHarmony"]] <- integratedObject@reductions[["harmony"]]
harmonyEmbeddings.matrix <- data.matrix(integratedObject@reductions[["postHarmony"]]@cell.embeddings)

#Once we have harmonized...iteratively layer queries (CMML) on original pcaspace
rdsFileList <- c('CMML1_LTB3966_2020-10-15', 'CMML2_LTB4121_2020-10-15', 'CMML3_LTB5109_2020-10-15', 'CMML4_4-J-003_postBMT_2020-10-15',
                 'CMML5_4-K-001_HMA_2020-10-15', 'CMML6_4-Q-001_HMA_2020-10-15', 'CMML7_5-E-001_HMA_2020-10-15',
                 'CMML8_5-H-001_2020-10-15', 'CMML9_SF-100109-106293_2020-10-15', 'CMML10_SF-100109-111451_2020-10-15',
                 'CMML11_SF-100109-110236_2020-10-15', 'CMML12_SF-140401-00158_2020-10-15', 'CMML13_SF-140602-00025_2020-10-15',
                 'CMML14_SF-120628-00475_2020-10-15', 'CMML15_SF-140722-00012_2020-10-15', 'CMML16_SF-130612-00056_2020-10-15',
                 'CMML17_4-S-001_HMA_2020-10-15','CMML18_2-V-001_2020-10-15', 'CMML19_SF-141010-00049_2020-10-15',
                 'CMML20_SF-161129-00158_2020-10-15', 'CMML21_6-AE-001_2020-10-15', 'CMML22_6-AC-001_2020-10-15',
                 'CMML24_SF-100109-101914_2020-10-15','CMML25_SF-120425-00035_2020-10-15', 'CMML23+26_6-AD-001_2020-10-15',
                 'CMML27_SF-120926-00014_2020-10-15','CMML28_SF-140318-00065_2020-10-15', 'CMML29_SF-140507-00419_2020-10-15',
                 'CMML30_SF-160268-00045_2020-10-15', 'CMML31_SF-160722-00003_2020-10-15', 'CMML32_SF-161123-00029_2020-10-15',
                 'CMML33_SF-130328-00016_2020-10-15', 'CMML34_SF-141104-00108_2020-10-15', 'CMML35_SF-141114-00033_2020-10-15',
                 'CMML36_SF-140925-00135_2020-10-15', 'CMML37_SF-140613-00036_2020-10-15', 'CMML38_SF-140804-00065_2020-10-15',
                 'CMML39_SF-150102-00008_2020-10-15', 'CMML40_SF-130709-00171_2020-10-15')

for (i in 1:length(rdsFileList)){
  
  #Define the individual seurat object for the sample currently being layered on
  cmml <- readRDS(paste('/Users/4472241/scCode/', rdsFileList[[i]], '.rds', sep=''))
  
  #Pre-process
  mito.genes <- grep(pattern = "^MT-", cmml@assays$RNA@counts@Dimnames[[1]], value = TRUE)
  percent.mito <- Matrix::colSums(x=GetAssayData(object=cmml, slot='counts')[mito.genes, ])/Matrix::colSums(x=GetAssayData(object=cmml, slot='counts'))
  cmml[['percent.mito']] <- percent.mito
  
  nFeatUpperN <- 5808.621 #Set to be consistent across all samples
  nFeatLowerN <- 200
  perMitoUpperN <- 0.25 
  
  # We filter out cells that have unique feature counts over 5808.621 or less than 200
  cmml <- subset(x=cmml, subset=nFeature_RNA > nFeatLowerN & nFeature_RNA < nFeatUpperN & percent.mito < perMitoUpperN)
  
  #Normalize and specify in metadata which individual
  cmml <- NormalizeData(cmml)
  cmml[["Identity"]] <- paste("cmml", i, sep='')
  
  #Scale according to reference (sevenHealthy) and project into pca space 
  cmml <- ScaleAndProjectPCA(cmml, sevenHealthy)
  
  #Merge data together with normal 7
  eightObject <- merge(sevenHealthy, cmml, merge.data = TRUE)
  
  #Get pca embeddings for cmml object 
  currentQueryEmbeddings <- cmml@reductions[['pca']]@cell.embeddings
    
  #Combine pca embeddings to input to harmony
  harmonyInput <- rbind(data.matrix(sevenHealthy@reductions[["harmony"]]@cell.embeddings), currentQueryEmbeddings)
  
  #Set metadata as query to differentiate between reference
  cmml[["Integration"]] <- "Query"
  
  #Run Harmony
  integratedEmbeddings <- HarmonyMatrix(harmonyInput,  eightObject@meta.data, "Integration", do_pca=FALSE,
                                        reference_values = "Reference",
                                        max.iter.harmony = 20,
                                        theta=2)
  
  #Merge cmml layered on with integratedObject (which will contain all samples)
  integratedObject <- merge(integratedObject, cmml, merge.data = TRUE)
  
  #Get only the new embeddings produced by running harmony
  onlyNewEmbeddings <- integratedEmbeddings[(sevenHealthy@assays[["RNA"]]@counts@Dim[[2]]+1):dim(integratedEmbeddings)[[1]],]
  
  #Add the new harmony embeddings from the current cmml sample to matrix containing all new harmony embeddings
  harmonyEmbeddings.matrix <- rbind(harmonyEmbeddings.matrix, onlyNewEmbeddings)
}

#******Outside for loop*******
#Add the harmony reduction of all samples to the integrated object
integratedObject[["postHarmony"]] <- CreateDimReducObject(embeddings = harmonyEmbeddings.matrix, 
                                                          key = "postHarm_", assay = DefaultAssay(integratedObject))

#Plot harmony and save
pdf(paste0(filedir, 'postHarmony39+7Integrated.pdf'), width = 12, height = 6)
print(DimPlot(integratedObject, reduction = 'postHarmony', group.by = c("Integration", 'Identity')))
dev.off()

#Add the feature.loadings (basis vector) to reduction for integrated object
integratedObject@reductions[["postHarmony"]]@feature.loadings <- sevenHealthy@reductions[["pca"]]@feature.loadings

#Run UMAP on integrated object
integratedObject <- RunUMAP(integratedObject, reduction = 'postHarmony', dims = 1:30)

#*******Optionally run using Leiden, though louvain is much faster in our implementation
#Install the R package for using leiden
#devtools::install_github("TomKellyGenetics/leiden")

#Compute SNN and cluster using louvain
integratedObject <- FindNeighbors(integratedObject, reduction = "postHarmony", dims = 1:50) %>% 
  FindClusters(resolution = 0.05)

#Install R accessible leidenalg and igraph
#reticulate::use_condaenv("/Users/4472241/anaconda3/envs/build-doc/", required = TRUE)

#Use SNN as input to leiden...write to csv for reading in python
#writeMM(data.frame(integratedObject@graphs[["RNA_snn"]]), '/Users/4472241/scCode/sparseAdjMatIntegratedObject_10_20_2020.txt')

#Use weighted adjacency matrix
#adjacencyMatrix <- integratedObject@graphs$RNA_snn
#snnGraph <- graph_from_adjacency_matrix(adjacencyMatrix, weighted = TRUE)
#Cluster_assignment <- leiden(snnGraph, resolution_parameter = 0.05, n_iterations=-1)
#integratedObject@meta.data$Leiden_assigned_clusters_res_0.05 <- Cluster_assignment

#Plot and save umap three different visualizations
pdf(paste0(filedir, '39+7IntegratedUMAP.pdf'), width = 12, height = 6)
p1 <- DimPlot(integratedObject, group.by = "Identity", reduction = "umap", label = FALSE) + NoLegend()
p2 <- DimPlot(integratedObject, group.by = "Integration", reduction = "umap") 
p3 <- DimPlot(integratedObject, group.by = "RNA_snn_res.0.05", reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()
print(p1)
print(p2)
print(p3)
dev.off()

#Save rds of fully integrated object using harmony
saveRDS(integratedObject, paste0(filedir,'fullHarmonyIntegratedProcessedObject.rds'))
#integratedObject <- readRDS('/Users/4472241/scCode/correctedFullIntegration/fullHarmonyIntegratedProcessedObject.rds')



#*********************EXPORT TO JULIA FOR DIVERSITY/HEATMAPS********************

# create temporary files for storing names and numbers
file1 <- "/Users/4472241/scCode/outputDataNames_20201111.csv"
file2 <- "/Users/4472241/scCode/outputData_20200722.csv"
filedir='/Users/4472241/scCode/correctedFullIntegration/'
date = paste0('_', Sys.Date(), '_')

# save the cluster identity for each individual cell
# note: you need to change the $RNA_snn_res.0.6 part if you identify a better resoultion for clustering your data 
cat(integratedObject@meta.data$RNA_snn_res.0.05, file=file2, sep=",\n")
# save the orig.ident per cell 
listNames <- integratedObject@meta.data$orig.ident
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
tabPerClus <- fulltab %>% group_by(cluster) %>% dplyr::count()
type <- sub("\\_.*","",fulltab$UMI)
fulltab <- cbind(fulltab, type)
A <- fulltab %>% group_by(cluster) %>% dplyr::count(type)
# saving matrix/table A as a CSV file that will later be read into Julia for diversity calculations
divout <- paste(filedir,"FullHarmonyIntegration_CellBreakdown_PerClusterPerType_res-0.05",date,".csv",sep="") 
write.csv(A, file=divout)

#Print the time it took to complete
print(proc.time()-ptm)

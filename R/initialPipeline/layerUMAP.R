library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
library(patchwork)
library(BiocNeighbors)
library(uwot)
library(data.table)
library(matrixStats)

#Clear workspace
rm(list=ls())

#Set directory to save to
filedir <- "/Users/4472241/scCode/finalLayerUMAP_Code+Output/allSamples/"
tag <- '_noHarmony_30Neighbors_AllSamples'

#Define function used to scale and comput pca
scaleRegressAndProjectPCA <- function(newObject, refObject, scale_max = 10, 
                                      vars.to.regress = c("nCount_RNA","percent.mito"),
                                      block.size = 1000, split.by = NULL,
                                      model.use = 'linear', use.umi = FALSE,
                                      do.scale = TRUE, do.center = TRUE){
  #Enter default settings for scaledata in seurat
  split.by = NULL
  model.use = 'linear'
  use.umi = FALSE
  do.scale = TRUE
  do.center = TRUE
  verbose = TRUE
  
  #Set latent.data to NULL initially
  latent.dataNew <- NULL
  latent.dataRef <- NULL
  
  #Temporary info
  #refObject <- oldSamples
  #newObject <- newSamples
  #vars.to.regress <- c("nCount_RNA","percent.mito")
  
  #Define chunk points function
  ChunkPoints <- function(dsize, csize) {
    return(vapply(
      X = 1L:ceiling(x = dsize / csize),
      FUN = function(i) {
        return(c(
          start = (csize * (i - 1L)) + 1L,
          end = min(csize * i, dsize)
        ))
      },
      FUN.VALUE = numeric(length = 2L)
    ))
  }
  
  RegressOutUsingRef <- function(
    data.exprRef,
    latent.dataRef,
    data.exprNew,
    latent.dataNew,
    features.regress,
    model.use = 'linear',
    use.umi = FALSE
  ) {
    
    # Check features.regress
    if (is.null(x = features.regress)) {
      features.regress <- 1:nrow(x = data.exprNew)
    }
    if (is.character(x = features.regress)) {
      features.regress <- intersect(x = features.regress, y = rownames(x = data.exprNew))
      if (length(x = features.regress) == 0) {
        stop("Cannot use features that are beyond the scope of data.expr")
      }
    } else if (max(features.regress) > nrow(x = data.exprNew)) {
      stop("Cannot use features that are beyond the scope of data.expr")
    }
    # Check dataset dimensions
    if (nrow(x = latent.dataNew) != ncol(x = data.exprNew)) {
      stop("Uneven number of cells between latent data and expression data")
    }
    
    # Create formula for regression using reference data
    vars.to.regress <- colnames(x = latent.dataRef)
    fmla <- paste('GENE ~', paste(vars.to.regress, collapse = '+'))
    fmla <- as.formula(object = fmla)
    
    # Make results matrix
    data.resid <- matrix(
      nrow = nrow(x = data.exprNew),
      ncol = ncol(x = data.exprNew)
    )
    
    #Define progress bar
    pb <- txtProgressBar(char = '=', style = 3, file = stderr())
    
    #Regress out using reference formula
    for (i in 1:length(x = features.regress)) {
      
      x <- features.regress[i]
      regression.matNew <- cbind(latent.dataNew, data.exprNew[x, ])
      if (model.use == "linear") {
        # Repeatedly call lm to fit model using reference, and then apply it to
        #the new cells
        regression.matRef <- cbind(latent.dataRef, data.exprRef[x,])
        colnames(x = regression.matRef) <- c(vars.to.regress, 'GENE')
        qr <- lm(fmla, data = regression.matRef, qr = TRUE)
        rm(regression.matRef)
      }
      colnames(x = regression.matNew) <- c(vars.to.regress, 'GENE')
      regression.matNew <- -qr[["coefficients"]][["(Intercept)"]]+data.exprNew[x,]-
        qr[["coefficients"]][["nCount_RNA"]]*latent.dataNew[,"nCount_RNA"]-
        qr[["coefficients"]][["percent.mito"]]*latent.dataNew[,"percent.mito"]
      data.resid[i, ] <- regression.matNew
      setTxtProgressBar(pb = pb, value = i / length(x = features.regress))
    }
    
    dimnames(x = data.resid) <- dimnames(x = data.exprNew)
    return(data.resid)
  }
  
  #****************FROM ScaleData.Seurat********************************
  ###Fill latent.data for New and Ref Object
  #New
  assayNew <- DefaultAssay(object = newObject)
  assay.dataNew <- GetAssayData(object = newObject, assay = assayNew)
  if (any(vars.to.regress %in% colnames(x = newObject[[]]))) {
    latent.dataNew <- newObject[[vars.to.regress[vars.to.regress %in% colnames(x = newObject[[]])]]]
  } else {
    latent.dataNew <- NULL
  }
  #Ref
  assayRef <- DefaultAssay(object = refObject)
  assay.dataRef <- GetAssayData(object = refObject, assay = assayRef)
  if (any(vars.to.regress %in% colnames(x = refObject[[]]))) {
    latent.dataRef <- refObject[[vars.to.regress[vars.to.regress %in% colnames(x = refObject[[]])]]]
  } else {
    latent.dataRef <- NULL
  }
  
  #***********FROM ScaleData.Default********************************
  #Keep only the variable features of reference
  features <- VariableFeatures(refObject)
  newObject.data <- assay.dataNew[features, , drop=FALSE]
  refObject.data <- assay.dataRef[features, , drop =FALSE]
  newObject.names <- dimnames(x = newObject.data)
  refObject.names <- dimnames(x = refObject.data)
  
  #If we have variables to regress
  #Order and assign cell names to rows
  latent.dataRef <- latent.dataRef[colnames(x = refObject.data), , drop = FALSE]
  latent.dataNew <- latent.dataNew[colnames(x = newObject.data), , drop = FALSE]
  rownames(x = latent.dataNew) <- colnames(x = newObject.data)
  rownames(x = latent.dataRef) <- colnames(x = refObject.data)
  
  #Update status by printing message
  message("Regressing out ", paste(vars.to.regress, collapse = ', '))
  
  #chunk.points is the same for both ref and new
  chunk.points <- ChunkPoints(dsize = nrow(x = newObject.data), csize = block.size)
  
  #split.by is true when we only regress.by.vars meaning split.cells is all cells
  split.by <- TRUE
  split.cellsRef <- split(x = colnames(x = refObject.data), f = split.by)
  split.cellsNew <- split(x = colnames(x = newObject.data), f = split.by)
  
  newObject.data <- RegressOutUsingRef(
    data.exprRef = refObject.data[, split.cellsRef[["TRUE"]], drop = FALSE],
    latent.dataRef = latent.dataRef[split.cellsRef[['TRUE']], , drop = FALSE],
    data.exprNew = newObject.data[, split.cellsNew[["TRUE"]], drop = FALSE],
    latent.dataNew = latent.dataNew[split.cellsNew[['TRUE']], , drop = FALSE],
    features.regress = features,
    model.use = "linear",
    use.umi = FALSE
  )
  
  refObject.data <- RegressOutUsingRef(
    data.exprRef = refObject.data[, split.cellsRef[["TRUE"]], drop = FALSE],
    latent.dataRef = latent.dataRef[split.cellsRef[['TRUE']], , drop = FALSE],
    data.exprNew = refObject.data[, split.cellsRef[["TRUE"]], drop = FALSE],
    latent.dataNew = latent.dataRef[split.cellsRef[['TRUE']], , drop = FALSE],
    features.regress = features,
    model.use = "linear",
    use.umi = FALSE
  )
  
  dimnames(x = newObject.data) <- newObject.names
  dimnames(x = refObject.data) <- refObject.names
  
  
  
  #**************Regular scale function without vars.to.regress*****************
  # inspired by https://www.r-bloggers.com/a-faster-scale-function/
  refMat <- refObject.data
  mat <- newObject.data
  
  #Find reference mean
  rmRef <- rowMeans2(x = refMat, na.rm = TRUE)
  
  #Find Reference SD
  rsdRef <- rowSds(refMat, center = rmRef)
  
  #Center against reference
  mat <- mat - rmRef
  
  #Scale against reference
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
  
  #return Seurat object with pca reduction included in reductions->pca
  return(newObject)
}

#Import the old samples, 31 patients and 7 healthy controls
oldSamples <- readRDS('/Users/4472241/scCode/Healthy-7pts+CMML-31pts-Cohort_CD34only+33kGenes_SeuratObject_thruStep9_res0.05-as-active_ReScaledAllGenes_2020-04-29.rds')

#Run initial umap, with the pca of old samples as input. Keep params the same as seurat default (note: not using seurat)
early.umap <- umap(oldSamples@reductions[['pca']]@cell.embeddings, n_neighbors = 30, metric = 'cosine', min_dist = 0.3, n_epochs = 200, scale = 'none', ret_model = TRUE)

#Make list of all the new samples to layer on top
newList <- c('CMML33_SF-130328-00016_2020-10-15','CMML34_SF-141104-00108_2020-10-15', 'CMML35_SF-141114-00033_2020-10-15',
'CMML36_SF-140925-00135_2020-10-15', 'CMML37_SF-140613-00036_2020-10-15', 'CMML38_SF-140804-00065_2020-10-15',
'CMML39_SF-150102-00008_2020-10-15', 'CMML40_SF-130709-00171_2020-10-15')

#set cutoffs
nFeatUpperN <- 5808.621 #Set to be consistent with old samples
nFeatLowerN <- 200
perMitoUpperN <- 0.25

#Read in new files, quality control, normalize, and merge with each other
count <- 0
for (i in newList){
  object <- readRDS(paste0('/Users/4472241/scCode/',i,'.rds'))

  #Find mito percentage in each cell
  mito.genes <- grep(pattern = "^MT-", object@assays$RNA@counts@Dimnames[[1]], value = TRUE)
  percent.mito <- Matrix::colSums(x=GetAssayData(object=object, slot='counts')[mito.genes, ])/Matrix::colSums(x=GetAssayData(object=object, slot='counts'))
  object[['percent.mito']] <- percent.mito

  #Filter out those that don't make the cut
  object <- subset(x=object, subset=nFeature_RNA > nFeatLowerN & percent.mito < perMitoUpperN & nFeature_RNA < nFeatUpperN)

  #Normalize
  object <- NormalizeData(object)
  
  #Merge all new samples
  if (count > 0){
    newSamples <- merge(newSamples, object, merge.data = TRUE)
  }
  else{
    newSamples <- object
  }
  count <- count+1
}

#Add metadata differentiating the old from new samples
newSamples[['OldNew']] <- 'New'
oldSamples[['OldNew']] <- 'Old'

#Merge the two seurat objects
all46 <- merge(oldSamples, newSamples)

#Scale and run pca on new samples in space of old samples for consistency
newSamples <- scaleRegressAndProjectPCA(newSamples, oldSamples)

#Correct the umap embeddings of the old samples to match Meghan's previously computed embeddings
early.umap$embedding <- oldSamples@reductions[['umap']]@cell.embeddings

#Combine the pca projections and add reduction to new seurat object (all46)
referenceEmbeddings <- oldSamples@reductions[['pca']]@cell.embeddings
all_pca_embeddings <- rbind(referenceEmbeddings, newSamples@reductions[['pca']]@cell.embeddings)
all46@reductions[['pca']] <- CreateDimReducObject(embeddings = all_pca_embeddings, 
                                                  key = "PCA_", assay = DefaultAssay(all46))

#Plot in pca space and save
pdf(paste0(filedir,'pcaCombined', tag, '.pdf'), width=10, height=6)
print(DimPlot(all46, reduction = 'pca', group.by = 'OldNew'))
dev.off()

#Run UMAP on the new dataset, projecting onto the 'early' umap
newSamples.umap = umap_transform(newSamples@reductions[['pca']]@cell.embeddings, early.umap)
#add row names
row.names(newSamples.umap) <- row.names(newSamples@reductions[['pca']]@cell.embeddings)
#add umap reduction to new samples seurat object
newSamples@reductions[['umap']] <- CreateDimReducObject(embeddings = newSamples.umap, 
                                                        key = "UMAP_", assay = DefaultAssay(newSamples))

#Combine the umap embeddings of new and old samples (computed using same fn.)
early.umap_embeddings <- data.frame(early.umap[['embedding']])
row.names(early.umap_embeddings) <- row.names(oldSamples@reductions[['pca']]@cell.embeddings)
newSample_embeddings <- data.frame(newSamples.umap)
colnames(early.umap_embeddings) <- c('X1', 'X2')
colnames(newSample_embeddings) <- c('X1', 'X2')
all_umap_embeddings <- data.matrix(rbind(early.umap_embeddings, newSample_embeddings))

#Add umap data to metadata of all46 object
all46@reductions[['umap']] <- CreateDimReducObject(embeddings = all_umap_embeddings, 
                                                   key = 'UMAP_', assay = DefaultAssay(all46))

#Plot umap of all samples and save
pdf(paste0(filedir, 'umapCombined', tag, '.pdf'), width=10, height=6)
print(DimPlot(all46, reduction = 'umap', group.by = "OldNew"))
dev.off()

#Now match nearest neighbors in pca space (first 50 dims) and assign clustering accordingly
new_to_old_index <- queryKNN(referenceEmbeddings[,1:50], newSamples@reductions[['pca']]@cell.embeddings[,1:50],
                             k=30, BNPARAM=KmknnParam())

## Define function to find mode of a list (mode of 30 neighbor clusters is used 
#to assign clusters)
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#Create function to copy the cluster assignments based on nearest neighbor assignments
copy2Query <- function(refClusters, queryIndex)
{
  refClusters <- oldClusters
  queryIndex <- new_to_old_index[['index']]
  lengthQuery <- dim(queryIndex)[[1]]
  output <- matrix(data=NA,nrow=lengthQuery,ncol=1)
  for (i in 1:lengthQuery){
    output[i,1] <- getmode(refClusters[queryIndex[i,]])
  }
  return(output)
}

#Run fn to assign new clusters and add all clusters to metadata
#Convert from factor to regular list
oldClusters <- as.numeric(levels(oldSamples@meta.data[['RNA_snn_res.0.05']]))[oldSamples@meta.data[['RNA_snn_res.0.05']]]
#Assign clusters
newClusters <- copy2Query(oldClusters, new_to_old_index[['index']])
#Convert both cluster assignments to dataframe and combine
newClustersdf <- data.frame('X' = newClusters)
oldClustersdf <- data.frame('X' = oldClusters)
clusterInfo <- rbind(oldClustersdf, newClustersdf)
#Convert back to factor and add to metadata
clusterInfo$X <- as.factor(clusterInfo$X)
all46[['clusters_noHarmony_res.0.05_30Neighbors']] <- clusterInfo$X

#Plot umap with clusters and save
pdf(paste0(filedir,'umapCombinedClustered', tag, '.pdf'), width=10, height=6)
print(DimPlot(all46, reduction = 'umap', group.by = 'clusters_noHarmony_res.0.05_30Neighbors'))
dev.off()

##Write out to csv for diversity calculations
## Notes, previously I did diversity calculations in R and generated plots, however, in the CMML dataset,
# as I got more samples I noticed that R wasn't plotting all of the diversity curves I was trying to visualize, 
# so I decided to export the data and calculated diversity in Julia
# create temporary files for storing names and numbers
file1 <- "/Users/4472241/scCode/layerUMAP+Clustering/outputDataNames_20201111.csv"
file2 <- "/Users/4472241/scCode/layerUMAP+Clustering/outputData_20201111.csv"
date <- paste0('_',Sys.Date(),'_')
# save the cluster identity for each individual cell
# note: you need to change the $RNA_snn_res.0.6 part if you identify a better resolution for clustering your data 
cat(all46@meta.data$clusters_noHarmony_res.0.05_30Neighbors, file=file2, sep=",\n")
# save the orig.ident per cell 
listNames <- all46@meta.data$orig.ident
write.table(data.frame(listNames),
            row.names=FALSE,
            col.names = FALSE, 
            file = file1,
            sep=",")
# merging the two files together in grouping to ultimately write out the number of cells per identity per cluster (so in this case, we'd end up with A being a table that has 5 samples x number of cluster idenfied number of rows and two columns )
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
divout <- paste(filedir,"CellBreakdown_PerClusterPerType_res-0.05_layeringOnNewSamples",tag,date,".csv",sep="") 
write.csv(A, file=divout)
write.csv(all46@reductions[['umap']]@cell.embeddings, paste0(filedir, "umap_EmbeddingsAll_layeredNewSamples_noHarmony", date, '.csv'))

#Save rds file of all the samples (optional)
saveRDS(all46, paste0(filedir, 'all_46_samples_umap_layered', tag, '.rds'))

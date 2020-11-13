library(Seurat)
library(matrixStats)
## Input two seurat objects. 
## refObject must be scaled and have a pca reduction.
## newObject should be normalized but does not need to be scaled.
## This function scales newObject, using the scaling in refObject, and
#then computes pca using the basis vector of the pca reduction used for refObject.
## scale_max sets the maximum value for the scaled newObject matrix. The default
#used in Seurat function ScaleData is 10
ScaleAndProjectPCA <- function(newObject, refObject, scale_max = 10){
  
  # inspired by https://www.r-bloggers.com/a-faster-scale-function/
  refMat <- data.matrix(refObject@assays[['RNA']]@data[VariableFeatures(refObject),])
  mat <- data.matrix(newObject@assays[['RNA']]@data[VariableFeatures(refObject),])
  
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

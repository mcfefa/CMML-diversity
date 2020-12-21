library(Seurat)
library(matrixStats)
library(future.apply)
## Input two seurat objects. 
## refObject must be scaled and have a pca reduction.
## newObject should be normalized but does not need to be scaled.
## This function scales newObject, regressing out the variables regressed out in the original
## scaling of the refObject ("nCount_RNA" and "percent.mito" in our case) and
#then computes pca using the basis vector of the pca reduction used for refObject.
## scale_max sets the maximum value for the scaled newObject matrix. The default
#used in Seurat function "ScaleData" is 10
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
  
  #Define chunk points function (from Seurat)
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
  
  #From Seurat originally but modified to regress out new data based on the old Data (reference)
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
  
  #Run function to regress out variables of new data using the regression of the old (ref) data
  newObject.data <- RegressOutUsingRef(
    data.exprRef = refObject.data[, split.cellsRef[["TRUE"]], drop = FALSE],
    latent.dataRef = latent.dataRef[split.cellsRef[['TRUE']], , drop = FALSE],
    data.exprNew = newObject.data[, split.cellsNew[["TRUE"]], drop = FALSE],
    latent.dataNew = latent.dataNew[split.cellsNew[['TRUE']], , drop = FALSE],
    features.regress = features,
    model.use = "linear",
    use.umi = FALSE
  )
  
  #Also regress out old data variables (as it would be done in Seurat "ScaleData")
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
  
  
  
  #**************Regular scale function post-regression*****************
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
  
  #********Run PCA***********
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

library(Seurat)

#Import the seurat object with all of our samples assigned to clusters
dir <- '/Users/ferrallm/Dropbox (UFL)/papers-in-progress/CMML-scRNAseq-Paper/figures/R data (Seurat objects)/'
allSeurat <- readRDS(paste0(dir, 'allSeurat_39+8_postStandardPipeline_withHarmony_withSeuratGeneScores_05-20-2021.rds'))

# Load pseudobulk matrix
bulkMatrix <- readRDS('./src/03b_pseudoBulk/bulkMatrix.rds')

columns <- t(as.vector(colnames(bulkMatrix)))
clus2genes <- read.csv("./src/03b_pseudoBulk/clus2.up.markers.csv",header=FALSE,sep=,)

## union of top 2000 variable genes & clus 2 up genes
clus2pseudo <- intersect(t(clus2genes),columns)

var.out.bool <- !colnames(bulkMatrix) %in% clus2pseudo
bulkMatrixCtrl <- bulkMatrix[,var.out.bool]

var.in.bool <- colnames(bulkMatrix) %in% clus2pseudo
bulkMatrixScore <- bulkMatrix[,var.in.bool]

## naive score
scores = rowMeans(bulkMatrixScore) - rowMeans(bulkMatrixCtrl)
scores

## naive score - random sample of 380 "control" highly variable genes
scores2 = rowMeans(bulkMatrixScore) - rowMeans(bulkMatrixCtrl[,sample(ncol(bulkMatrixCtrl),5)])
scores2

## export data to csv
df <- data.frame(scores,scores2)
write.csv(df, "./src/03b_pseudoBulk/scores.csv", row.names=TRUE)

## export gene list/signature used 
write.csv(clus2pseudo, "./src/03b_pseudoBulk/scores_genelist.csv", row.names=TRUE)

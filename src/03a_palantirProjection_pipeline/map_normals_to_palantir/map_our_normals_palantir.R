# Map our normals to setty rep 1
library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
library(patchwork)
library(BiocNeighbors)

#WORKFLOW: INTEGRATE THE THREE REPLICATES FROM THE SETTY PAPER, THEN MAP REPS 1 AND 2
#TO THE EMBEDDINGS FROM REP3. PROCEED AS NORMAL INTEGRATING EACH CMML SAMPLE AND 
#FINDING THE NEAREST NEIGHBOR TO FIND EMBEDDINGS

### Use Palantir rep3 to match embeddings from all 39 samples (normal +cmml + hma)
###Plot here and then export embeddings and metadata to csv for analysis in python

rm(list=ls())

reps123 <- readRDS(paste0('/Users/Brian/scCode/mapToRep1/reps123Integrated.rds'))
branchProbs_reps123 <- readRDS(paste0('/Users/Brian/scCode/mapToRep1/branchProbsAll3.rds'))
diffPot_reps123 <- readRDS(paste0('/Users/Brian/scCode/mapToRep1/diffPotential.rds'))
branchProbs_reps123$diffPot <- diffPot_reps123$`Diff. Potential`
palantirEmbeddings <- data.frame("X" = reps123@meta.data[["Palantir Embeddings"]][["X"]], 
                                 "Y" = reps123@meta.data[["Palantir Embeddings"]][["Y"]])

allSeurat <- readRDS('/Users/Brian/scCode/allSeurat_withsingleR_withBias_03092022.rds')
normals <- subset(allSeurat, subset = tech != "MCC")

# Re-run PCA with the slightly reduced common gene set (463 down to 452 genes)
reps123 <- NormalizeData(reps123)
reps123 <- ScaleData(reps123)
reps123 <- RunPCA(reps123, dims = 1:30)

source('/Users/Brian/scCode/densityPlot_redblue.R')

normal.list <- SplitObject(normals, split.by = "orig.ident")

for (object in normal.list) {
  # Normalize and scale object
  object <- NormalizeData(object)
  object <- ScaleData(object)
  
  # Find anchors to project their metadata to ours
  anchors <- FindTransferAnchors(reference = reps123, 
                                 query = object,
                                 dims = 1:30, 
                                 reference.reduction = "pca")
  
  # Project the pseudotime/diffPot variables from Setty to nat imm
  branch_Probs_mat <- t(as.matrix(branchProbs_reps123))
  colnames(branch_Probs_mat) <- reps123@assays$RNA@counts@Dimnames[[2]]
  predictions_branchProbs <- TransferData(anchorset = anchors, refdata = branch_Probs_mat,
                                          dims = 1:30, k.weight = 30)
  pseudotime_add_metadata <- as.data.frame(t(as.matrix(predictions_branchProbs@data)))
  object <- AddMetaData(object, metadata = pseudotime_add_metadata)
  
  # Project the palantir Embeddings from Setty to nat imm
  palEmbed_mat <- t(as.matrix(palantirEmbeddings))
  colnames(palEmbed_mat) <- reps123@assays$RNA@counts@Dimnames[[2]]
  predictions_palEmbeddings <- TransferData(anchorset = anchors, refdata = palEmbed_mat,
                                            dims = 1:30, k.weight = 2)
  palEmbeddings_add_metadata <- as.data.frame(t(as.matrix(predictions_palEmbeddings@data)))
  palEmbeddings_add_metadata <- palEmbeddings_add_metadata[which(!is.na(palEmbeddings_add_metadata[,1])),]
  object <- AddMetaData(object, metadata = palEmbeddings_add_metadata)
  
  palEmbeddings <- palEmbeddings_add_metadata[which(!is.na(palEmbeddings_add_metadata[,1])),]
  densityPlot_redblue(palEmbeddings, palantirEmbeddings, object@meta.data[["orig.ident"]][1], "/Users/Brian/scCode/palantirMapToRep1_normals/")
}

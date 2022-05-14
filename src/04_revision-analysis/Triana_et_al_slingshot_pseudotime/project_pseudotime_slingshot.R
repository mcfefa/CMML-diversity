# Project pseudotime from healthy reference (Healthy.rds from Nat. Imm. paper)
# onto our cells. Demonstrate that myelocyte/cDC lineage is overrepresented in 
# mono bias patients compared to other
library(Seurat)
library(ggplot2)

rm(list = ls())

# Import our integrated data
allSeurat <- readRDS('/Users/Brian/scCode/singleR/seuratObject_allSamples_withsingleRType.rds')

# Import the normal 
reference_healthy <- readRDS('/Users/Brian/scCode/nat_imm_sc_ref_maps/Healthy.rds')

# Subset to CD34+
reference_healthy_cd34_plus <- subset(x = reference_healthy, subset = CD34plus == TRUE)
# Rename to reference, this is our reference for mapping onto
reference <- reference_healthy_cd34_plus

# Subset to the same genes
genes_reference <- reference@assays$RNA@counts@Dimnames[[1]]
genes_allSeurat <- allSeurat@assays$RNA@counts@Dimnames[[1]]

# Fix some genes that weren't common due to "." being "-" in allSeurat
genes_not_common <- genes_reference[which(!genes_reference %in% genes_allSeurat)]
genes_fixed <- gsub("[.]", "-", genes_not_common)
reference@assays$RNA@counts@Dimnames[[1]][which(!genes_reference %in% genes_allSeurat)] <- genes_fixed
genes_reference <- reference@assays$RNA@counts@Dimnames[[1]]
genes_common <- genes_reference[which(genes_reference %in% genes_allSeurat)]
write.csv(genes_common, '/Users/Brian/scCode/nat_imm_sc_ref_maps/common_cell_type_genes.csv', row.names=FALSE)

table(genes_common %in% genes_allSeurat) # 452 common genes

# Subset to  common genes
#allSeurat
subset.matrix <- allSeurat@assays$RNA@counts[genes_common, ] # Pull the raw expression matrix from the original Seurat object containing only the genes of interest
allSeurat_subset <- CreateSeuratObject(subset.matrix) # Create a new Seurat object with just the genes of interest
metadata_ours <- allSeurat@meta.data[, c("orig.ident", "tech", "percent.mito", "clusterResolution_0.05", "consensusSingleR")] # Pull the identities from the original Seurat object as a data.frame
allSeurat_subset <- AddMetaData(object = allSeurat_subset, metadata = metadata_ours) # Add the idents to the meta.data slot

#reference
subset.matrix <- reference@assays$RNA@counts[genes_common, ] # Pull the raw expression matrix from the original Seurat object containing only the genes of interest
ref_subset <- CreateSeuratObject(subset.matrix) # Create a new Seurat object with just the genes of interest
metadata_theirs <- reference@meta.data[, c("orig.ident", "nCount_RNA" , "nFeature_RNA", 
                                           "ct", "RNA_snn_res.0.9")]
pseudotime_mat <- reference@meta.data[, c("Erythroid", "Megakaryocte" , "Myelocytes", 
                                          "B.cells", "cDC")]
pseudotime_mat[is.na(pseudotime_mat)] <- 0

# Add mofa umap reduction as two distinct metadata
mofa_metadata1 <- reference@reductions$MOFAUMAP@cell.embeddings[,"MOFAUMAP_1"]
mofa_metadata2 <- reference@reductions$MOFAUMAP@cell.embeddings[,"MOFAUMAP_2"]
ref_subset <- AddMetaData(object = ref_subset, metadata = metadata_theirs)
ref_subset <- AddMetaData(object = ref_subset, metadata = pseudotime_mat)
ref_subset <- AddMetaData(object = ref_subset, metadata = mofa_metadata1, col.name = "MOFA_1")
ref_subset <- AddMetaData(object = ref_subset, metadata = mofa_metadata2, col.name = "MOFA_2")

# Re-run PCA with the slightly reduced common gene set (463 down to 452 genes)
ref_subset <- NormalizeData(ref_subset)
ref_subset <- ScaleData(ref_subset)
ref_subset <- RunPCA(ref_subset, features = genes_common, dims = 1:30)

allSeurat_subset <- NormalizeData(allSeurat_subset)
allSeurat_subset <- ScaleData(allSeurat_subset)

# Save the reference seurat object with a PCA projection from the common genes with our seurat object
#saveRDS(ref_subset, '/Users/Brian/scCode/nat_imm_sc_ref_maps/nat_imm_reference_seurat_obj_03092022.rds')
ref_subset <- readRDS('/Users/Brian/scCode/nat_imm_sc_ref_maps/nat_imm_reference_seurat_obj_03092022.rds')

# Save our scaled seurat object ready to be used for finding anchors
#saveRDS(allSeurat_subset, '/Users/Brian/scCode/nat_imm_sc_ref_maps/allSeurat_subset_postScale_03092022.rds')
allSeurat_subset <- readRDS('/Users/Brian/scCode/nat_imm_sc_ref_maps/allSeurat_subset_postScale_03092022.rds')
genes_common <- read.csv('/Users/Brian/scCode/nat_imm_sc_ref_maps/common_cell_type_genes.csv')$x

pseudotime_mat <- ref_subset@meta.data[, c("Erythroid", "Megakaryocte" , "Myelocytes", 
                                          "B.cells", "cDC")]

# Find anchors to project their metadata to ours
anchors <- FindTransferAnchors(reference = ref_subset, 
                                         query = allSeurat_subset,
                                         dims = 1:30, 
                                         reference.reduction = "pca", features = genes_common)

# Project the pseudotime variables from theirs to ours (might need to be re-written)
predictions_pseudotime <- TransferData(anchorset = anchors, refdata = t(as.matrix(pseudotime_mat)),
                            dims = 1:30)

# Add the projected metadata to our seurat object
pseudotime_add_metadata <- as.data.frame(t(as.matrix(predictions_pseudotime@data)))
allSeurat_subset <- AddMetaData(allSeurat_subset, metadata = pseudotime_add_metadata)

VlnPlot(allSeurat_subset, features = "Myelocytes", group.by = "orig.ident", pt.size = 0) + NoLegend()

# Map metadata from orig.ident to Bias group
normals <- unique(allSeurat_subset$orig.ident)[1:8]
mono_bias <- c("6AD001", "SF14031800065", "SF14072200012", "SF14111400033", "2V001", 
               "SF100109101914", "SF100109106293", "SF12062800475", "SF16112300029")
mep_bias <- c("LTB3966", "SF14050700419", "SF15010200008", "SF14061300036", "SF100109111451",
              "SF12042500035", "SF14101000049", "SF16026800045", "SF16072200003")
allSeurat_subset$bias <- "Normal-like"
allSeurat_subset$bias[allSeurat_subset$orig.ident %in% mono_bias] <- "Mono bias"
allSeurat_subset$bias[allSeurat_subset$orig.ident %in% mep_bias] <- "MEP bias"
allSeurat_subset$bias[allSeurat_subset$orig.ident %in% normals] <- "Normal"

color_scale <- c("#65CB01", "#32CEFF", "black", "#E70F24")

pdf('/Users/Brian/scCode/project_pseudotime_nat_imm/Monocyte_colors_fixed.pdf')
VlnPlot(allSeurat_subset, features = "Myelocytes", group.by = "bias", pt.size = 0, 
        cols = color_scale[c(2, 4, 1, 3)], sort = T) + NoLegend() +
  ggtitle("Monocyte Pseudotime") + xlab("Bias") + theme(axis.title.x=element_blank())
dev.off()

pdf('/Users/Brian/scCode/project_pseudotime_nat_imm/Erythroid_colors_fixed..pdf')
VlnPlot(allSeurat_subset, features = "Erythroid", group.by = "bias", pt.size = 0,
        cols = color_scale[c(1, 4, 3, 2)], sort = T) + NoLegend() +
  ggtitle("Erythroid Pseudotime") + xlab("Bias") + theme(axis.title.x=element_blank())
dev.off()

pdf('/Users/Brian/scCode/project_pseudotime_nat_imm/cDC_colors_fixed..pdf')
VlnPlot(allSeurat_subset, features = "cDC", group.by = "bias", pt.size = 0,
        cols = color_scale[c(4, 1, 3, 2)], sort = T) + NoLegend() +
  ggtitle("cDC Pseudotime") + xlab("Bias") + theme(axis.title.x=element_blank())
dev.off()

pdf('/Users/Brian/scCode/project_pseudotime_nat_imm/B_cell_colors_fixed..pdf')
VlnPlot(allSeurat_subset, features = "B.cells", group.by = "bias", pt.size = 0,
        cols = color_scale[c(3, 4, 1, 2)], sort = T) + NoLegend() +
  ggtitle("B cell Pseudotime") + xlab("Bias") + theme(axis.title.x=element_blank())
dev.off()

pdf('/Users/Brian/scCode/project_pseudotime_nat_imm/Megakaryocyte_colors_fixed..pdf')
VlnPlot(allSeurat_subset, features = "Megakaryocte", group.by = "bias", pt.size = 0,
        cols = color_scale[c(4, 3, 1, 2)], sort = T) + NoLegend() +
  ggtitle("Megakaryocyte Pseudotime") + xlab("Bias") + theme(axis.title.x=element_blank())
dev.off()

  

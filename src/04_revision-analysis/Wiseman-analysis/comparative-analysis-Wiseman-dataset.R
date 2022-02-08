library(Seurat)
library(Matrix)
library(stringr)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(RColorBrewer)
library(gplots)

# sessionInfo()
# R version 4.1.2 (2021-11-01)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Big Sur 11.6.2
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base
# 
# other attached packages:
#   [1] gplots_3.1.1       RColorBrewer_1.1-2 forcats_0.5.1      dplyr_1.0.7        purrr_0.3.4
# [6] readr_2.1.2        tidyr_1.1.4        tibble_3.1.5       tidyverse_1.3.1    patchwork_1.1.1
# [11] ggplot2_3.3.5      stringr_1.4.0      Matrix_1.3-4       SeuratObject_4.0.2 Seurat_4.0.5
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15            colorspace_2.0-2      deldir_1.0-5          ellipsis_0.3.2
# [5] ggridges_0.5.3        fs_1.5.0              rstudioapi_0.13       spatstat.data_2.1-0
# [9] leiden_0.3.9          listenv_0.8.0         ggrepel_0.9.1         fansi_0.5.0
# [13] lubridate_1.8.0       xml2_1.3.3            codetools_0.2-18      splines_4.1.2
# [17] polyclip_1.10-0       jsonlite_1.7.2        broom_0.7.12          ica_1.0-2
# [21] cluster_2.1.2         dbplyr_2.1.1          png_0.1-7             uwot_0.1.10
# [25] shiny_1.7.1           sctransform_0.3.2     spatstat.sparse_2.0-0 compiler_4.1.2
# [29] httr_1.4.2            backports_1.4.1       assertthat_0.2.1      fastmap_1.1.0
# [33] lazyeval_0.2.2        cli_3.0.1             later_1.3.0           htmltools_0.5.2
# [37] tools_4.1.2           igraph_1.2.6          gtable_0.3.0          glue_1.4.2
# [41] RANN_2.6.1            reshape2_1.4.4        Rcpp_1.0.7            scattermore_0.7
# [45] cellranger_1.1.0      vctrs_0.3.8           nlme_3.1-153          lmtest_0.9-38
# [49] globals_0.14.0        rvest_1.0.2           mime_0.12             miniUI_0.1.1.1
# [53] lifecycle_1.0.1       irlba_2.3.3           gtools_3.9.2          goftest_1.2-3
# [57] future_1.22.1         MASS_7.3-54           zoo_1.8-9             scales_1.1.1
# [61] spatstat.core_2.3-0   hms_1.1.1             promises_1.2.0.1      spatstat.utils_2.2-0
# [65] parallel_4.1.2        reticulate_1.22       pbapply_1.5-0         gridExtra_2.3
# [69] rpart_4.1-15          stringi_1.7.5         caTools_1.18.2        bitops_1.0-7
# [73] rlang_0.4.11          pkgconfig_2.0.3       matrixStats_0.61.0    lattice_0.20-45
# [77] ROCR_1.0-11           tensor_1.5            htmlwidgets_1.5.4     cowplot_1.1.1
# [81] tidyselect_1.1.1      parallelly_1.28.1     RcppAnnoy_0.0.19      plyr_1.8.6
# [85] magrittr_2.0.1        R6_2.5.1              generics_0.1.0        DBI_1.1.2
# [89] pillar_1.6.3          haven_2.4.3           withr_2.4.2           mgcv_1.8-38
# [93] fitdistrplus_1.1-6    survival_3.2-13       abind_1.4-5           future.apply_1.8.1
# [97] modelr_0.1.8          crayon_1.4.1          KernSmooth_2.23-20    utf8_1.2.2
# [101] spatstat.geom_2.3-0   plotly_4.10.0         tzdb_0.2.0            grid_4.1.2
# [105] readxl_1.3.1          data.table_1.14.2     reprex_2.0.1          digest_0.6.28
# [109] xtable_1.8-4          httpuv_1.6.3          munsell_0.5.0         viridisLite_0.4.0

rm(list = ls())

setwd("~/GitHub/CMML-diversity")
datadir <-'~/GitHub/CMML-diversity/data/Wiseman/'
dir <- '~/GitHub/CMML-diversity/src/04_revision-analysis/Wiseman-analysis/'
date = "_02-04-2022"

############## IMPORT WISEMAN SAMPLES ##################################

W1.data <- Read10X(data.dir = paste(datadir,"W1_BC278/",sep=""))
W1 <- CreateSeuratObject(counts = W1.data, project = "BC278")
W1 <- RenameCells(object=W1, add.cell.id="BC278")

W2.data <- Read10X(data.dir = paste(datadir,"W2_BC416/",sep=""))
W2 <- CreateSeuratObject(counts = W2.data, project = "BC416")
W2 <- RenameCells(object=W2, add.cell.id="BC416")

W3.data <- Read10X(data.dir = paste(datadir,"W3_BC543/",sep=""))
W3 <- CreateSeuratObject(counts = W3.data, project = "BC543")
W3 <- RenameCells(object=W3, add.cell.id="BC543")

W4a.data <- Read10X(data.dir = paste(datadir,"W4a_BC572/",sep=""))
W4a <- CreateSeuratObject(counts = W4a.data, project = "BC572")
W4a <- RenameCells(object=W4a, add.cell.id="BC572")

W4b.data <- Read10X(data.dir = paste(datadir,"W4b_BC572/",sep=""))
W4b <- CreateSeuratObject(counts = W4b.data, project = "BC572")
W4b <- RenameCells(object=W4b, add.cell.id="BC572")

W5.data <- Read10X(data.dir = paste(datadir,"W5_BC746/",sep=""))
W5 <- CreateSeuratObject(counts = W5.data, project = "BC746")
W5 <- RenameCells(object=W5, add.cell.id="BC746")

W6.data <- Read10X(data.dir = paste(datadir,"W6_BC776/",sep=""))
W6 <- CreateSeuratObject(counts = W6.data, project = "BC776")
W6 <- RenameCells(object=W6, add.cell.id="BC776")

W7.data <- Read10X(data.dir = paste(datadir,"W7_BC786/",sep=""))
W7 <- CreateSeuratObject(counts = W7.data, project = "BC786")
W7 <- RenameCells(object=W7, add.cell.id="BC786")

W8.data <- Read10X(data.dir = paste(datadir,"W8_HV1/",sep=""))
W8 <- CreateSeuratObject(counts = W8.data, project = "HV1")
W8 <- RenameCells(object=W8, add.cell.id="HV1")

W9.data <- Read10X(data.dir = paste(datadir,"W9_HV2/",sep=""))
W9 <- CreateSeuratObject(counts = W9.data, project = "HV2")
W9 <- RenameCells(object=W9, add.cell.id="HV2")

W10.data <- Read10X(data.dir = paste(datadir,"W10_HV3/",sep=""))
W10 <- CreateSeuratObject(counts = W10.data, project = "HV3")
W10 <- RenameCells(object=W10, add.cell.id="HV3")

WisemanCMML <- merge(x=W1, y=c(W2,W3,W4a,W4b,W5,W6,W7), add.cell.ids=c("BC278","BC416","BC543","BC572","BC572","BC746","BC776","BC786"), project="WisemanCMML")
saveRDS(WisemanCMML,paste(datadir,"Wiseman-CMML-Cohort_Raw-No-Cutoffs-Applied",date,".rds",sep=""))

WisemanNormal <- merge(x=W8, y=c(W9,W10), add.cell.ids=c("HV1","HV2","HV3"), project="WisemanNormal")
saveRDS(WisemanNormal,paste(datadir,"Wiseman-Normal-Cohort_Raw-No-Cutoffs-Applied",date,".rds",sep=""))

# Check Dimesions
levels(WisemanCMML)
# [1] "BC278" "BC416" "BC543" "BC572" "BC746" "BC776" "BC786"
dim(WisemanCMML)
# [1] 36601  6468

levels(WisemanNormal)
# [1] "HV1" "HV2" "HV3"
dim(WisemanNormal)
# [1] 36601  1966

#Add to metadata to separate normals from cmml samples
WisemanCMML[['orig.ident']] <- WisemanCMML@active.ident
WisemanCMML[['tech']] <- 'Wiseman'

WisemanNormal[['orig.ident']] <- WisemanNormal@active.ident
WisemanNormal[['tech']] <- 'Wiseman'

#Get summary of nFeatures
summary(WisemanCMML@meta.data$nFeature_RNA)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 40    1105    1986    1974    2602    6071 

summary(WisemanNormal@meta.data$nFeature_RNA)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  29    1967    2270    2271    2599    5071 

#Calculate mito percentage
mito.genes <- grep(pattern = "^MT-", WisemanCMML@assays$RNA@counts@Dimnames[[1]], value = TRUE)
percent.mito.cmml <- Matrix::colSums(x=GetAssayData(object=WisemanCMML, slot='counts')[mito.genes, ])/Matrix::colSums(x=GetAssayData(object=WisemanCMML, slot='counts'))
WisemanCMML[['percent.mito']] <- percent.mito.cmml

percent.mito.normal <- Matrix::colSums(x=GetAssayData(object=WisemanNormal, slot='counts')[mito.genes, ])/Matrix::colSums(x=GetAssayData(object=WisemanNormal, slot='counts'))
WisemanNormal[['percent.mito']] <- percent.mito.normal


#Set standard cutoffs across all CMML samples (same as nFeatLowerN_norm, diff mito, diff upperFeatN)
nFeatUpper_CMML <- mean(WisemanCMML@meta.data$nFeature_RNA, na.rm=TRUE) + 2*sd(WisemanCMML@meta.data$nFeature_RNA, na.rm=TRUE) # 4072
nFeatLower_CMML <- 450 
perMitoUpper_CMML <- 0.25 

nFeatUpper_norm <- mean(WisemanNormal@meta.data$nFeature_RNA, na.rm=TRUE) + 2*sd(WisemanNormal@meta.data$nFeature_RNA, na.rm=TRUE) # 3645
nFeatLower_norm <- 450 
perMitoUpper_norm <- 0.05 

# Plot the violin and scatter plot of number of features with the cutoffs shown
# CMML number of features
pdf(paste(dir, "Wiseman_CMML_All_nFeatRNA", date, ".pdf",sep=""), width = 11, height = 3)
vPlot <- VlnPlot(object = WisemanCMML, features = c("nFeature_RNA"), pt.size = 0, combine = T, group.by = "orig.ident") + 
  scale_y_continuous(limits = c(0,9000)) & 
  geom_hline(yintercept = nFeatUpper_CMML) & geom_hline(yintercept = nFeatLower_CMML)
print(vPlot)
dev.off()

pdf(paste(dir, "Wiseman_Normal_All_nFeatRNA", date, ".pdf",sep=""), width = 11, height = 3)
vPlot <- VlnPlot(object = WisemanNormal, features = c("nFeature_RNA"), pt.size = 0, combine = T, group.by = "orig.ident") + 
  scale_y_continuous(limits = c(0,9000)) & 
  geom_hline(yintercept = nFeatUpper_norm) & geom_hline(yintercept = nFeatLower_norm)
print(vPlot)
dev.off()

# Plot filters based on mitochondrial content
# CMML mito plot
pdf(paste(dir, "Wiseman_CMML_mito", date,".pdf",sep=""), width = 10, height = 3)
vPlot <- VlnPlot(object = WisemanCMML, features = c("percent.mito"), pt.size = 0, group.by = "orig.ident") +
  scale_y_continuous(limits = c(0,1)) & 
  geom_hline(yintercept = perMitoUpper_CMML)
print(vPlot)
dev.off()

pdf(paste(dir, "Wiseman_Normal_mito", date,".pdf",sep=""), width = 10, height = 3)
vPlot <- VlnPlot(object = WisemanNormal, features = c("percent.mito"), pt.size = 0, group.by = "orig.ident") +
  scale_y_continuous(limits = c(0,1)) & 
  geom_hline(yintercept = perMitoUpper_norm)
print(vPlot)
dev.off()

# subset data based on cutoffs
WisemanCMML.qc <- subset(x=WisemanCMML, subset=nFeature_RNA > nFeatLower_CMML & nFeature_RNA < nFeatUpper_CMML & percent.mito < perMitoUpper_CMML)

dim(WisemanCMML) # [1] 36601  6468
dim(WisemanCMML.qc) # [1] 36601  5813
# # Keep about 90% of cells, 10% fall outside bounds

WisemanNormal.qc <- subset(x=WisemanNormal, subset=nFeature_RNA > nFeatLower_norm & nFeature_RNA < nFeatUpper_norm & percent.mito < perMitoUpper_norm)

dim(WisemanNormal) # [1] 36601  1966
dim(WisemanNormal.qc) # [1] 36601  1768
# # Keep about 90% of cells, 10% fall outside bounds

saveRDS(WisemanCMML.qc, paste0(dir, "Wiseman_CMML_postQC",date,".rds",sep=""))

##### Merge both post-qc and combine to give us our seurat object for future use
Wiseman <- merge(WisemanNormal.qc, WisemanCMML.qc)
saveRDS(Wiseman, paste0(dir, "Wiseman_postQC_CMML7+healthy3",date,".rds",sep=""))

##### Post-QC Processing #####################################
# Log normalize 
Wiseman <- NormalizeData(Wiseman, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features 
Wiseman <- FindVariableFeatures(Wiseman, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes and plot
top10 <- head(VariableFeatures(Wiseman), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Wiseman)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf(paste(dir, "Wiseman_HVG_labeled-top-10", date,".pdf",sep=""))
vPlot <- plot1 + plot2
print(vPlot)
dev.off()

# Scale data (only use HVG), regressing out effects of nCountRNA and percent.mito
Wiseman <- ScaleData(Wiseman, features = VariableFeatures(Wiseman), vars.to.regress = c("nCount_RNA","percent.mito"))

# Run, plot and save PCA
Wiseman <- RunPCA(Wiseman, features = VariableFeatures(object = Wiseman))
pdf(paste0(dir, "Wiseman_PCA_allSeurat_standard_seurat_pipeline_preHarmony", date, ".pdf",sep=""))
DimPlot(Wiseman, reduction = "pca", split.by = "tech", pt.size = 0.0001)
dev.off()

# Determine dimensionality of dataset
pdf(paste0(dir, "Wiseman_PCA_Elbow-Plot_allSeurat_standard_seurat_pipeline_preHarmony", date, ".pdf",sep=""))
ElbowPlot(Wiseman, ndims = 50)
dev.off()

saveRDS(Wiseman, paste0(dir, "Wiseman_postQC_CMML7+healthy3_thruPCA",date,".rds",sep=""))

## Skipping Harmony integration & clustering here, focused on validating the results from Ferrall-Fairbanks et al

# Run, plot, and save UMAP projections (based on PCA)
Wiseman <- RunUMAP(Wiseman, dims = 1:50)
pdf(paste0(dir, "Wiseman_UMAP_withoutHarmony",date,".pdf"))
DimPlot(Wiseman, reduction = "umap", group.by = "orig.ident")
dev.off()

## FEATURE PLOT - CD120b
feature_vec <- c("TNFRSF1B")
pdf(paste0(dir, "Wiseman_UMAP_CD120b",date,".pdf",sep=""))
FeaturePlot(Wiseman, features = feature_vec)
dev.off()

## check FcgR expression per Reviewer comment
options(future.globals.maxSize = 128000*1024*1024)
feature_vec <- c("FCGR2A")
pdf(paste0(dir,"Wiseman_UMAP_FcgR",date,".pdf",sep=""))
FeaturePlot(Wiseman, features = feature_vec)
dev.off()


#<--------------------- here 


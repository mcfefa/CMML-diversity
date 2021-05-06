library(Seurat)
library(Matrix)
library(stringr)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(RColorBrewer)
library(gplots)

#Redo this sessionInfo() call with minimum packages after restarting R kernel
#sessionInfo() #Output:
#R version 4.0.2 (2020-06-22)
#Platform: x86_64-apple-darwin17.0 (64-bit)
#Running under: macOS Catalina 10.15.7

#Matrix products: default
#BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
#LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

#locale:
#  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

#attached base packages:
#  [1] splines   stats4    parallel  stats     graphics  grDevices utils     datasets 
#[9] methods   base     

#other attached packages:
#  [1] Matrix_1.3-2                data.table_1.14.0           copykat_1.0.4              
#[4] maftools_2.4.12             raster_3.4-5                sp_1.4-5                   
#[7] BiocNeighbors_1.6.0         harmony_1.0                 Rcpp_1.0.6                 
#[10] gplots_3.1.1                RColorBrewer_1.1-2          forcats_0.5.1              
#[13] stringr_1.4.0               dplyr_1.0.5                 purrr_0.3.4                
#[16] readr_1.4.0                 tidyr_1.1.3                 tibble_3.1.0               
#[19] tidyverse_1.3.0             patchwork_1.1.1             ggplot2_3.3.3              
#[22] SeuratObject_4.0.0          Seurat_4.0.1                pheatmap_1.0.12            
#[25] devtools_2.3.2              usethis_2.0.1               deepSNV_1.34.1             
#[28] VariantAnnotation_1.34.0    Rsamtools_2.4.0             VGAM_1.1-5                 
#[31] Biostrings_2.56.0           XVector_0.28.0              SummarizedExperiment_1.18.2
#[34] DelayedArray_0.14.1         matrixStats_0.58.0          Biobase_2.48.0             
#[37] GenomicRanges_1.40.0        GenomeInfoDb_1.24.2         IRanges_2.22.2             
#[40] S4Vectors_0.26.1            BiocGenerics_0.34.0         mitoClone_0.1.0            

#loaded via a namespace (and not attached):
#  [1] utf8_1.2.1               reticulate_1.18          tidyselect_1.1.0        
#[4] RSQLite_2.2.5            AnnotationDbi_1.50.3     htmlwidgets_1.5.3       
#[7] grid_4.0.2               BiocParallel_1.22.0      Rtsne_0.15              
#[10] munsell_0.5.0            codetools_0.2-18         ica_1.0-2               
#[13] future_1.21.0            miniUI_0.1.1.1           withr_2.4.1             
#[16] colorspace_2.0-0         rstudioapi_0.13          ROCR_1.0-11             
#[19] tensor_1.5               listenv_0.8.0            labeling_0.4.2          
#[22] GenomeInfoDbData_1.2.3   polyclip_1.10-0          MCMCpack_1.5-0          
#[25] farver_2.1.0             bit64_4.0.5              rprojroot_2.0.2         
#[28] coda_0.19-4              parallelly_1.24.0        vctrs_0.3.7             
#[31] generics_0.1.0           BiocFileCache_1.12.1     R6_2.5.0                
#[34] bitops_1.0-7             spatstat.utils_2.1-0     cachem_1.0.4            
#[37] assertthat_0.2.1         dlm_1.1-5                promises_1.2.0.1        
#[40] scales_1.1.1             gtable_0.3.0             globals_0.14.0          
#[43] conquer_1.0.2            processx_3.5.0           goftest_1.2-2           
#[46] mcmc_0.9-7               MatrixModels_0.5-0       rlang_0.4.10            
#[49] rtracklayer_1.48.0       lazyeval_0.2.2           spatstat.geom_2.0-1     
#[52] broom_0.7.5              reshape2_1.4.4           abind_1.4-5             
#[55] modelr_0.1.8             GenomicFeatures_1.40.1   backports_1.2.1         
#[58] httpuv_1.5.5             tools_4.0.2              ellipsis_0.3.1          
#[61] spatstat.core_2.0-0      sessioninfo_1.1.1        ggridges_0.5.3          
#[64] plyr_1.8.6               progress_1.2.2           zlibbioc_1.34.0         
#[67] RCurl_1.98-1.3           ps_1.6.0                 prettyunits_1.1.1       
#[70] rpart_4.1-15             openssl_1.4.3            deldir_0.2-10           
#[73] pbapply_1.4-3            cowplot_1.1.1            zoo_1.8-9               
#[76] haven_2.3.1              ggrepel_0.9.1            cluster_2.1.1           
#[79] fs_1.5.0                 magrittr_2.0.1           scattermore_0.7         
#[82] SparseM_1.81             reprex_1.0.0             lmtest_0.9-38           
#[85] RANN_2.6.1               parallelDist_0.2.4       fitdistrplus_1.1-3      
#[88] pkgload_1.2.0            hms_1.0.0                mime_0.10               
#[91] xtable_1.8-4             XML_3.99-0.6             readxl_1.3.1            
#[94] gridExtra_2.3            testthat_3.0.2           compiler_4.0.2          
#[97] biomaRt_2.44.4           KernSmooth_2.23-18       crayon_1.4.1            
#[100] htmltools_0.5.1.1        segmented_1.3-4          mgcv_1.8-34             
#[103] later_1.1.0.1            RcppParallel_5.1.4       lubridate_1.7.10        
#[106] DBI_1.1.1                dbplyr_2.1.0             MASS_7.3-53.1           
#[109] rappdirs_0.3.3           cli_2.4.0                igraph_1.2.6            
#[112] pkgconfig_2.0.3          GenomicAlignments_1.24.0 plotly_4.9.3            
#[115] spatstat.sparse_2.0-0    xml2_1.3.2               rvest_1.0.0             
#[118] callr_3.6.0              digest_0.6.27            sctransform_0.3.2       
#[121] RcppAnnoy_0.0.18         spatstat.data_2.1-0      cellranger_1.1.0        
#[124] leiden_0.3.7             uwot_0.1.10              kernlab_0.9-29          
#[127] curl_4.3                 quantreg_5.85            gtools_3.8.2            
#[130] shiny_1.6.0              lifecycle_1.0.0          nlme_3.1-152            
#[133] jsonlite_1.7.2           desc_1.3.0               viridisLite_0.4.0       
#[136] askpass_1.1              BSgenome_1.56.0          fansi_0.4.2             
#[139] pillar_1.5.1             lattice_0.20-41          Rhtslib_1.20.0          
#[142] fastmap_1.1.0            httr_1.4.2               pkgbuild_1.2.0          
#[145] survival_3.2-10          glue_1.4.2               remotes_2.2.0           
#[148] png_0.1-7                bit_4.0.4                mixtools_1.2.0          
#[151] stringi_1.5.3            blob_1.2.1               caTools_1.18.2          
#[154] memoise_2.0.0            irlba_2.3.3              future.apply_1.7.0

rm(list = ls())

dir <- '/Users/4472241/scCode/redoEntirePipeline_05-05-2021/'
date = "_05-06-2021"

###### Switch to Normals ###########################

#Import Hua Pt 3 (previously excluded)
raw.HuaPt3 <- Read10X(data.dir = paste0(dir, 'GRCh38_HuaPt3'))
huaPt3 <- CreateSeuratObject(counts = raw.HuaPt3, project = "HuaPt3", orig.ident = 'HuaPt3')
#Add the tech/dataset it comes from to metadata
huaPt3[['tech']] <- 'Hua'

#Import the rest of the normals
sevenNormal <- readRDS(paste0(dir, 'Normal_CD34only_Zheng+Setty+Hua_UnprocessedCohort_SeuratObj_20200321.rds'))
levels(sevenNormal) #[1] "CD34"     "HuaPt1"   "HuaPt2"   "HuaPt4"   "SettyPt1" "SettyPt2" "SettyPt3"
dim(sevenNormal) #[1] 33694 77751
dim(huaPt3) #[1] 33694  5131

#Split the normals up by sample
sevenNormal[['orig.ident']] <- sevenNormal@active.ident
seven.split <- SplitObject(sevenNormal, split.by = "orig.ident")

#Extract count matrix and re-make seurat object for each sample
new.seven.split <- list()
for (i in 1:length(seven.split)){
  #Get name from list
  name <- names(seven.split)[i]
  #Define which tech/dataset this comes from
  if (name %in% c("SettyPt1", "SettyPt2", "SettyPt3")){
    tech <- "Setty"
  }else if (name %in% 'CD34'){
    tech <- "Zheng"
    #Rename to Zheng
    name <- "Zheng"
  }else{
    tech <- "Hua"
  }
  
  #Get raw counts, assign metadata to separate the samples/datasets
  raw.counts <- seven.split[[i]]@assays[['RNA']]@counts
  seurat <- CreateSeuratObject(counts = raw.counts, project = name, orig.ident = name)
  seurat[['tech']] <- tech
  new.seven.split[[i]] <- seurat
  names(new.seven.split)[i] <- name
}

#Combine all samples from same dataset to respective seurat object
setty <- new.seven.split[["SettyPt1"]]
setty <- merge(setty, new.seven.split[['SettyPt2']])
setty <- merge(setty, new.seven.split[['SettyPt3']])
zheng <- new.seven.split[['Zheng']]
hua <- new.seven.split[["HuaPt1"]]
hua <- merge(hua, new.seven.split[["HuaPt2"]])
hua <- merge(hua, huaPt3)
hua <- merge(hua, new.seven.split[["HuaPt4"]])

summary(setty@meta.data$nFeature_RNA)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 12     605    1622    1579    2378    7412 
summary(zheng@meta.data$nFeature_RNA)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 394     959    1293    1300    1571    3812 
summary(hua@meta.data$nFeature_RNA)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 404     994    1377    1491    1810    5195

#Set standard cutoffs across all normals
nFeatLowerN_norm <- 450
perMitoUpperN_norm <- 0.05

#Set cutoffs for Setty (upper # features cutoff will be different for each tech/dataset)
mito.genes <- grep(pattern = "^MT-", setty@assays$RNA@counts@Dimnames[[1]], 
                   value = TRUE)
percent.mito <- Matrix::colSums(x=GetAssayData(object=setty, slot='counts')[mito.genes, ])/Matrix::colSums(x=GetAssayData(object=setty, slot='counts'))
setty[['percent.mito']] <- percent.mito

nFeatUpperN_setty <- mean(setty@meta.data$nFeature_RNA, na.rm=TRUE) + 2*sd(setty@meta.data$nFeature_RNA, na.rm=TRUE) # [1] 3755.362

#Set cutoffs for zheng (upper # features cutoff will be different for each dataset (Setty, Hua, Zheng))
mito.genes <- grep(pattern = "^MT-", zheng@assays$RNA@counts@Dimnames[[1]], 
                   value = TRUE)
percent.mito <- Matrix::colSums(x=GetAssayData(object=zheng, slot='counts')[mito.genes, ])/Matrix::colSums(x=GetAssayData(object=zheng, slot='counts'))
zheng[['percent.mito']] <- percent.mito

nFeatUpperN_zheng <- mean(zheng@meta.data$nFeature_RNA, na.rm=TRUE) + 2*sd(zheng@meta.data$nFeature_RNA, na.rm=TRUE) # [1] 2213.11

#Set cutoffs for hua (upper # features cutoff will be different for each dataset (Setty, Hua, Zheng))
mito.genes <- grep(pattern = "^MT-", hua@assays$RNA@counts@Dimnames[[1]], 
                   value = TRUE)
percent.mito <- Matrix::colSums(x=GetAssayData(object=hua, slot='counts')[mito.genes, ])/Matrix::colSums(x=GetAssayData(object=hua, slot='counts'))
hua[['percent.mito']] <- percent.mito

nFeatUpperN_hua <- mean(hua@meta.data$nFeature_RNA, na.rm=TRUE) + 2*sd(hua@meta.data$nFeature_RNA, na.rm=TRUE) # [1] 2773.66

# Plot the violin and scatter plot of number of features with the cutoffs shown
# Setty # features
pdf(paste(dir, "setty_nFeatRNA", date, ".pdf",sep=""))
#Ignore message about scale, plotting range to be same across all samples
vPlot <- VlnPlot(object = setty, features = c("nFeature_RNA"), pt.size = 0, combine = T, group.by = "orig.ident") + 
  scale_y_continuous(limits = c(0,6000)) & 
  geom_hline(yintercept = nFeatUpperN_setty) & geom_hline(yintercept = nFeatLowerN_norm)
print(vPlot)
nGenes <- new.seven.split[["SettyPt1"]]@meta.data$nFeature_RNA
nGenes_sort <- data.frame("nFeatRNA" = sort(nGenes), "X" = 1:length(nGenes))
scatPlot <- ggplot(nGenes_sort, aes(x = X, y = nFeatRNA)) + geom_point() +
  ggtitle("SettyPt1") &
  geom_hline(yintercept = nFeatUpperN_hua) & geom_hline(yintercept = nFeatLowerN_norm)
print(scatPlot)
nGenes <- new.seven.split[["SettyPt2"]]@meta.data$nFeature_RNA
nGenes_sort <- data.frame("nFeatRNA" = sort(nGenes), "X" = 1:length(nGenes))
scatPlot <- ggplot(nGenes_sort, aes(x = X, y = nFeatRNA)) + geom_point() +
  ggtitle("SettyPt2") &
  geom_hline(yintercept = nFeatUpperN_hua) & geom_hline(yintercept = nFeatLowerN_norm)
print(scatPlot)
nGenes <- new.seven.split[["SettyPt3"]]@meta.data$nFeature_RNA
nGenes_sort <- data.frame("nFeatRNA" = sort(nGenes), "X" = 1:length(nGenes))
scatPlot <- ggplot(nGenes_sort, aes(x = X, y = nFeatRNA)) + geom_point() +
  ggtitle("SettyPt3") &
  geom_hline(yintercept = nFeatUpperN_hua) & geom_hline(yintercept = nFeatLowerN_norm)
print(scatPlot)
dev.off()

# Zheng # features
pdf(paste(dir, "zheng_nFeatRNA", date, ".pdf",sep=""))
#Ignore message about scale, plotting range to be same across all samples
vPlot <- VlnPlot(object = zheng, features = c("nFeature_RNA"), pt.size = 0, combine = T)  + 
  scale_y_continuous(limits = c(0,6000)) & 
  geom_hline(yintercept = nFeatUpperN_zheng) & geom_hline(yintercept = nFeatLowerN_norm)
print(vPlot)
nGenes <- new.seven.split[["Zheng"]]@meta.data$nFeature_RNA
nGenes_sort <- data.frame("nFeatRNA" = sort(nGenes), "X" = 1:length(nGenes))
scatPlot <- ggplot(nGenes_sort, aes(x = X, y = nFeatRNA)) + geom_point() +
  ggtitle("Zheng") &
  geom_hline(yintercept = nFeatUpperN_hua) & geom_hline(yintercept = nFeatLowerN_norm)
print(scatPlot)
dev.off()

# Hua # features
pdf(paste(dir, "hua_nFeatRNA", date, ".pdf",sep=""))
#Ignore message about scale, plotting range to be same across all samples
vPlot <- VlnPlot(object = hua, features = c("nFeature_RNA"), pt.size = 0, combine = T, group.by = "orig.ident") + 
  scale_y_continuous(limits = c(0,6000)) & 
  geom_hline(yintercept = nFeatUpperN_hua) & geom_hline(yintercept = nFeatLowerN_norm)
print(vPlot)
nGenes <- new.seven.split[["HuaPt1"]]@meta.data$nFeature_RNA
nGenes_sort <- data.frame("nFeatRNA" = sort(nGenes), "X" = 1:length(nGenes))
scatPlot <- ggplot(nGenes_sort, aes(x = X, y = nFeatRNA)) + geom_point() + 
  ggtitle("HuaPt1") &
  geom_hline(yintercept = nFeatUpperN_hua) & geom_hline(yintercept = nFeatLowerN_norm)
print(scatPlot)
nGenes <- new.seven.split[["HuaPt2"]]@meta.data$nFeature_RNA
nGenes_sort <- data.frame("nFeatRNA" = sort(nGenes), "X" = 1:length(nGenes))
scatPlot <- ggplot(nGenes_sort, aes(x = X, y = nFeatRNA)) + geom_point() + 
  ggtitle("HuaPt2") &
  geom_hline(yintercept = nFeatUpperN_hua) & geom_hline(yintercept = nFeatLowerN_norm)
print(scatPlot)
nGenes <- huaPt3@meta.data$nFeature_RNA
nGenes_sort <- data.frame("nFeatRNA" = sort(nGenes), "X" = 1:length(nGenes))
scatPlot <- ggplot(nGenes_sort, aes(x = X, y = nFeatRNA)) + geom_point() + 
  ggtitle("HuaPt3") &
  geom_hline(yintercept = nFeatUpperN_hua) & geom_hline(yintercept = nFeatLowerN_norm)
print(scatPlot)
nGenes <- new.seven.split[["HuaPt4"]]@meta.data$nFeature_RNA
nGenes_sort <- data.frame("nFeatRNA" = sort(nGenes), "X" = 1:length(nGenes))
scatPlot <- ggplot(nGenes_sort, aes(x = X, y = nFeatRNA)) + geom_point() + 
  ggtitle("HuaPt4") &
  geom_hline(yintercept = nFeatUpperN_hua) & geom_hline(yintercept = nFeatLowerN_norm)
print(scatPlot)
dev.off()

# Plot filters based on mitochondrial content
# Setty mito plot
pdf(paste(dir, "Setty_mito", date,".pdf",sep=""))
vPlot <- VlnPlot(object = setty, features = c("percent.mito"), pt.size = 0, group.by = "orig.ident") +
  scale_y_continuous(limits = c(0,1)) & 
  geom_hline(yintercept = perMitoUpperN_norm)
print(vPlot)
dev.off()

# Zheng mito plot
pdf(paste(dir, "zheng_mito", date,".pdf",sep=""))
vPlot <- VlnPlot(object = zheng, features = c("percent.mito"), pt.size = 0, group.by = "orig.ident")+
  scale_y_continuous(limits = c(0,1)) & 
  geom_hline(yintercept = perMitoUpperN_norm)
print(vPlot)
dev.off()

# Hua mito plot
pdf(paste(dir, "Hua_mito", date,".pdf",sep=""))
vPlot <- VlnPlot(object = hua, features = c("percent.mito"), pt.size = 0, group.by = "orig.ident")+ 
  scale_y_continuous(limits = c(0,1)) & 
  geom_hline(yintercept = perMitoUpperN_norm)
print(vPlot)
dev.off()

#Filter based on cutoffs
setty.qc <- subset(x=setty, subset=nFeature_RNA > nFeatLowerN_norm & nFeature_RNA < nFeatUpperN_setty & percent.mito < perMitoUpperN_norm)
zheng.qc <- subset(x = zheng, subset = nFeature_RNA > nFeatLowerN_norm & nFeature_RNA < nFeatUpperN_zheng & percent.mito < perMitoUpperN_norm)
hua.qc <- subset(x=hua, subset=nFeature_RNA > nFeatLowerN_norm & nFeature_RNA < nFeatUpperN_hua & percent.mito < perMitoUpperN_norm)

#See how many we filter
dim(setty) #[1] 33694 41331
dim(setty.qc) #[1] 33694 25041
dim(zheng) #[1] 33694  9262
dim(zheng.qc) #[1] 33694  8799
dim(hua) #[1] 33694 32289
dim(hua.qc) #[1] 33694 29832

#Combine all of the normals together
healthy8.qc <- merge(setty.qc, zheng.qc)
healthy8.qc <- merge(healthy8.qc, hua.qc)

saveRDS(healthy8.qc, paste0(dir, 'postQC_healthy8', date, '.rds'))

############## SWITCH TO CMML SAMPLES ##################################

#Import the seurat object for the first ~40 samples (1 repeat so 39)
first39 <- readRDS(paste0(dir, 'first40/CMML_first40_UnprocessedCohort_SeuratObj_20200930.rds'))

first39 <- first40
levels(first39)
dim(first39)

#Add to metadata to separate normals from cmml/moffitt samples
first39[['orig.ident']] <- first39@active.ident
first39[['tech']] <- 'MCC'

#Get summary of all
summary(first39@meta.data$nFeature_RNA)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 17    1067    2418    2496    3654    9742

#Calculate mito percentage
mito.genes <- grep(pattern = "^MT-", first39@assays$RNA@counts@Dimnames[[1]], value = TRUE)
percent.mito <- Matrix::colSums(x=GetAssayData(object=first39, slot='counts')[mito.genes, ])/Matrix::colSums(x=GetAssayData(object=first39, slot='counts'))
first39[['percent.mito']] <- percent.mito

#Set standard cutoffs across all CMML samples (same as nFeatLowerN_norm???, diff mito, diff upperFeatN)
nFeatUpper_CMML <- mean(first39@meta.data$nFeature_RNA, na.rm=TRUE) + 2*sd(first39@meta.data$nFeature_RNA, na.rm=TRUE) # 5808.62 (consistent with Meghan)
nFeatLower_CMML <- 450 
perMitoUpper_CMML <- 0.25 

# Plot the violin and scatter plot of number of features with the cutoffs shown
# CMML number of features
pdf(paste(dir, "cmml_All_nFeatRNA", date, ".pdf",sep=""), width = 11, height = 3)
#Ignore message about scale, plotting range to be same across all samples
vPlot <- VlnPlot(object = first39, features = c("nFeature_RNA"), pt.size = 0, combine = T, group.by = "orig.ident") + 
  scale_y_continuous(limits = c(0,9000)) & 
  geom_hline(yintercept = nFeatUpper_CMML) & geom_hline(yintercept = nFeatLower_CMML)
print(vPlot)
dev.off()

# Plot filters based on mitochondrial content
# CMML mito plot
pdf(paste(dir, "CMML_mito", date,".pdf",sep=""), width = 10, height = 3)
vPlot <- VlnPlot(object = first39, features = c("percent.mito"), pt.size = 0, group.by = "orig.ident") +
  scale_y_continuous(limits = c(0,1)) & 
  geom_hline(yintercept = perMitoUpper_CMML)
print(vPlot)
dev.off()

# subset data based on cutoffs
first39.qc <- subset(x=first39, subset=nFeature_RNA > nFeatLower_CMML & nFeature_RNA < nFeatUpper_CMML & percent.mito < perMitoUpper_CMML)

dim(first39) #[1]  33694 182189
dim(first39.qc) #[1]  33694 137578
# # Keep about 75% of cells, 25% fall outside bounds

saveRDS(first39.qc, paste0(dir, 'postQC_CMML39', date, '.rds'))

##### Merge both post-qc and combine to give us our seurat object for future use
combined_39_8 <- merge(healthy8.qc, first39.qc)
saveRDS(combined_39_8, paste0(dir, 'postQC_CMML39+healthy8', date, '.rds'))


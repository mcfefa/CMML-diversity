library(Seurat)
library(Matrix)
library(stringr)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(RColorBrewer)
library(gplots)

#sessionInfo() output after restarting R (06-10-2021)
#R version 4.0.2 (2020-06-22)
#Platform: x86_64-apple-darwin17.0 (64-bit)
#Running under: macOS Catalina 10.15.7

#Matrix products: default
#BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
#LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

#locale:
#  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] gplots_3.1.1       RColorBrewer_1.1-2 forcats_0.5.1      dplyr_1.0.5        purrr_0.3.4       
#[6] readr_1.4.0        tidyr_1.1.3        tibble_3.1.1       tidyverse_1.3.1    patchwork_1.1.1   
#[11] ggplot2_3.3.3      stringr_1.4.0      SeuratObject_4.0.1 Seurat_4.0.1       Matrix_1.3-3      

#loaded via a namespace (and not attached):
#  [1] Rtsne_0.15             colorspace_2.0-1       deldir_0.2-10          ellipsis_0.3.2        
#[5] ggridges_0.5.3         XVector_0.28.0         fs_1.5.0               GenomicRanges_1.40.0  
#[9] rstudioapi_0.13        spatstat.data_2.1-0    leiden_0.3.7           listenv_0.8.0         
#[13] ggrepel_0.9.1          lubridate_1.7.10       fansi_0.4.2            xml2_1.3.2            
#[17] codetools_0.2-18       splines_4.0.2          polyclip_1.10-0        jsonlite_1.7.2        
#[21] broom_0.7.6            ica_1.0-2              dbplyr_2.1.1           cluster_2.1.2         
#[25] png_0.1-7              uwot_0.1.10            shiny_1.6.0            sctransform_0.3.2     
#[29] spatstat.sparse_2.0-0  compiler_4.0.2         httr_1.4.2             backports_1.2.1       
#[33] assertthat_0.2.1       fastmap_1.1.0          lazyeval_0.2.2         cli_2.5.0             
#[37] later_1.2.0            htmltools_0.5.1.1      tools_4.0.2            igraph_1.2.6          
#[41] gtable_0.3.0           glue_1.4.2             GenomeInfoDbData_1.2.3 RANN_2.6.1            
#[45] reshape2_1.4.4         Rcpp_1.0.6             scattermore_0.7        cellranger_1.1.0      
#[49] vctrs_0.3.8            nlme_3.1-152           lmtest_0.9-38          globals_0.14.0        
#[53] rvest_1.0.0            mime_0.10              miniUI_0.1.1.1         lifecycle_1.0.0       
#[57] irlba_2.3.3            gtools_3.8.2           goftest_1.2-2          future_1.21.0         
#[61] zlibbioc_1.34.0        MASS_7.3-54            zoo_1.8-9              scales_1.1.1          
#[65] spatstat.core_2.1-2    hms_1.0.0              promises_1.2.0.1       spatstat.utils_2.1-0  
#[69] parallel_4.0.2         reticulate_1.20        pbapply_1.4-3          gridExtra_2.3         
#[73] rpart_4.1-15           stringi_1.5.3          S4Vectors_0.26.1       caTools_1.18.2        
#[77] BiocGenerics_0.34.0    GenomeInfoDb_1.24.2    rlang_0.4.11           pkgconfig_2.0.3       
#[81] matrixStats_0.58.0     bitops_1.0-7           lattice_0.20-44        ROCR_1.0-11           
#[85] tensor_1.5             htmlwidgets_1.5.3      cowplot_1.1.1          tidyselect_1.1.1      
#[89] parallelly_1.25.0      RcppAnnoy_0.0.18       plyr_1.8.6             magrittr_2.0.1        
#[93] R6_2.5.0               IRanges_2.22.2         generics_0.1.0         DBI_1.1.1             
#[97] haven_2.4.1            pillar_1.6.0           withr_2.4.2            mgcv_1.8-35           
#[101] fitdistrplus_1.1-3     survival_3.2-11        abind_1.4-5            RCurl_1.98-1.3        
#[105] future.apply_1.7.0     modelr_0.1.8           crayon_1.4.1           KernSmooth_2.23-20    
#[109] utf8_1.2.1             spatstat.geom_2.1-0    plotly_4.9.3           readxl_1.3.1          
#[113] grid_4.0.2             data.table_1.14.0      reprex_2.0.0           digest_0.6.27         
#[117] xtable_1.8-4           httpuv_1.6.0           stats4_4.0.2           munsell_0.5.0         
#[121] viridisLite_0.4.0 

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


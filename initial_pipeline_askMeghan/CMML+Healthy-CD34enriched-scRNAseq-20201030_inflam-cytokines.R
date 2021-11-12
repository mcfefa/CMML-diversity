## Run interactively on the Cluster
## qsub -I -q bigmemQ -l nodes=1:ppn=16 -l walltime=120:00:00 
## module load python/3.7.2
## module load R/3.6.0
## export R_LIBS_USER=$HOME/apps/R:$R_LIB_USER     #<---- if installing in home directory --- didn't do
## R

## ADAPTED: CMML+Healthy-CD34enriched-scRNAseq-20200731_MALAT1+classify.R to look at cyokines expression across UMAP for PQ Grant

#### 10/30/2020
## on compute-3-5 with 14 ppns

library(future)
library(ggplot2)
library(bigmemory)
library(Seurat)
library(data.table)
library(dplyr)
library(rdetools)
#library(clustree)
sessionInfo()

# R version 3.6.0 (2019-04-26)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS release 6.10 (Final)

# Matrix products: default
# BLAS:   /share/apps/R-3.6.0/lib64/R/lib/libRblas.so
# LAPACK: /share/apps/R-3.6.0/lib64/R/lib/libRlapack.so

# Random number generation:
#  RNG:     Mersenne-Twister
#  Normal:  Inversion
#  Sample:  Rounding

# locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
#  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
#  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
#  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base

# other attached packages:
# [1] rdetools_1.0      dplyr_0.8.3       data.table_1.12.8 clustree_0.4.2
# [5] ggraph_2.0.0      Seurat_3.1.2      bigmemory_4.5.36  ggplot2_3.2.1
# [9] future_1.16.0

# loaded via a namespace (and not attached):
#   [1] TH.data_1.0-10      bigmemory.sri_0.1.3 Rtsne_0.15
#   [4] colorspace_1.4-1    ggridges_0.5.2      leiden_0.3.3
#   [7] listenv_0.8.0       farver_2.0.3        npsurv_0.4-0
#  [10] graphlayouts_0.5.0  ggrepel_0.8.1       mvtnorm_1.0-12
#  [13] codetools_0.2-16    splines_3.6.0       R.methodsS3_1.7.1
#  [16] mnormt_1.5-5        lsei_1.2-0          TFisher_0.2.0
#  [19] polyclip_1.10-0     jsonlite_1.6        ica_1.0-2
#  [22] cluster_2.1.0       png_0.1-7           R.oo_1.23.0
#  [25] uwot_0.1.5          ggforce_0.3.1       sctransform_0.2.1
#  [28] compiler_3.6.0      httr_1.4.1          assertthat_0.2.1
#  [31] Matrix_1.2-18       lazyeval_0.2.2      tweenr_1.0.1
#  [34] htmltools_0.4.0     tools_3.6.0         rsvd_1.0.2
#  [37] igraph_1.2.4.2      gtable_0.3.0        glue_1.3.1
#  [40] RANN_2.6.1          reshape2_1.4.3      Rcpp_1.0.3
#  [43] Biobase_2.44.0      vctrs_0.2.2         multtest_2.40.0
#  [46] gdata_2.18.0        ape_5.3             nlme_3.1-143
#  [49] gbRd_0.4-11         lmtest_0.9-37       stringr_1.4.0
#  [52] globals_0.12.5      lifecycle_0.1.0     irlba_2.3.3
#  [55] gtools_3.8.1        MASS_7.3-51.5       zoo_1.8-7
#  [58] scales_1.1.0        tidygraph_1.1.2     parallel_3.6.0
#  [61] sandwich_2.5-1      RColorBrewer_1.1-2  reticulate_1.14
#  [64] pbapply_1.4-2       gridExtra_2.3       stringi_1.4.5
#  [67] mutoss_0.1-12       plotrix_3.7-7       caTools_1.18.0
#  [70] BiocGenerics_0.30.0 bibtex_0.4.2.2      Rdpack_0.11-1
#  [73] SDMTools_1.1-221.2  rlang_0.4.4         pkgconfig_2.0.3
#  [76] bitops_1.0-6        lattice_0.20-38     ROCR_1.0-7
#  [79] purrr_0.3.3         htmlwidgets_1.5.1   cowplot_1.0.0
#  [82] tidyselect_1.0.0    RcppAnnoy_0.0.14    plyr_1.8.5
#  [85] magrittr_1.5        R6_2.4.1            gplots_3.0.1.2
#  [88] multcomp_1.4-12     pillar_1.4.3        withr_2.1.2
#  [91] sn_1.5-4            fitdistrplus_1.0-14 survival_3.1-8
#  [94] tibble_2.1.3        future.apply_1.4.0  tsne_0.1-3
#  [97] crayon_1.3.4        KernSmooth_2.23-16  plotly_4.9.1
# [100] viridis_0.5.1       grid_3.6.0          metap_1.3
# [103] digest_0.6.23       tidyr_1.0.2         numDeriv_2016.8-1.1
# [106] R.utils_2.9.2       RcppParallel_4.4.4  stats4_3.6.0
# [109] munsell_0.5.0       viridisLite_0.3.0

plan("multiprocess", workers=14)
setwd("/share/lab_padron/Meghan/scRNAseq/CMML/")
#setwd("/Volumes/lab_padron/Meghan/scRNAseq/CMML/analysis")

# #options(future.globals.maxSize=134217728000) ### for 128GB * 1024*1024
# options(future.globals.maxSize=536870912000) ###  512GB: 536870912000
options(future.globals.maxSize=402653184000) ###  384GB: 536870912000
# #options(future.globals.maxSize=171966464000) ### for 164GB * 1024*1024

filedir <- "./analysis/PQ3_H7+C31-Cohort_CD34only+33kGenes_"
date <- "_2020-10-30"

FullSetMultiRes <- readRDS("./analysis/09_C31+N7_CD34only_regularPipeline/Healthy-7pts+CMML-31pts-Cohort_CD34only+33kGenes_SeuratObject_thruStep9_res0.05-as-active_ReScaledAllGenes_2020-04-29.rds")

### VISUALIZING MALAT 1 EXPRESSION IN THESE CELLS

cytokines = c("BMP4","BTC","C3","CAMP","CAT","CCL13","CCL14","CCL18","CCL24","CCL8","CMTM1","CMTM5","CTGF","CXCL12","CXCL2","CXCL3","CYR61","DEFA3","EBI3","EPGN","EREG","FGF10","FGF13","FGF16","FGF20","GDF3","GMFB","GMFG","IGF1","IL15","IL1A","IL1B","IL33","IL5","IL6","IL6ST","IL8","LTB","LTBP3","NAMPT","NPY","NRG1","NRG4","OGN","OSGIN1","OXT","PF4V1","POMC","RLN2","S100A6","SBDS","SCT","SECTM1","SEMA4B","SEMA6D","SPP1","TFGB1","TNC","TNFSF10","TNFSF13","TNFSF4","VEGFB","VIP")


cytokines2 = c("CCL13","CCL14","CCL18","CCL24","CCL8","CXCL12","CXCL2","CXCL3","CXCL8","GMFB","GMFG","IGF1","IL15","IL1A","IL1B","IL33","IL5","IL6","IL6ST","S100A6","SCT","TNFSF10","TNFSF13","TNFSF4","VEGFB")

shortList = c("IL5","IL6","IL6ST","CXCL8")

## cytokines
pdf(paste(filedir,"VlnPlot_shortList",date,".pdf",sep=""),width=16,height=8)
dPlot <- VlnPlot(FullSetMultiRes, features = shortList, group.by='orig.ident', pt.size=0)
print(dPlot)
dev.off()

pdf(paste(filedir,"RidgePlot_shortList",date,".pdf",sep=""),width=16,height=8)
dPlot <- RidgePlot(FullSetMultiRes, features = shortList, group.by='orig.ident', pt.size=0)
print(dPlot)
dev.off()

pdf(paste(filedir,"RidgePlot_shortList_byCluster",date,".pdf",sep=""),width=16,height=8)
dPlot <- RidgePlot(FullSetMultiRes, features = shortList, group.by='cluster', pt.size=0)
print(dPlot)
dev.off()

pdf(paste(filedir,"VlnPlot_shortList_byCluster",date,".pdf",sep=""),width=16,height=8)
dPlot <- VlnPlot(FullSetMultiRes, features = shortList, group.by='cluster', pt.size=0)
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_shortList",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = shortList)
print(dPlot)
dev.off()


### INDIVIDUAL FEATURE PLOTS

cytokines = c("","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","SEMA6D","SPP1","TFGB1","TNC","TNFSF10","TNFSF13","TNFSF4","VEGFB","VIP")

pdf(paste(filedir,"FeaturePlot_SEMA4B",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("SEMA4B"))
print(dPlot)
dev.off()
#<---
pdf(paste(filedir,"FeaturePlot_SECTM1",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("SECTM1"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_SCT",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("SCT"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_SBDS",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("SBDS"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_S100A6",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("S100A6"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_RLN2",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("RLN2"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_POMC",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("POMC"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_PF4V1",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("PF4V1"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_OXT",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("OXT"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_OSGIN1",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("OSGIN1"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_OGN",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("OGN"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_NRG4",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("NRG4"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_NRG1",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("NRG1"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_NPY",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("NPY"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_NAMPT",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("NAMPT"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_LTBP3",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("LTBP3"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_LTB",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("LTB"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_IL6ST",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("IL6ST"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_IL33",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("IL33"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_IL1B",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("IL1B"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_IL1A",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("IL1A"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_IGF1",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("IGF1"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_GMFG",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("GMFG"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_GMFB",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("GMFB"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_GDF3",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("GDF3"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_FGF20",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("FGF20"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_FGF16",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("FGF16"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_FGF13",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("FGF13"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_FGF10",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("FGF10"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_EREG",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("EREG"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_EPGN",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("EPGN"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_EBI3",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("EBI3"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_DEFA3",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("DEFA3"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_CYR61",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("CYR61"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_CXCL3",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("CXCL3"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_CXCL2",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("CXCL2"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_CXCL12",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("CXCL12"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_CTGF",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("CTGF"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_CMTM5",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("CMTM5"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_CMTM18",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("CMTM1"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_CCL8",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("CCL8"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_CCL24",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("CCL24"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_CCL18",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("CCL18"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_CCL14",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("CCL14"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_CCL13",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("CCL13"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_CAT",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("CAT"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_CAMP",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("CAMP"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_C3",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("C3"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_BTC",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("BTC"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_BMP4",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("BMP4"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_IL8",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("CXCL8"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_IL5",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("IL5"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_IL6",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("IL6"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_IL15",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("IL15"))
print(dPlot)
dev.off()

pdf(paste(filedir,"FeaturePlot_IL5RA",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes, features = c("IL5RA"))
print(dPlot)
dev.off()


## try converting to loom file to open in python
# https://satijalab.org/seurat/v3.1/conversion_vignette.html
# pbmc.loom <- as.loom(pbmc, filename = "../output/pbmc3k.loom", verbose = FALSE)
# pbmc.loom

## Error in if (size == Inf) { : missing value where TRUE/FALSE needed
## troublshooting error in making loom object: https://github.com/satijalab/seurat/issues/1451
install.packages("devtools") 
devtools::install_github(repo = 'mojaveazure/loomR', ref = 'develop') 
devtools::install_github(repo = "hhoeflin/hdf5r")

library(loomR)

## second troubleshooting attempt: https://github.com/mojaveazure/loomR/issues/40
obj <- FullSetMultiRes
for(j in 1:ncol(obj@meta.data)){
	if(is.factor(obj@meta.data[,j]) == T){
	obj@meta.data[,j] = as.character(obj@meta.data[,j]) # Force the variable type to be character
	obj@meta.data[,j][is.na(obj@meta.data[,j])] <- "N.A"
}
if(is.character(obj@meta.data[,j]) == T){
	obj@meta.data[,j][is.na(obj@meta.data[,j])] <- "N.A"
 }
}

FullSetMultiRes.loom <- as.loom(obj, filename = "./analysis/C31+N7_Cohort_Exported_att3_2020-07-31.loom", verbose = FALSE)

FullSetMultiRes.loom
# Filename: /share/lab_padron/Meghan/scRNAseq/CMML/analysis/C31+N7_Cohort_Exported_att3_2020-07-31.loom
# Access type: H5F_ACC_RDWR
# Attributes: version, chunks, LOOM_SPEC_VERSION, assay, last_modified
# Listing:
#        name    obj_type   dataset.dims dataset.type_class
#   col_attrs   H5I_GROUP           <NA>               <NA>
#  col_graphs   H5I_GROUP           <NA>               <NA>
#      layers   H5I_GROUP           <NA>               <NA>
#      matrix H5I_DATASET 174714 x 33694          H5T_FLOAT
#   row_attrs   H5I_GROUP           <NA>               <NA>
#  row_graphs   H5I_GROUP           <NA>               <NA>

FullSetMultiRes.loom$close_all()

### next week pick up by trying to load this into python and pick up the scoring pipeline from that nature paper

#<-------- HERE 


# CMML SAMPLE IDS
# First 8: LTB3966, LTB4121, LTB5109 (6169), 4J003, 4K001, 4Q001, 5E001, 5H001
# Second 8: SF100109106293, SF100109111451, SF100109110236, SF14040100158,
#           SF14060200025, SF12062800475, SF14072200012, SF13061200056
# Third 8: 4S001, 2V001, SF14101000049, SF16112900158, 6AE001, 6AC001, 
#          6AD001, SF100109101914
# Fourth 8: SF12042500035, SF12091900043, SF12092600014, SF14031800065, 
#           SF14050700419, SF16026800045, SF16072200003, SF16112300029

##### Step 1: Load all the data & Quality Control for CMML

# CMMLpts <- readRDS("./analysis/CMML_first32_UnprocessedCohort_SeuratObj_20191223.rds")
# dim(CMMLpts)
# # [1]  33694 145320

# levels(CMMLpts)
# # [1] "2V001"          "4J003"          "4K001"          "4Q001"
# # [5] "4S001"          "5E001"          "5H001"          "6AC001"
# # [9] "6AD001"         "6AE001"         "LTB3966"        "LTB4121"
# #[13] "LTB6169"        "SF100109101914" "SF100109106293" "SF100109110236"
# #[17] "SF100109111451" "SF12042500035"  "SF12062800475"  "SF12091900043"
# #[21] "SF12092600014"  "SF13061200056"  "SF14031800065"  "SF14040100158"
# #[25] "SF14050700419"  "SF14060200025"  "SF14072200012"  "SF14101000049"
# #[29] "SF16026800045"  "SF16072200003"  "SF16112300029"  "SF16112900158"

# ## realized these are from the same patient at the same timepoint, so nameing them the same
# CMMLptsU <- RenameIdents(CMMLpts, 'SF12091900043' = '6AD001')

# levels(CMMLptsU)
# # [1] "6AD001"         "2V001"          "4J003"          "4K001"         
# # [5] "4Q001"          "4S001"          "5E001"          "5H001"         
# # [9] "6AC001"         "6AE001"         "LTB3966"        "LTB4121"       
# #[13] "LTB6169"        "SF100109101914" "SF100109106293" "SF100109110236"
# #[17] "SF100109111451" "SF12042500035"  "SF12062800475"  "SF12092600014" 
# #[21] "SF13061200056"  "SF14031800065"  "SF14040100158"  "SF14050700419" 
# #[25] "SF14060200025"  "SF14072200012"  "SF14101000049"  "SF16026800045" 
# #[29] "SF16072200003"  "SF16112300029"  "SF16112900158" 

# rm(CMMLpts)
# CMMLpts <- CMMLptsU
# rm(CMMLptsU)

# date = "_2020-03-21"

# mito.genes <- grep(pattern = "^MT-", CMMLpts@assays$RNA@counts@Dimnames[[1]], value = TRUE)
# percent.mito <- Matrix::colSums(x=GetAssayData(object=CMMLpts, slot='counts')[mito.genes, ])/Matrix::colSums(x=GetAssayData(object=CMMLpts, slot='counts'))
# CMMLpts[['percent.mito']] <- percent.mito

# ## Descriptive stats:
# summary(CMMLpts@meta.data$nFeature_RNA)
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# #   17    1155    2499    2559    3711    9742
# mean(CMMLpts@meta.data$nFeature_RNA, na.rm=TRUE) + 2*sd(CMMLpts@meta.data$nFeature_RNA, na.rm=TRUE) # 5873.477
# nFeatUpper <- mean(CMMLpts@meta.data$nFeature_RNA, na.rm=TRUE) + 2*sd(CMMLpts@meta.data$nFeature_RNA, na.rm=TRUE) 
# nFeatLower <- 200 
# perMitoUpper <- 0.25 

# # subset data based on cutoffs
# CMMLpts <- subset(x=CMMLpts, subset=nFeature_RNA > nFeatLower & nFeature_RNA < nFeatUpper & percent.mito < perMitoUpper)

# dim(CMMLpts)
# # [1]  33694 114939
# # Keep about 80% of cells, 20% fall outside bounds
# # not changes to the metrics with renaming (as expected)

# ##### Step 2: Load all the data & Quality Control for Normal & Merge of Datasets
# ## ================ LOGIC/WORK-FLOW ============================
# ### need to create raw for all the different normals - create different rds files for different cohorts
# ### then combine them
# ### apply cutoffs uniformly across all patient samples
# ### then merge with CMML patients

# HealthyFirst4 <- readRDS("./analysis/Normal_CD34only_UnprocessedCohort_SeuratObj_20190920.rds")

# HuaPt1.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/Hua-CD34-dataset/BM1/GRCh38/")
# HuaPt1 <- CreateSeuratObject(counts=HuaPt1.data, project="HuaPt1")
# HuaPt1 <- RenameCells(object=HuaPt1, add.cell.id="HuaPt1")
# dim(HuaPt1)
# # [1] 33694 10787
# HuaPt2.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/Hua-CD34-dataset/BM2/GRCh38/")
# HuaPt2 <- CreateSeuratObject(counts=HuaPt2.data, project="HuaPt2")
# HuaPt2 <- RenameCells(object=HuaPt2, add.cell.id="HuaPt2")
# dim(HuaPt2)
# # [1] 33694 11256
# # **** removed this patient because of previous runs, this was far outside the range of other diversity values ****
# #HuaPt3.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/Hua-CD34-dataset/BM3/GRCh38/")
# #HuaPt3 <- CreateSeuratObject(counts=HuaPt3.data, project="HuaPt3")
# #HuaPt3 <- RenameCells(object=HuaPt3, add.cell.id="HuaPt3")
# #dim(HuaPt3)
# ## [1] 33694  5131
# HuaPt4.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/Hua-CD34-dataset/BM4/GRCh38/")
# HuaPt4 <- CreateSeuratObject(counts=HuaPt4.data, project="HuaPt4")
# HuaPt4 <- RenameCells(object=HuaPt4, add.cell.id="HuaPt4")
# dim(HuaPt4)
# # [1] 33694  5115

# ### cut down extra genes introduced by Setty dataset
# dim(HealthyFirst4)
# # [1] 74882 50593
# a <- HealthyFirst4@assays$RNA@counts@Dimnames[[1]] # length(a) [1] 74882
# b <- HuaPt1@assays$RNA@counts@Dimnames[[1]] # length(b) [1] 33694
# c <- intersect(a,b) # length(c) [1] 33694
# HealthyFirst4CUT <- HealthyFirst4[c,]
# dim(HealthyFirst4CUT)
# # [1] 33694 50593

# HealthyA <- merge(x=HealthyFirst4CUT, y=HuaPt1, add.cell.id2 = "HuaPt1", merge.data=TRUE, project = "NormalCohort")
# HealthyB <- merge(x=HealthyA, y=HuaPt2,  add.cell.id2 = "HuaPt2", merge.data=TRUE, project = "NormalCohort")
# Healthy7Pts <- merge(x=HealthyB, y=HuaPt4,  add.cell.id2 = "HuaPt4", merge.data=TRUE, project = "NormalCohort")

# dim(Healthy7Pts)
# # [1] 33694 77751

# saveRDS(Healthy7Pts, file = "./analysis/Normal_CD34only_Zheng+Setty+Hua_UnprocessedCohort_SeuratObj_20200321.rds")

# mito.genes <- grep(pattern = "^MT-", Healthy7Pts@assays$RNA@counts@Dimnames[[1]], value = TRUE)
# percent.mito <- Matrix::colSums(x=GetAssayData(object=Healthy7Pts, slot='counts')[mito.genes, ])/Matrix::colSums(x=GetAssayData(object=Healthy7Pts, slot='counts'))
# Healthy7Pts[['percent.mito']] <- percent.mito

# ## Descriptive stats:
# summary(Healthy7Pts@meta.data$nFeature_RNA)
# ### OLD - 8 patients
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# #   12     908    1436    1514    2046    7412 
# # upper cutoff: # [1] 3281.558
# ### 7 patients
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# #     12     896    1455    1523    2078    7412

# mean(Healthy7Pts@meta.data$nFeature_RNA, na.rm=TRUE) + 2*sd(Healthy7Pts@meta.data$nFeature_RNA, na.rm=TRUE) # [1] 3326.408
# nFeatUpperN <- mean(Healthy7Pts@meta.data$nFeature_RNA, na.rm=TRUE) + 2*sd(Healthy7Pts@meta.data$nFeature_RNA, na.rm=TRUE) 
# nFeatLowerN <- 200
# perMitoUpperN <- 0.05 

# # We filter out cells that have unique feature counts over 2*std+mean or less than 200
# Healthy7Pts <- subset(x=Healthy7Pts, subset=nFeature_RNA > nFeatLowerN & nFeature_RNA < nFeatUpperN & percent.mito < perMitoUpperN)

# dim(Healthy7Pts)
# #### OLD 8 patient numbers
# # [1] 33694 64567
# ## Keep about 77.9% of the original dataset
# #### After removing 1 normal: 
# # [1] 33694 59775

# saveRDS(Healthy7Pts, file = "./analysis/Normal_CD34only_Zheng+Setty+Hua_QCprocessedCohort_SeuratObj_20200321.rds")

# #### COMBINE COHORTS
# FullSet <- merge(x = Healthy7Pts, y = CMMLpts,  add.cell.id2 = "CMMLpts", merge.data=TRUE, project = "FullCohort")
# dim(FullSet)
# #### 32+8 patient cohort
# # [1]  33694 179506
# #### 31+7 patient cohort
# # [1]  33694 174714

# pdf(paste(filedir,"ViolinPlot_nFeatures",date,".pdf",sep=""))
# vPlot <- VlnPlot(object = FullSet, features = c("nFeature_RNA"))
# print(vPlot)
# dev.off() 

# pdf(paste(filedir,"ViolinPlot_nFeatures_noPts",date,".pdf",sep=""))
# vPlot <- VlnPlot(object = FullSet, features = c("nFeature_RNA"), pt.size=0)
# print(vPlot)
# dev.off() 

# pdf(paste(filedir,"ViolinPlot_nFeatures_noPts_noLeg",date,".pdf",sep=""))
# vPlot <- VlnPlot(object = FullSet, features = c("nFeature_RNA"), pt.size=0) + theme(legend.position = "none")
# print(vPlot)
# dev.off() 

# pdf(paste(filedir,"ViolinPlot_perMito",date,".pdf",sep=""))
# vPlot <- VlnPlot(object = FullSet, features = c("percent.mito"))
# print(vPlot)
# dev.off() 

# pdf(paste(filedir,"ViolinPlot_perMito_noPts",date,".pdf",sep=""))
# vPlot <- VlnPlot(object = FullSet, features = c("percent.mito"), pt.size=0)
# print(vPlot)
# dev.off() 

# pdf(paste(filedir,"ViolinPlot_perMito_noPts_noLeg",date,".pdf",sep=""))
# vPlot <- VlnPlot(object = FullSet, features = c("percent.mito"), pt.size=0) + theme(legend.position = "none")
# print(vPlot)
# dev.off() 

# pdf(paste(filedir,"ScatterPlot_RNAcount-vs-mito",date,".pdf",sep=""))
# fsPlot <- FeatureScatter(object=FullSet, feature1 = "nCount_RNA", feature2 = "percent.mito")
# print(fsPlot)
# dev.off()

# pdf(paste(filedir,"ScatterPlot_RNAcount-vs-mito_noLeg",date,".pdf",sep=""))
# fsPlot <- FeatureScatter(object=FullSet, feature1 = "nCount_RNA", feature2 = "percent.mito") + theme(legend.position = "none")
# print(fsPlot)
# dev.off()

# pdf(paste(filedir,"ScatterPlot_nGenes-vs-nRNA",date,".pdf",sep=""))
# fsPlot <- FeatureScatter(object=FullSet, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# print(fsPlot)
# dev.off()

# pdf(paste(filedir,"ScatterPlot_nGenes-vs-nRNA_noLeg",date,".pdf",sep=""))
# fsPlot <- FeatureScatter(object=FullSet, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")  + theme(legend.position = "none")
# print(fsPlot)
# dev.off()

# FullSet$original <- FullSet$orig.ident
# FullSet$orig.ident <- NULL
# FullSet$orig.ident <- FullSet@active.ident
# saveRDS(FullSet, file = paste(filedir,"SeuratObj_thrStep2_combined",date,".rds",sep=""))

# #### Step 3: Normalize
# FullSet <- NormalizeData(object=FullSet, normalization.method="LogNormalize", scale.factor=1e4)
# dim(FullSet)
# # [1]  33694 174714
# saveRDS(FullSet, file = paste(filedir,"SeuratObject_thruStep3",date,".rds",sep=""))

# ##### Step 4: Detection of variable features across the single cells
# FullSet <- FindVariableFeatures(object=FullSet, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf), do.plot=TRUE)
# length(x = VariableFeatures(object=FullSet))  
# ## number of variable features detected after filtering, normalization: 4210
# dim(FullSet)
# # [1]  33694 174714

# saveRDS(FullSet, file = paste(filedir,"SeuratObject_thruStep4",date,".rds",sep=""))

# ##### Step 5: Scaling the data and removing unwanted sources of variation 
# FullSet <- ScaleData(object=FullSet, vars.to.regress = c("nCount_RNA","percent.mito")) 
# saveRDS(FullSet, file = paste(filedir,"SeuratObject_thruStep5",date,".rds",sep=""))

# ##### Step 6: Perform linear dimensional reduction 
# FullSet <- RunPCA(object=FullSet, features=VariableFeatures(object=FullSet), verbose = TRUE, npcs=200)

# # Examine and visualize PCA results a few different ways
# print(x=FullSet[['pca']], dims = 1:5, nfeatures = 5, projected = FALSE)

# pdf(paste(filedir,"PCA-Dim-Loadings",date,".pdf",sep=""))
# vdlPlot <- VizDimLoadings(object=FullSet, dims = 1:2)
# print(vdlPlot)
# dev.off()

# pdf(paste(filedir,"PCA-Dim-Plot",date,".pdf",sep=""))
# dPlot <- DimPlot(object=FullSet)
# print(dPlot)
# dev.off()

# # ProjectDim scores each feature in the dataset (including features not included in the PCA) based on their  
# # correlation with the calculated components. Though we don't use this further here, it can be used to identify  
# # markers that are strongly correlated with cellular heterogeneity, but may not have passed through variable 
# # feature selection. The results of the projected PCA can be explored by setting `projected = TRUE` in the 
# # functions above
# FullSet <- ProjectDim(object=FullSet)

# pdf(paste(filedir,"PCA-Dim-Heatmaps-1thr12",date,".pdf",sep=""))
# dhPlot <- DimHeatmap(object=FullSet, dims = 1:9, cells = 500, balanced = TRUE)
# print(dhPlot)
# dev.off()

# saveRDS(FullSet, file = paste(filedir,"SeuratObject_thruStep6",date,".rds",sep=""))

# ###### Step 7: Find Clusters
# FullSet <- FindNeighbors(object=FullSet, dims = 1:100)
# FullSet <- FindClusters(object=FullSet, resolution = 0.4)

# saveRDS(FullSet, file = paste(filedir,"SeuratObject_thruStep7_Res0p4only",date,".rds",sep=""))

# file1 <- "./analysis/CMMLoutputDataNames_31C+7Npts_CD34only_res0p4_20200321.csv";
# file2 <- "./analysis/CMMLoutputData_31C+7Npts_CD34only_res0p4_20200321.csv";

# cat(FullSet@meta.data$RNA_snn_res.0.4, file=file2, sep=",\n")

# listNames <- FullSet@meta.data$orig.ident
# write.table(data.frame(listNames),
#           row.names=FALSE,
#           col.names = FALSE, 
#           file = file1,
#           sep=",")

# mydat1 <- read.csv(file2)
# mydat2 <- read.csv(file1)
# fulldat <- cbind(mydat2[1],mydat1[1])
# fulltab <- as.data.table(fulldat)
# # name table columns
# names(fulltab)[1] <- paste("UMI")
# names(fulltab)[2] <- paste("cluster")
# # group data based on clusters
# group_by(fulltab, cluster)
# # create a table counting unqiue UMIs/cells per cluster
# tabPerClus <- fulltab %>% group_by(cluster) %>% count()
# type <- sub("\\_.*","",fulltab$UMI)
# fulltab <- cbind(fulltab, type)

# divout <-  paste(filedir,"CellBreakdown_PerClusterPerType_res-0p4",date,".csv",sep="") 
# A <- fulltab %>% group_by(cluster) %>% count(type)
# write.csv(A, file=divout)

# range <- c(0, 0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6)
# for (res in range) {
#   FullSet <- FindNeighbors(object=FullSet, dims = 1:100)
#   FullSet <- FindClusters(object=FullSet, resolution = res, n.start=100)
# }

# pdf(paste(filedir,"clustree_res-scan",date,".pdf",sep=""), width=18, height=36)
# tPlot <- clustree(FullSet, prefix = 'RNA_snn_res.')
# print(tPlot)
# dev.off()

# saveRDS(FullSet, paste(filedir,"SeuratObject_thruStep7_multiRes",date,".rds",sep=""))

# ## trying to implement Leiden clustering 

# FullSet2 <- FullSet
# FullSet2 <- FindNeighbors(object=FullSet2, dims = 1:100)
# FullSet2 <- FindClusters(object=FullSet2, resolution = 0.025)

# # 0%   10   20   30   40   50   60   70   80   90   100%
# # [----|----|----|----|----|----|----|----|----|----|
# # **************************************************|
# # Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck

# # Number of nodes: 174714
# # Number of edges: 7745341

# # Running Louvain algorithm...
# # Maximum modularity in 10 random starts: 0.9896
# # Number of communities: 15
# # Elapsed time: 51 seconds
# # 4 singletons identified. 11 final clusters.


# pdf(paste(filedir,"UMAP_res0.025_byClusters",date,".pdf",sep=""), width=12, height=8)
# dPlot <- DimPlot(FullSet2, reduction.use="umap") + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.025_clustersByType",date,".pdf",sep=""), width=12, height=8)
# dPlot <- DimPlot(FullSet2, reduction.use="umap", group.by='orig.ident') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.025_clustersByType_noLegend",date,".pdf",sep=""), width=12, height=8)
# dPlot <- DimPlot(FullSet2, reduction.use="umap", group.by='orig.ident') + theme(legend.position = "none") + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# ## failed because didn't have leidenalg installed/needs newer R
# ## will try again once FullSet3 finishes
# # to try exited
# module load anaconda/3
# conda activate /home/4467528/MegPython
# # conda forge installed leidenalg
# module load R/3.6.0
# R
# # run first few to add libraries, multiprocess, filedir, etc
# FullSet2 <- readRDS(paste(filedir,"SeuratObject_thruStep8_UMAP_2020-03-23.rds",sep=""))
# #FullSet2 <- FindNeighbors(object=FullSet2, dims = 1:100)
# library(reticulate)
# #use_condaenv("/home/4467528/MegPython")
# reticulate::py_install(packages ='igraph')
# reticulate::py_install(packages ='leidenalg')
# if (!requireNamespace("devtools"))
#     install.packages("devtools")
# devtools::install_github("TomKellyGenetics/leiden")
# install.packages("leiden")
# library("leiden")
# FullSet2 <- FindClusters(object=FullSet2, algorithm=4)

# use_condaenv(condaenv="Renv", conda="/home/4467528/MegPython")
# ## can't get this to work, going to put on the backburner now

# # tried on a separate node => only 8 nodes, 384 GB
# # filedir = "./analysis/Healthy-7pts+CMML-31pts-Cohort_CD34only+33kGenes_"
# # date <- "_2020-03-23"
# # FullSet3 <- readRDS(paste(filedir,"SeuratObject_thruStep8_UMAP_2020-03-23.rds",sep=""))
# FullSet3 <- FindNeighbors(object=FullSet3, k.param=150)
# FullSet3 <- FindClusters(object=FullSet3, resolution = 0.025)
# date <- "_2020-03-25"

# # 0%   10   20   30   40   50   60   70   80   90   100%
# # [----|----|----|----|----|----|----|----|----|----|
# # **************************************************|
# # Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck

# # Number of nodes: 174714
# # Number of edges: 41982074

# # Running Louvain algorithm...
# # Maximum modularity in 10 random starts: 0.9829
# # Number of communities: 6
# # Elapsed time: 415 seconds

# pdf(paste(filedir,"UMAP_k150_res0.025_byClusters",date,".pdf",sep=""), width=12, height=8)
# dPlot <- DimPlot(FullSet3, reduction.use="umap") + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_k150_res0.025_clustersByType",date,".pdf",sep=""), width=12, height=8)
# dPlot <- DimPlot(FullSet3, reduction.use="umap", group.by='orig.ident') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_k150_res0.025_clustersByType_noLegend",date,".pdf",sep=""), width=12, height=8)
# dPlot <- DimPlot(FullSet3, reduction.use="umap", group.by='orig.ident') + theme(legend.position = "none") + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# file1 <- "./analysis/CMMLoutputDataNames_31C+7Npts_CD34only_res0.025_k150_20200325.csv";
# file2 <- "./analysis/CMMLoutputData_31C+7Npts_CD34only_res0.025_k150_20200325.csv";

# cat(FullSet3@meta.data$RNA_snn_res.0.025, file=file2, sep=",\n")

# listNames <- FullSet3@meta.data$orig.ident
# write.table(data.frame(listNames),
#           row.names=FALSE,
#           col.names = FALSE, 
#           file = file1,
#           sep=",")

# mydat1 <- read.csv(file2)
# mydat2 <- read.csv(file1)
# fulldat <- cbind(mydat2[1],mydat1[1])
# fulltab <- as.data.table(fulldat)
# # name table columns
# names(fulltab)[1] <- paste("UMI")
# names(fulltab)[2] <- paste("cluster")
# # group data based on clusters
# group_by(fulltab, cluster)
# # create a table counting unqiue UMIs/cells per cluster
# tabPerClus <- fulltab %>% group_by(cluster) %>% count()
# type <- sub("\\_.*","",fulltab$UMI)
# fulltab <- cbind(fulltab, type)

# divout <-  paste(filedir,"CellBreakdown_PerClusterPerType_res-0.025_k150",date,".csv",sep="") 
# A <- fulltab %>% group_by(cluster) %>% count(type)
# write.csv(A, file=divout)

# range <- c(0, 0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6)
# for (res in range) {
#   FullSet3 <- FindNeighbors(object=FullSet3, k.param=150)
#   FullSet3 <- FindClusters(object=FullSet3, resolution = res, n.start=100)
# }

# pdf(paste(filedir,"clustree_res-scan_k150",date,".pdf",sep=""), width=18, height=36)
# tPlot <- clustree(FullSet3, prefix = 'RNA_snn_res.')
# print(tPlot)
# dev.off()

# saveRDS(FullSet3, paste(filedir,"SeuratObject_thruStep8_multiRes-k150",date,".rds",sep=""))

# date <- "_2020-03-26"
# FullSet3 <- RunUMAP(object=FullSet3, reduction = "pca", dims = 1:100, n.neighbors=150)

# # Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
# # To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
# # This message will be shown once per session
# # 14:33:42 UMAP embedding parameters a = 0.9922 b = 1.112
# # 14:33:44 Read 174714 rows and found 100 numeric columns
# # 14:33:44 Using Annoy for neighbor search, n_neighbors = 150
# # 14:33:44 Building Annoy index with metric = cosine, n_trees = 50
# # 0%   10   20   30   40   50   60   70   80   90   100%
# # [----|----|----|----|----|----|----|----|----|----|
# # **************************************************|
# # 14:36:48 Writing NN index file to temp file /scratch/2494762.colossus.local/RtmpU6yHoW/filead2e9427df1
# # 14:36:48 Searching Annoy index using 8 threads, search_k = 15000
# # 14:38:10 Annoy recall = 100%
# # 14:38:13 Commencing smooth kNN distance calibration using 8 threads
# # 14:40:09 Initializing from normalized Laplacian + noise
# # 14:43:22 Commencing optimization for 200 epochs, with 36548038 positive edges
# # 0%   10   20   30   40   50   60   70   80   90   100%
# # [----|----|----|----|----|----|----|----|----|----|
# # **************************************************|
# # 14:58:22 Optimization finished

# date <- "_2020-03-27"
# FullSet3 <- FindClusters(object=FullSet3, resolution = 0.025)

# pdf(paste(filedir,"UMAP_v3_k150_res0.025_byClusters",date,".pdf",sep=""), width=12, height=8)
# dPlot <- DimPlot(FullSet3, reduction.use="umap") + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_v3_k150_res0.025_clustersByType",date,".pdf",sep=""), width=12, height=8)
# dPlot <- DimPlot(FullSet3, reduction.use="umap", group.by='orig.ident') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()


# # range <- c(0, 0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6)

# file1 <- "./analysis/CMMLoutputDataNames_31C+7Npts_CD34only_res0.6_k150_20200327.csv";
# file2 <- "./analysis/CMMLoutputData_31C+7Npts_CD34only_res0.6_k150_20200327.csv";

# cat(FullSet3@meta.data$RNA_snn_res.0.6, file=file2, sep=",\n")

# listNames <- FullSet3@meta.data$orig.ident
# write.table(data.frame(listNames),
#           row.names=FALSE,
#           col.names = FALSE, 
#           file = file1,
#           sep=",")

# mydat1 <- read.csv(file2)
# mydat2 <- read.csv(file1)
# fulldat <- cbind(mydat2[1],mydat1[1])
# fulltab <- as.data.table(fulldat)
# # name table columns
# names(fulltab)[1] <- paste("UMI")
# names(fulltab)[2] <- paste("cluster")
# # group data based on clusters
# group_by(fulltab, cluster)
# # create a table counting unqiue UMIs/cells per cluster
# tabPerClus <- fulltab %>% group_by(cluster) %>% count()
# type <- sub("\\_.*","",fulltab$UMI)
# fulltab <- cbind(fulltab, type)

# divout <-  paste(filedir,"CellBreakdown_PerClusterPerType_res-0.6_k150",date,".csv",sep="") 
# A <- fulltab %>% group_by(cluster) %>% count(type)
# write.csv(A, file=divout)
# rm(mydat1,mydat2,fulldat,fulltab,tabPerClus,divout,A)

########## OLD ANALYSIS 

# OPTIONAL STEP - JACK-STRAW ANALYSIS
# Note, previously on 9/5 didn't complete because it errored out due to meeting the global.maxSize limit, 
# so increased for this iteration. Checked compute node 1-20 has 528GB, so can increase to 512 GB limit
#FullSet <- JackStraw(object=FullSet, num.replicate=500, dims=200)
#FullSet <- ScoreJackStraw(object=FullSet, dims=1:200)

#pdf("./analysis/Healthy-4pts+CMML-16pts-Cohort_CD34only_JackStrawPlot_100dim_2019-09-20.pdf")
#jsPlot <- JackStrawPlot(object=FullSet, dims = 1:100)
#print(jsPlot)
#dev.off()

#pdf("./analysis/Healthy-4pts+CMML-16pts-Cohort_CD34only_JackStrawPlot_nolegened_100dim_2019-09-20.pdf")
#jsPlot <- JackStrawPlot(object=FullSet, dims = 1:100)+theme(legend.position="none")
#print(jsPlot)
#dev.off()

#pdf("./analysis/Healthy-4pts+CMML-16pts-Cohort_CD34only_JackStraw_ElbowPlot_100dim_2019-09-20.pdf") 
#ePlot <- ElbowPlot(object=FullSet, dims= 1:100)
#print(ePlot)
#dev.off()

#saveRDS(FullSet, "./analysis/Healthy-4pts+CMML-16pts-Cohort_CD34only_SeuratObject_thruStep7+JS_20190920.rds")

###### Step 8: Visualization with UMAP
# FullSet <- RunUMAP(object=FullSet, reduction = "pca", dims = 1:100, n.neighbors=175)
# # Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
# # To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
# # This message will be shown once per session
# # 19:06:14 UMAP embedding parameters a = 0.9922 b = 1.112
# # 19:06:15 Read 174714 rows and found 100 numeric columns
# # 19:06:15 Using Annoy for neighbor search, n_neighbors = 30
# # 19:06:15 Building Annoy index with metric = cosine, n_trees = 50
# # 0%   10   20   30   40   50   60   70   80   90   100%
# # [----|----|----|----|----|----|----|----|----|----|
# # **************************************************|
# # 19:07:38 Writing NN index file to temp file /scratch/2490682.colossus.local/RtmpjCAOF2/file32948023324
# # 19:07:38 Searching Annoy index using 16 threads, search_k = 3000
# # 19:08:07 Annoy recall = 100%
# # 19:08:08 Commencing smooth kNN distance calibration using 16 threads
# # 19:08:21 Initializing from normalized Laplacian + noise
# # 19:08:49 Commencing optimization for 200 epochs, with 8082714 positive edges
# # 0%   10   20   30   40   50   60   70   80   90   100%
# # [----|----|----|----|----|----|----|----|----|----|
# # **************************************************|
# # 19:12:42 Optimization finished

# FullSet <- FindNeighbors(object=FullSet, dims = 1:50)
# FullSet <- FindClusters(object=FullSet, resolution = 0.4)

# pdf(paste(filedir,"UMAP_byClusters",date,".pdf",sep=""), width=12, height=8)
# dPlot <- DimPlot(FullSet, reduction.use="umap") + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_clustersByType",date,".pdf",sep=""), width=12, height=8)
# dPlot <- DimPlot(FullSet, reduction.use="umap", group.by='orig.ident') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_clustersByType_noLegend",date,".pdf",sep=""), width=12, height=8)
# dPlot <- DimPlot(FullSet, reduction.use="umap", group.by='orig.ident') + theme(legend.position = "none") + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# pdf(paste(filedir,"UMAP_byClusters_noLims",date,".pdf",sep=""), width=12, height=8)
# dPlot <- DimPlot(FullSet, reduction.use="umap")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_clustersByType_noLims",date,".pdf",sep=""), width=12, height=8)
# dPlot <- DimPlot(FullSet, reduction.use="umap", group.by='orig.ident')
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_clustersByType_noLegend_noLims",date,".pdf",sep=""), width=12, height=8)
# dPlot <- DimPlot(FullSet, reduction.use="umap", group.by='orig.ident') + theme(legend.position = "none")
# print(dPlot)
# dev.off()

# postscript(file=paste(filedir,"UMAP-bySampleID",date,"_forAdobe.eps",sep=""), width=12, height=8, horizontal = TRUE, onefile = FALSE)
# DimPlot(FullSet, reduction.use="umap", group.by='orig.ident') + xlim(-15, 15) + ylim(-15, 15)
# dev.off()

# postscript(file=paste(filedir,"UMAP-bySampleID_noLims",date,"_forAdobe.eps",sep=""), width=12, height=8, horizontal = TRUE, onefile = FALSE)
# DimPlot(FullSet, reduction.use="umap", group.by='orig.ident')
# dev.off()

# date <- "_2020-03-23"
# saveRDS(FullSet, paste(filedir,"SeuratObject_thruStep8_UMAP",date,".rds",sep=""))

# ### try with larger neighborhood
# FullSet <- RunUMAP(object=FullSet, reduction = "pca", dims = 1:100)
# # 

# pdf(paste(filedir,"UMAP_v2_byClusters",date,".pdf",sep=""), width=12, height=8)
# dPlot <- DimPlot(FullSet, reduction.use="umap") + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_v2_clustersByType",date,".pdf",sep=""), width=12, height=8)
# dPlot <- DimPlot(FullSet, reduction.use="umap", group.by='orig.ident') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_v2_clustersByType_noLegend",date,".pdf",sep=""), width=12, height=8)
# dPlot <- DimPlot(FullSet, reduction.use="umap", group.by='orig.ident') + theme(legend.position = "none") + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# pdf(paste(filedir,"UMAP_v2_byClusters_noLims",date,".pdf",sep=""), width=12, height=8)
# dPlot <- DimPlot(FullSet, reduction.use="umap")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_v2_clustersByType_noLims",date,".pdf",sep=""), width=12, height=8)
# dPlot <- DimPlot(FullSet, reduction.use="umap", group.by='orig.ident')
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_v2_clustersByType_noLegend_noLims",date,".pdf",sep=""), width=12, height=8)
# dPlot <- DimPlot(FullSet, reduction.use="umap", group.by='orig.ident') + theme(legend.position = "none")
# print(dPlot)
# dev.off()

# postscript(file=paste(filedir,"UMAP-v2-bySampleID",date,"_forAdobe.eps",sep=""), width=12, height=8, horizontal = TRUE, onefile = FALSE)
# DimPlot(FullSet, reduction.use="umap", group.by='orig.ident') + xlim(-15, 15) + ylim(-15, 15)
# dev.off()

# postscript(file=paste(filedir,"UMAP-v2-bySampleID_noLims",date,"_forAdobe.eps",sep=""), width=12, height=8, horizontal = TRUE, onefile = FALSE)
# DimPlot(FullSet, reduction.use="umap", group.by='orig.ident')
# dev.off()

# saveRDS(FullSet, paste(filedir,"SeuratObject_thruStep8_UMAP_v2",date,".rds",sep=""))

# ### can't really tell a difference in UMAP, so going back to original one


# rm(FullSet)
# date <- "_2020-03-30"
# FullSet <- readRDS(paste(filedir,"SeuratObject_thruStep8_UMAP_2020-03-23.rds",sep=""))
# # to confirm have the right resolution of clusters before differential gene expression
# FullSet <- FindClusters(object=FullSet, resolution = 0.05)

# # 0%   10   20   30   40   50   60   70   80   90   100%
# # [----|----|----|----|----|----|----|----|----|----|
# # **************************************************|
# # Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck

# # Number of nodes: 174714
# # Number of edges: 6515448

# # Running Louvain algorithm...
# # Maximum modularity in 10 random starts: 0.9853
# # Number of communities: 18
# # Elapsed time: 57 seconds
# # 1 singletons identified. 17 final clusters.

# pdf(paste(filedir,"UMAP_res0.05_defaultK_byClusters",date,".pdf",sep=""), width=12, height=8)
# dPlot <- DimPlot(FullSet, reduction.use="umap") + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_clustersByType",date,".pdf",sep=""), width=12, height=8)
# dPlot <- DimPlot(FullSet, reduction.use="umap", group.by='orig.ident') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_clustersByType_noLegend",date,".pdf",sep=""), width=12, height=8)
# dPlot <- DimPlot(FullSet, reduction.use="umap", group.by='orig.ident') + theme(legend.position = "none") + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# pdf(paste(filedir,"UMAP_res0.05_defaultK_byClusters_noLims",date,".pdf",sep=""), width=12, height=8)
# dPlot <- DimPlot(FullSet, reduction.use="umap")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_clustersByType_noLims",date,".pdf",sep=""), width=12, height=8)
# dPlot <- DimPlot(FullSet, reduction.use="umap", group.by='orig.ident')
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_clustersByType_noLegend_noLims",date,".pdf",sep=""), width=12, height=8)
# dPlot <- DimPlot(FullSet, reduction.use="umap", group.by='orig.ident') + theme(legend.position = "none")
# print(dPlot)
# dev.off()

# postscript(file=paste(filedir,"UMAP-res0.05_defaultK-bySampleID",date,"_forAdobe.eps",sep=""), width=12, height=8, horizontal = TRUE, onefile = FALSE)
# DimPlot(FullSet, reduction.use="umap", group.by='orig.ident') + xlim(-15, 15) + ylim(-15, 15)
# dev.off()

# postscript(file=paste(filedir,"UMAP-res0.05_defaultK-bySampleID_noLims",date,"_forAdobe.eps",sep=""), width=12, height=8, horizontal = TRUE, onefile = FALSE)
# DimPlot(FullSet, reduction.use="umap", group.by='orig.ident')
# dev.off()


### May also want to try Leiden clustering as well
# https://rdrr.io/cran/Seurat/man/FindClusters.html
# https://github.com/satijalab/seurat/pull/1051
# https://github.com/satijalab/seurat/pull/1060/files/86989d147164e44dfb82e0e893418e0ecb1003c3

###########################################################################################################
################# PLOTS TO BE MADE AFTER ANALYSES WAS RUN #################################################
###########################################################################################################

########################################################################################
########### MARKER ANALYSIS: Proteins of Interest
########################################################################################

# ########### EXPRESSION ANALYSIS FOR PETER
# # known mutations, expression across single cells
# mutMarkers <- c("NRAS","TET2","EZH2","ASXL1","CBL","DNMT3A","GATA2","SRSF2","ETV6","BCORL1","KIT","IDH2","KRAS","RUNX1")
# pdf(paste(filedir,"mutationGeneExp_DotPlot",date,".pdf",sep=""))
# dPlot <- DotPlot(FullSet, features = rev(mutMarkers), cols = c("blue", "red"), dot.scale = 8, group.by='orig.ident') + RotatedAxis()
# print(dPlot)
# dev.off()

# CMMLonly <- subset(x=FullSet, subset = (orig.ident=="2V001" | orig.ident=="4J003" |orig.ident=="4K001" |orig.ident=="4Q001" | orig.ident=="4S001" | orig.ident=="5E001" |orig.ident=="5H001" | orig.ident=="6AC001" | orig.ident=="6AD001" | orig.ident=="6AE001" | orig.ident=="LTB3966" | orig.ident=="LTB4121" | orig.ident=="6169" | orig.ident=="SF100109101914" | orig.ident=="SF100109106293" | orig.ident=="SF100109110236" | orig.ident=="SF100109111451" | orig.ident=="SF12062800475" | orig.ident=="SF13061200056" | orig.ident=="SF14040100158" | orig.ident=="SF14060200025" | orig.ident=="SF14072200012" | orig.ident=="SF14101000049" | orig.ident=="SF16112900158" | orig.ident=="SF12042500035" | orig.ident=="SF12092600014" | orig.ident=="SF14031800065" | orig.ident=="SF14050700419" | orig.ident=="SF16026800045" | 
# orig.ident=="SF16072200003" | orig.ident=="SF16112300029"))

# pdf(paste(filedir,"mutationGeneExp_CMMLonly_DotPlot",date,".pdf",sep=""))
# dPlot <- DotPlot(CMMLonly, features = rev(mutMarkers), cols = c("blue", "red"), dot.scale = 8, group.by='orig.ident') + RotatedAxis()
# print(dPlot)
# dev.off()

# # hematopoietic markers from PALANTIR
# hematopMarkers <- c("CD34","TAL1","GATA2","GATA1","KLF1","SPI1","EGR1","IRF8","CEBPA","ELANE","JUN","TCF7","EBF1","MALAT1","CREB1")
# pdf(paste(filedir,"hematopoeiticMarkers_DotPlot",date,".pdf",sep=""))
# dPlot <- DotPlot(FullSet, features = rev(hematopMarkers), cols = c("blue", "red"), dot.scale = 8, group.by='orig.ident') + RotatedAxis()
# print(dPlot)
# dev.off()

# # Fuco Enzymes 
# markers.to.plot <- c("FUK","GMDS","TSTA3", "FUT1","FUT2","FUT3","FUT4","FUT5","FUT6","FUT7","FUT8","FUT9","FUT10","FUT11","POFUT1","POFUT2")
# pdf(paste(filedir,"fucoEnzy_DotPlot",date,".pdf",sep=""))
# dPlot <- DotPlot(FullSet, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, group.by='orig.ident') + RotatedAxis()
# print(dPlot)
# dev.off()

# # Abhishek's cytokine receptor panels
# panelA.markers.to.plot <- c("HAVCR2","CXCR1","CXCR2","CSF1R","CCR2","CXCR4","IFNGR1","IL3RA","TLR4","KIT","IL15RA","IL2RG","IL6R","IL5RA","CSF3R","TLR2","FLT3","MPL","TNFRSF1B","IL18R1","TNFRSF1A","CSF2RA")
# pdf(paste(filedir,"cytoRecept_PanelA_DotPlot",date,".pdf",sep=""))
# dPlot <- DotPlot(FullSet, features = rev(panelA.markers.to.plot), cols = c("blue", "red"), dot.scale = 8, group.by='orig.ident') + RotatedAxis()
# print(dPlot)
# dev.off()

# panelB.markers.to.plot <- c("F11R","LAIR1","TFRC", "CD33","CD36","CD180","PLAUR","CD244") ### LILRA not recognized
# pdf(paste(filedir,"cytoRecept_PanelB_DotPlot",date,".pdf",sep=""))
# dPlot <- DotPlot(FullSet, features = rev(panelB.markers.to.plot), cols = c("blue", "red"), dot.scale = 8, group.by='orig.ident') + RotatedAxis()
# print(dPlot)
# dev.off()


# # Cysteine Cathepsins
# cts.markers.to.plot <- c("CTSK","CTSS","CTSL","CTSV")
# pdf(paste(filedir,"CTS_DotPlot",date,".pdf",sep=""))
# dPlot <- DotPlot(FullSet, features = rev(cts.markers.to.plot), cols = c("blue", "red"), dot.scale = 8, group.by='orig.ident') + RotatedAxis()
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"CTS_VlnPlot_CTSS",date,".pdf",sep=""))
# plot1 <- VlnPlot(FullSet, features = "CTSS", group.by='orig.ident', pt.size = 0, combine = FALSE)
# print(plot1)
# dev.off()
# pdf(paste(filedir,"CTS_VlnPlot_CTSS_noLegend",date,".pdf",sep=""))
# plot1 <- VlnPlot(FullSet, features = "CTSS", group.by='orig.ident', pt.size = 0, combine = FALSE)+ theme(legend.position = "none")
# print(plot1)
# dev.off()
# pdf(paste(filedir,"CTS_VlnPlot_CTSL",date,".pdf",sep=""))
# plot1 <- VlnPlot(FullSet, features = "CTSL", group.by='orig.ident', pt.size = 0, combine = FALSE)
# print(plot1)
# dev.off()
# pdf(paste(filedir,"CTS_VlnPlot_CTSL_noLegend",date,".pdf",sep=""))
# plot1 <- VlnPlot(FullSet, features = "CTSL", group.by='orig.ident', pt.size = 0, combine = FALSE) 
# print(plot1)
# dev.off()


# ### PRINT OUT LOWER RESOULTION DATA CLUSTERING DATA
# # Range of resolutions: 
# # range <- c(0, 0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6)

# file1 <- "./analysis/CMMLoutputDataNames_31C+7Npts_CD34only_res0.6_20200325.csv";
# file2 <- "./analysis/CMMLoutputData_31C+7Npts_CD34only_res0.6_20200325.csv";

# cat(FullSet@meta.data$RNA_snn_res.0.6, file=file2, sep=",\n")

# listNames <- FullSet@meta.data$orig.ident
# write.table(data.frame(listNames),
#           row.names=FALSE,
#           col.names = FALSE, 
#           file = file1,
#           sep=",")

# mydat1 <- read.csv(file2)
# mydat2 <- read.csv(file1)
# fulldat <- cbind(mydat2[1],mydat1[1])
# fulltab <- as.data.table(fulldat)
# # name table columns
# names(fulltab)[1] <- paste("UMI")
# names(fulltab)[2] <- paste("cluster")
# # group data based on clusters
# group_by(fulltab, cluster)
# # create a table counting unqiue UMIs/cells per cluster
# tabPerClus <- fulltab %>% group_by(cluster) %>% count()
# type <- sub("\\_.*","",fulltab$UMI)
# fulltab <- cbind(fulltab, type)

# divout <-  paste(filedir,"CellBreakdown_PerClusterPerType_res-0.6",date,".csv",sep="") 
# A <- fulltab %>% group_by(cluster) %>% count(type)
# write.csv(A, file=divout)

# rm(file1,file2,listNames,mydat1,mydat2,fulldat,fulltab,tabPerClus,type,divout,A)

# date <- "_2020-04-01"

# file1 <- "./analysis/CMMLoutputDataNames_31C+7Npts_CD34only_res0.05_check_20200401.csv";
# file2 <- "./analysis/CMMLoutputData_31C+7Npts_CD34only_res0.05_check_20200401.csv";

# cat(FullSet@meta.data$RNA_snn_res.0.05, file=file2, sep=",\n")

# listNames <- FullSet@meta.data$orig.ident
# write.table(data.frame(listNames),
#           row.names=FALSE,
#           col.names = FALSE, 
#           file = file1,
#           sep=",")

# mydat1 <- read.csv(file2)
# mydat2 <- read.csv(file1)
# fulldat <- cbind(mydat2[1],mydat1[1])
# fulltab <- as.data.table(fulldat)
# # name table columns
# names(fulltab)[1] <- paste("UMI")
# names(fulltab)[2] <- paste("cluster")
# # group data based on clusters
# group_by(fulltab, cluster)
# # create a table counting unqiue UMIs/cells per cluster
# tabPerClus <- fulltab %>% group_by(cluster) %>% count()
# type <- sub("\\_.*","",fulltab$UMI)
# fulltab <- cbind(fulltab, type)

# divout <-  paste(filedir,"CellBreakdown_PerClusterPerType_res-0.05_check",date,".csv",sep="") 
# A <- fulltab %>% group_by(cluster) %>% count(type)
# write.csv(A, file=divout)

# data_embeddings2 <- as.data.frame(as.matrix(Embeddings(object=signatureSet2[["umap"]])))
# fwrite(x = data_embeddings2, file = "./analysis/UMAP-Embeddings_Clus5+8+9+10+15_2020-04-01.csv", row.names=TRUE, col.names=TRUE)

# data_embeddings3 <- as.data.frame(as.matrix(Embeddings(object=signatureSet[["umap"]])))
# fwrite(x = data_embeddings3, file = "./analysis/UMAP-Embeddings_Clus5+8_2020-04-01.csv", row.names=TRUE, col.names=TRUE)

# # to re-calc diversity for signatureSet
# file1 <- "./analysis/CMMLoutputDataNames_31C+7Npts_CD34only_res0.05_SigSet1_20200401.csv";
# file2 <- "./analysis/CMMLoutputData_31C+7Npts_CD34only_res0.05_SigSet1_20200401.csv";

# cat(signatureSet@meta.data$RNA_snn_res.0.05, file=file2, sep=",\n")

# listNames <- signatureSet@meta.data$orig.ident
# write.table(data.frame(listNames),
#           row.names=FALSE,
#           col.names = FALSE, 
#           file = file1,
#           sep=",")

# mydat1 <- read.csv(file2)
# mydat2 <- read.csv(file1)
# fulldat <- cbind(mydat2[1],mydat1[1])
# fulltab <- as.data.table(fulldat)
# # name table columns
# names(fulltab)[1] <- paste("UMI")
# names(fulltab)[2] <- paste("cluster")
# # group data based on clusters
# group_by(fulltab, cluster)
# # create a table counting unqiue UMIs/cells per cluster
# tabPerClus <- fulltab %>% group_by(cluster) %>% count()
# type <- sub("\\_.*","",fulltab$UMI)
# fulltab <- cbind(fulltab, type)

# divout <-  paste(filedir,"CellBreakdown_PerClusterPerType_res-0.05_SigSet1",date,".csv",sep="") 
# A <- fulltab %>% group_by(cluster) %>% count(type)
# write.csv(A, file=divout)

# # to recalc diversity for signatureSet2
# file1 <- "./analysis/CMMLoutputDataNames_31C+7Npts_CD34only_res0.05_SigSet2_20200402.csv";
# file2 <- "./analysis/CMMLoutputData_31C+7Npts_CD34only_res0.05_SigSet2_20200401.csv";

# cat(signatureSet2@meta.data$RNA_snn_res.0.05, file=file2, sep=",\n")

# listNames <- signatureSet2@meta.data$orig.ident
# write.table(data.frame(listNames),
#           row.names=FALSE,
#           col.names = FALSE, 
#           file = file1,
#           sep=",")

# mydat1 <- read.csv(file2)
# mydat2 <- read.csv(file1)
# fulldat <- cbind(mydat2[1],mydat1[1])
# fulltab <- as.data.table(fulldat)
# # name table columns
# names(fulltab)[1] <- paste("UMI")
# names(fulltab)[2] <- paste("cluster")
# # group data based on clusters
# group_by(fulltab, cluster)
# # create a table counting unqiue UMIs/cells per cluster
# tabPerClus <- fulltab %>% group_by(cluster) %>% count()
# type <- sub("\\_.*","",fulltab$UMI)
# fulltab <- cbind(fulltab, type)

# divout <-  paste(filedir,"CellBreakdown_PerClusterPerType_res-0.05_SigSet2",date,".csv",sep="") 
# A <- fulltab %>% group_by(cluster) %>% count(type)
# write.csv(A, file=divout)

## Discovered slightly different clustering with new resolution, so going back to original 0.05 resoultion used when calculating stats/groupings

# note loading data after Step 7, so will also re-run UMAP 
# FullSetMultiRes <- readRDS("./analysis/Healthy-7pts+CMML-31pts-Cohort_CD34only+33kGenes_SeuratObject_thruStep7_multiRes_2020-03-21.rds")

# levels(FullSetMultiRes@active.ident)
# # [1] "0"  "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14"
# # [16] "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29"
# # [31] "30" "31" "32" "33" "34" "35" "36" "37" "38" "39" "40" "41" "42" "43" "44"
# # [46] "45" "46" "47" "48" "49" "50" "51" "52"

# FullSetMultiRes@active.ident <- FullSetMultiRes@meta.data$RNA_snn_res.0.05

# levels(FullSetMultiRes@active.ident)
# #  [1] "0"  "1"  "10" "11" "12" "13" "14" "15" "16" "17" "2"  "3"  "4"  "5"  "6" 
# # [16] "7"  "8"  "9"
# # this seems to have worked

# # also stashed in meta data as "cluster.ident"; confirmed with levels
# FullSetMultiRes$cluster.ident <- FullSetMultiRes@active.ident 

# saveRDS(FullSetMultiRes, paste(filedir,"SeuratObject_thruStep7_res0.05-as-active",date,".rds",sep=""))


## Step 8: VISULATION WITH UMAP
# FullSetMultiRes <- RunUMAP(object=FullSetMultiRes, reduction = "pca", dims = 1:100)

# # 16:01:08 UMAP embedding parameters a = 0.9922 b = 1.112
# # 16:01:09 Read 174714 rows and found 100 numeric columns
# # 16:01:09 Using Annoy for neighbor search, n_neighbors = 30
# # 16:01:09 Building Annoy index with metric = cosine, n_trees = 50
# # 0%   10   20   30   40   50   60   70   80   90   100%
# # [----|----|----|----|----|----|----|----|----|----|
# # **************************************************|
# # 16:02:28 Writing NN index file to temp file /scratch/2495825.colossus.local/RtmpKcRFco/fileb747719a0364
# # 16:02:28 Searching Annoy index using 16 threads, search_k = 3000
# # 16:03:00 Annoy recall = 100%
# # 16:03:01 Commencing smooth kNN distance calibration using 16 threads
# # 16:03:12 Initializing from normalized Laplacian + noise
# # 16:03:38 Commencing optimization for 200 epochs, with 8082714 positive edges
# # 0%   10   20   30   40   50   60   70   80   90   100%
# # [----|----|----|----|----|----|----|----|----|----|
# # **************************************************|
# # 16:07:28 Optimization finished

# pdf(paste(filedir,"UMAP_res0.05_defaultK_byClusters",date,".pdf",sep=""), width=12, height=8)
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_clustersByType",date,".pdf",sep=""), width=12, height=8)
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap") + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_clustersByType_noLegend",date,".pdf",sep=""), width=12, height=8)
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap") + theme(legend.position = "none") + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# pdf(paste(filedir,"UMAP_res0.05_defaultK_byClusters_noLims",date,".pdf",sep=""), width=12, height=8)
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", group.by='cluster')
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_clustersByType_noLims",date,".pdf",sep=""), width=12, height=8)
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_clustersByType_noLegend_noLims",date,".pdf",sep=""), width=12, height=8)
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap") + theme(legend.position = "none")
# print(dPlot)
# dev.off()

# postscript(file=paste(filedir,"UMAP-res0.05_defaultK-bySampleID",date,"_forAdobe.eps",sep=""), width=12, height=8, horizontal = TRUE, onefile = FALSE)
# DimPlot(FullSetMultiRes, reduction.use="umap") + xlim(-15, 15) + ylim(-15, 15)
# dev.off()

# postscript(file=paste(filedir,"UMAP-res0.05_defaultK-bySampleID_noLims",date,"_forAdobe.eps",sep=""), width=12, height=8, horizontal = TRUE, onefile = FALSE)
# DimPlot(FullSetMultiRes, reduction.use="umap")
# dev.off()

# ### HIGHLIGHTED UMAP PLOTS

# ## re-defining whichcells, stash clusters into $cluster, then make $orig.ident to @active.ident
# FullSetMultiRes$cluster <- FullSetMultiRes@active.ident
# FullSetMultiRes@active.ident <- FullSetMultiRes$orig.ident

# n1 <- WhichCells(FullSetMultiRes, idents = c("CD34"))
# n2 <- WhichCells(FullSetMultiRes, idents = c("SettyPt1"))
# n3 <- WhichCells(FullSetMultiRes, idents = c("SettyPt2"))
# n4 <- WhichCells(FullSetMultiRes, idents = c("SettyPt3"))
# n5 <- WhichCells(FullSetMultiRes, idents = c("HuaPt1"))
# n6 <- WhichCells(FullSetMultiRes, idents = c("HuaPt2"))
# n7 <- WhichCells(FullSetMultiRes, idents = c("HuaPt4"))

# Normals <- c(n1,n2,n3,n4,n5,n6,n7)

# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_Normal_All-7pts",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=Normals) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()

# l1 <- WhichCells(FullSetMultiRes, idents = c("SF16072200003"))
# l2 <- WhichCells(FullSetMultiRes, idents = c("4K001"))
# l3 <- WhichCells(FullSetMultiRes, idents = c("SF16112300029"))
# l4 <- WhichCells(FullSetMultiRes, idents = c("5E001"))
# l5 <- WhichCells(FullSetMultiRes, idents = c("SF14050700419"))
# l6 <- WhichCells(FullSetMultiRes, idents = c("SF12042500035"))
# l7 <- WhichCells(FullSetMultiRes, idents = c("LTB3966"))
# l8 <- WhichCells(FullSetMultiRes, idents = c("SF100109101914"))
# l9 <- WhichCells(FullSetMultiRes, idents = c("4S001"))
# l10 <- WhichCells(FullSetMultiRes, idents = c("4Q001"))
# l11 <- WhichCells(FullSetMultiRes, idents = c("2V001"))
# l12 <- WhichCells(FullSetMultiRes, idents = c("LTB4121"))

# LowDeltaQDiversity <- c(l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12)

# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_CMML_Low-deltaqD-Diversity",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=LowDeltaQDiversity) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()

# h1 <- WhichCells(FullSetMultiRes, idents = c("SF100109111451"))
# h2 <- WhichCells(FullSetMultiRes, idents = c("LTB6169"))
# h3 <- WhichCells(FullSetMultiRes, idents = c("SF14031800065"))
# h4 <- WhichCells(FullSetMultiRes, idents = c("SF100109106293"))
# h5 <- WhichCells(FullSetMultiRes, idents = c("4J003"))
# h6 <- WhichCells(FullSetMultiRes, idents = c("6AC001"))
# h7 <- WhichCells(FullSetMultiRes, idents = c("SF14040100158"))
# h8 <- WhichCells(FullSetMultiRes, idents = c("SF16026800045"))
# h9 <- WhichCells(FullSetMultiRes, idents = c("6AE001"))
# h10 <- WhichCells(FullSetMultiRes, idents = c("SF16112900158"))
# h11 <- WhichCells(FullSetMultiRes, idents = c("5H001"))
# h12 <- WhichCells(FullSetMultiRes, idents = c("SF100109110236"))
# h13 <- WhichCells(FullSetMultiRes, idents = c("SF12092600014"))
# h14 <- WhichCells(FullSetMultiRes, idents = c("SF14101000049"))
# h15 <- WhichCells(FullSetMultiRes, idents = c("SF12062800475"))
# h16 <- WhichCells(FullSetMultiRes, idents = c("SF14072200012"))
# h17 <- WhichCells(FullSetMultiRes, idents = c("SF13061200056"))
# h18 <- WhichCells(FullSetMultiRes, idents = c("6AD001"))
# h19 <- WhichCells(FullSetMultiRes, idents = c("SF14060200025"))

# HighDeltaQDiversity <- c(h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12,h13,h14,h15,h16,h17,h18,h19)

# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_CMML_High-deltaqD-Diversity",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=HighDeltaQDiversity) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()

# pdf(paste(filedir,"UMAP_res0.05_defaultK_CMML_High-deltaqD-diversity_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=HighDeltaQDiversity, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# pdf(paste(filedir,"UMAP_res0.05_defaultK_CMML_Low-deltaqD-diversity_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=LowDeltaQDiversity, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# pdf(paste(filedir,"UMAP_res0.05_defaultK_Normal_All-7-Pts_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=Normals, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# ## NORMAL
# # ZhengCells <- n1
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_Normal_Zheng-Pt",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=n1) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Normal_Zheng-Pt_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=n1, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# SettyPt1Cells <- n2
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_Normal_Setty-Pt1",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=SettyPt1Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Normal_Setty-Pt1_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=SettyPt1Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# SettyPt2Cells <- n3
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_Normal_Setty-Pt2",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=SettyPt2Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Normal_Setty-Pt2_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=SettyPt2Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# SettyPt3Cells <- n4
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_Normal_Setty-Pt3",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=SettyPt3Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Normal_Setty-Pt3_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=SettyPt3Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# HuaPt1Cells <- n5
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_Normal_Hua-Pt1",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=HuaPt1Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Normal_Hua-Pt1_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=HuaPt1Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# HuaPt2Cells <- n6
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_Normal_Hua-Pt2",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=HuaPt2Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Normal_Hua-Pt2_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=HuaPt2Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# HuaPt4Cells <- n7
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_Normal_Hua-Pt4",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=HuaPt4Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Normal_Hua-Pt4_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=HuaPt4Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# ## CMML

# Pt2V001Cells <- l11
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_CMML_2-V-001",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=Pt2V001Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_CMML_2-V-001_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=Pt2V001Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# Pt4J003Cells <- h5
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_CMML_4-J-003",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=Pt4J003Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_CMML_4-J-003_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=Pt4J003Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# Pt4K001Cells <- l2
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_CMML_4-K-001",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=Pt4K001Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_CMML_4-K-001_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=Pt4K001Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# Pt4Q001Cells <- l10
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_CMML_4-Q-001",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=Pt4Q001Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_CMML_4-Q-001_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=Pt4Q001Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# Pt4S001Cells <- l9
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_CMML_4-S-001",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=Pt4S001Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_CMML_4-S-001_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=Pt4S001Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# Pt5E001Cells <- l4
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_CMML_5-E-001",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=Pt5E001Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_CMML_5-E-001_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=Pt5E001Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# Pt5H001Cells <- h11
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_CMML_5-H-001",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=Pt5H001Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_CMML_5-H-001_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=Pt5H001Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# Pt6AC001Cells <- h6
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_CMML_6-AC-001",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=Pt6AC001Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_CMML_6-AC-001_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=Pt6AC001Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# Pt6AD001Cells <- h18
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_CMML_6-AD-001",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=Pt6AD001Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_CMML_6-AD-001_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=Pt6AD001Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# Pt6AE001Cells <- h9
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_CMML_6-AE-001",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=Pt6AE001Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_CMML_6-AE-001_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=Pt6AE001Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# LTB3966Cells <- l7
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_CMML_LTB3966",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=LTB3966Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_CMML_LTB3966_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=LTB3966Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# LTB4121Cells <- l12
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_CMML_LTB4121",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=LTB4121Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_CMML_LTB4121_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=LTB4121Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# LTB5109Cells <- h2
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_CMML_LTB5109",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=LTB5109Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_CMML_LTB5109_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=LTB5109Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# SF100109101914Cells <- l8
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_CMML_SF-100109-101914",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=SF100109101914Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_CMML_SF-100109-101914_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=SF100109101914Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# SF100109106293Cells <- h4
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_CMML_SF-100109-106293",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=SF100109106293Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_CMML_SF-100109-106293_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=SF100109106293Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# SF100109110236Cells <- h12
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_CMML_SF-100109-110236",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=SF100109110236Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_CMML_SF-100109-110236_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=SF100109110236Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# SF100109111451Cells <- h1
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_CMML_SF-100109-111451",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=SF100109111451Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_CMML_SF-100109-111451_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=SF100109111451Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# SF12042500035Cells <- l6
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_CMML_SF-120425-00035",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=SF12042500035Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_CMML_SF-120425-00035_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=SF12042500035Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# SF12062800475Cells <- h15
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_CMML_SF-120628-00475",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=SF12062800475Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_CMML_SF-120628-00475_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=SF12062800475Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# SF12092600014Cells <- h13
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_CMML_SF-120926-00014",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=SF12092600014Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_CMML_SF-120926-00014_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=SF12092600014Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# SF13061200056Cells <- h17
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_CMML_SF-130612-00056",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=SF13061200056Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_CMML_SF-130612-00056_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=SF13061200056Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# SF14031800065Cells <- h3
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_CMML_SF-140318-00065",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=SF14031800065Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_CMML_SF-140318-00065_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=SF14031800065Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# SF14040100158Cells <- h7
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_CMML_SF-140401-00158",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=SF14040100158Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_CMML_SF-140401-00158_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=SF14040100158Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# SF14050700419Cells <- l5
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_CMML_SF-140507-00419",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=SF14040100158Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_CMML_SF-140507-00419_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=SF14050700419Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# SF14060200025Cells <- h19
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_CMML_SF-140602-00025",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=SF14060200025Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_CMML_SF-140602-00025_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=SF14060200025Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# SF14072200012Cells <- h16
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_CMML_SF-140722-00012",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=SF14072200012Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_CMML_SF-140722-00012_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=SF14072200012Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# SF14101000049Cells <- h14
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_CMML_SF-141010-00049",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=SF14101000049Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_CMML_SF-141010-00049_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=SF14101000049Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# SF16026800045Cells <- h8
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_CMML_SF-160268-00045",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=SF16026800045Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_CMML_SF-160268-00045_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=SF16026800045Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# SF16072200003Cells <- l1
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_CMML_SF-160722-00003",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=SF16072200003Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_CMML_SF-160722-00003_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=SF16072200003Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# SF16112300029Cells <- l3
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_CMML_SF-161123-00029",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=SF16112300029Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_CMML_SF-161123-00029_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=SF16112300029Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()

# SF16112900158Cells <- h10
# pdf(paste(filedir,"UMAP_res0.05_defaultK_Highlighted_CMML_SF-161129-00158",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction.use="umap", cells.highlight=SF16112900158Cells) + xlim(-15, 15) + ylim(-15, 15)+ theme(legend.position = "none")
# print(dPlot)
# dev.off()
# pdf(paste(filedir,"UMAP_res0.05_defaultK_CMML_SF-161129-00158_Only_byClusters",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, cells=SF16112900158Cells, reduction.use="umap", group.by='cluster') + xlim(-15, 15) + ylim(-15, 15)
# print(dPlot)
# dev.off()


##### PICKING UP 7/22/2020

#### FOR PLOTTING CELL CYCLE MARKERS
## Seurat tutorial: https://satijalab.org/seurat/v3.1/cell_cycle_vignette.html 
## Cell Cycle Genes seem to be loaded with Seurat
s.genes <- cc.genes$s.genes
#  [1] "MCM5"     "PCNA"     "TYMS"     "FEN1"     "MCM2"     "MCM4"    
#  [7] "RRM1"     "UNG"      "GINS2"    "MCM6"     "CDCA7"    "DTL"     
# [13] "PRIM1"    "UHRF1"    "MLF1IP"   "HELLS"    "RFC2"     "RPA2"    
# [19] "NASP"     "RAD51AP1" "GMNN"     "WDR76"    "SLBP"     "CCNE2"   
# [25] "UBR7"     "POLD3"    "MSH2"     "ATAD2"    "RAD51"    "RRM2"    
# [31] "CDC45"    "CDC6"     "EXO1"     "TIPIN"    "DSCC1"    "BLM"     
# [37] "CASP8AP2" "USP1"     "CLSPN"    "POLA1"    "CHAF1B"   "BRIP1"   
# [43] "E2F8"    
g2m.genes <- cc.genes$g2m.genes
#  [1] "HMGB2"   "CDK1"    "NUSAP1"  "UBE2C"   "BIRC5"   "TPX2"    "TOP2A"  
#  [8] "NDC80"   "CKS2"    "NUF2"    "CKS1B"   "MKI67"   "TMPO"    "CENPF"  
# [15] "TACC3"   "FAM64A"  "SMC4"    "CCNB2"   "CKAP2L"  "CKAP2"   "AURKB"  
# [22] "BUB1"    "KIF11"   "ANP32E"  "TUBB4B"  "GTSE1"   "KIF20B"  "HJURP"  
# [29] "CDCA3"   "HN1"     "CDC20"   "TTK"     "CDC25C"  "KIF2C"   "RANGAP1"
# [36] "NCAPD2"  "DLGAP5"  "CDCA2"   "CDCA8"   "ECT2"    "KIF23"   "HMMR"   
# [43] "AURKA"   "PSRC1"   "ANLN"    "LBR"     "CKAP5"   "CENPE"   "CTCF"   
# [50] "NEK2"    "G2E3"    "GAS2L3"  "CBX5"    "CENPA" 


FullSetMultiRes_CC <- CellCycleScoring(FullSetMultiRes, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# Warning: The following features are not present in the object: MLF1IP, not searching for symbol synonyms

# head(FullSetMultiRes_CC[[]])

pdf(paste(filedir,"RidgePlot_CellCycle_UMAP_res0.05_defaultK",date,".pdf",sep=""))
dPlot <- RidgePlot(FullSetMultiRes_CC, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
print(dPlot)
dev.off()

FullSetMultiRes_CC <- RunPCA(FullSetMultiRes_CC, features = c(s.genes, g2m.genes), dim=1:100)
# Warning in irlba(A = t(x = object), nv = npcs, ...) :
#   You're computing too large a percentage of total singular values, use a standard svd instead.
# PC_ 1 
# Positive:  MKI67, TOP2A, NUSAP1, UBE2C, BIRC5, AURKB, CENPF, TPX2, GTSE1, HMGB2 
#       TYMS, HMMR, CDK1, CCNB2, CDC20, CDCA3, CDCA8, CENPE, CENPA, CKS1B 
#       DLGAP5, NUF2, CKAP2L, RRM2, NDC80, KIF11, TUBB4B, SMC4, KIF23, KIF2C 
# Negative:  TIPIN, UNG, G2E3, POLA1, MSH2, GAS2L3, CASP8AP2, RPA2, UBR7, CHAF1B 
#       POLD3, BRIP1, MCM6, CDCA7, LBR, RFC2, MCM2, HN1, EXO1, PRIM1 
#       CCNE2, RANGAP1, BLM, MCM5, SLBP, CTCF, DTL, WDR76, ECT2, FAM64A 
# PC_ 2 
# Positive:  PCNA, GINS2, MCM5, MCM4, HELLS, TYMS, UNG, MCM6, FEN1, CLSPN 
#       MCM2, CDC45, CDC6, NASP, SLBP, UHRF1, DTL, CDCA7, RAD51, WDR76 
#       GMNN, RFC2, ATAD2, CHAF1B, USP1, RRM1, RAD51AP1, RRM2, CCNE2, DSCC1 
# Negative:  CENPA, CDC20, HMMR, DLGAP5, CENPE, UBE2C, AURKA, CCNB2, GTSE1, CDCA3 
#       TOP2A, TPX2, PSRC1, CENPF, FAM64A, CDCA8, NEK2, KIF23, AURKB, BUB1 
#       CDCA2, HN1, KIF2C, BIRC5, NUF2, CKAP5, TUBB4B, MKI67, NCAPD2, GAS2L3 
# PC_ 3 
# Positive:  RRM2, E2F8, ATAD2, AURKB, HJURP, MKI67, KIF11, NUSAP1, NDC80, TACC3 
#       CCNE2, RAD51AP1, ANLN, CLSPN, CKAP2L, TMPO, UHRF1, CDK1, CTCF, SMC4 
#       BLM, BRIP1, USP1, CDCA2, WDR76, HELLS, CDCA3, KIF23, DSCC1, GTSE1 
# Negative:  UNG, CDC20, CCNB2, RPA2, MCM6, HN1, MSH2, MCM2, MCM5, GINS2 
#       RANGAP1, HMMR, DLGAP5, AURKA, NASP, CKS2, CENPA, MCM4, CENPF, BUB1 
#       FAM64A, CENPE, CKAP5, NEK2, KIF20B, LBR, PRIM1, RFC2, TIPIN, CKS1B 
# PC_ 4 
# Positive:  HN1, CDK1, TUBB4B, CDC45, CDC6, CKS1B, RAD51, EXO1, FEN1, GMNN 
#       RRM2, RAD51AP1, DSCC1, KIF2C, RPA2, TTK, HJURP, UNG, AURKA, KIF23 
#       UBE2C, RFC2, CKAP2L, NDC80, E2F8, DTL, CDCA3, PCNA, ANLN, CDC25C 
# Negative:  CDCA7, CBX5, CTCF, LBR, ANP32E, UHRF1, USP1, SMC4, NASP, HELLS 
#       NCAPD2, TMPO, HMGB2, CKAP5, CENPF, CASP8AP2, MKI67, CKAP2, G2E3, MCM2 
#       ATAD2, GTSE1, PRIM1, CENPE, MCM4, CENPA, POLA1, HMMR, NUSAP1, MSH2 
# PC_ 5 
# Positive:  LBR, ANP32E, RPA2, SLBP, HN1, CTCF, RANGAP1, BIRC5, TUBB4B, RRM2 
#       UBR7, CCNB2, CCNE2, TIPIN, HMGB2, CDC20, NCAPD2, MKI67, RRM1, CENPF 
#       UHRF1, RFC2, FEN1, SMC4, CDC6, E2F8, GMNN, MSH2, TMPO, RAD51 
# Negative:  CDCA7, PRIM1, HELLS, G2E3, PSRC1, AURKA, CBX5, CKAP2, MCM2, BRIP1 
#       BLM, TTK, ANLN, KIF23, CDCA2, UNG, MCM5, CKS2, NDC80, MCM4 
#       CDC25C, GAS2L3, CKS1B, MCM6, NEK2, RAD51AP1, UBE2C, CENPA, CDK1, KIF2C 
# Warning message:
# In PrepDR(object = object, features = features, verbose = verbose) :
#   The following 1 features requested have not been scaled (running reduction without them): MLF1IP


## regressing out cell cycle score differences during data scaling -- alternative to minimize impact on downstream 
FullSetMultiRes_CC$CC.Difference <- FullSetMultiRes_CC$S.Score - FullSetMultiRes_CC$G2M.Score

FullSetMultiRes_CC <- ScaleData(FullSetMultiRes_CC, vars.to.regress = "CC.Difference", features = rownames(FullSetMultiRes_CC))
# Regressing out CC.Difference
# Centering and scaling data matrix

FullSetMultiRes_CC <- RunPCA(FullSetMultiRes_CC, features = c(s.genes, g2m.genes), verbose = TRUE, npcs=200)
# Warning in irlba(A = t(x = object), nv = npcs, ...) :
#   You're computing too large a percentage of total singular values, use a standard svd instead.

# Warning in irlba(A = t(x = object), nv = npcs, ...) :
#   did not converge--results might be invalid!; try increasing work or maxit
# PC_ 1 
# Positive:  TYMS, TOP2A, MKI67, UBE2C, NUSAP1, PCNA, BIRC5, CENPF, AURKB, TPX2 
#       CKS1B, CLSPN, RRM2, FEN1, CCNB2, CDK1, GTSE1, CDC20, HMMR, NASP 
#       GINS2, HMGB2, RAD51AP1, CDCA3, TUBB4B, GMNN, CDCA8, CKAP2L, MCM4, MCM5 
# Negative:  G2E3, GAS2L3, LBR, POLA1, PSRC1, FAM64A, RANGAP1, CASP8AP2, ECT2, UBR7 
#       TIPIN, BRIP1, HN1, CCNE2, CTCF, CDC25C, POLD3, EXO1, MSH2, BLM 
#       CHAF1B, CKAP5, ANLN, RPA2, PRIM1, CKAP2, E2F8, TACC3, TTK, NCAPD2 
# PC_ 2 
# Positive:  RRM2, ATAD2, E2F8, HMGB2, UHRF1, AURKB, NUSAP1, MKI67, CCNE2, KIF11 
#       USP1, HJURP, HELLS, NDC80, SMC4, RAD51AP1, CLSPN, CKAP2L, TMPO, TACC3 
#       ANLN, GTSE1, CDK1, BIRC5, CDCA2, WDR76, PCNA, CDCA3, DSCC1, CTCF 
# Negative:  HN1, UNG, CDC20, RPA2, MCM6, CCNB2, MCM5, MSH2, RANGAP1, AURKA 
#       MCM2, HMMR, DLGAP5, GINS2, TUBB4B, NASP, CKS2, RFC2, CENPF, FAM64A 
#       MCM4, CKAP5, CENPA, BUB1, CENPE, KIF20B, CKS1B, LBR, TIPIN, NEK2 
# PC_ 3 
# Positive:  CDK1, HN1, TUBB4B, KIF23, RRM2, CDC45, HJURP, AURKA, KIF2C, UBE2C 
#       CDCA3, EXO1, RAD51, NDC80, CDC6, CKAP2L, RAD51AP1, E2F8, TTK, DSCC1 
#       CKS1B, FEN1, CDC25C, ANLN, GMNN, CDCA2, RPA2, TOP2A, CDCA8, NEK2 
# Negative:  CDCA7, LBR, CBX5, ANP32E, CTCF, UHRF1, HMGB2, SMC4, USP1, HELLS 
#       NASP, TMPO, NCAPD2, CKAP5, CENPF, MCM2, CKAP2, MCM4, CASP8AP2, G2E3 
#       PRIM1, MCM6, RANGAP1, MKI67, MSH2, CKS2, KIF20B, MCM5, BIRC5, GINS2 
# PC_ 4 
# Positive:  HN1, LBR, TUBB4B, CTCF, ANP32E, RANGAP1, TACC3, TMPO, CDC6, E2F8 
#       RRM2, FEN1, CDC45, RAD51, RFC2, RPA2, CBX5, EXO1, CKS1B, DTL 
#       CLSPN, GMNN, RRM1, PCNA, KIF20B, POLD3, WDR76, MCM6, RAD51AP1, SMC4 
# Negative:  CENPA, HMMR, PSRC1, AURKA, CDCA7, CDC20, CENPE, DLGAP5, GTSE1, NEK2 
#       BUB1, PRIM1, UHRF1, CDCA2, CDCA3, UBE2C, FAM64A, CDCA8, USP1, NUF2 
#       HELLS, GAS2L3, CCNB2, TPX2, KIF23, TOP2A, AURKB, KIF2C, CENPF, CDC25C 
# PC_ 5 
# Positive:  SLBP, HMGB2, CDC20, CCNB2, BIRC5, CCNE2, GMNN, PCNA, CKS2, UBR7 
#       CDC6, CKAP2, TYMS, CHAF1B, ANP32E, DSCC1, DLGAP5, CKS1B, UHRF1, PRIM1 
#       NUF2, FEN1, HMMR, DTL, SMC4, TIPIN, CDC45, RPA2, CENPF, USP1 
# Negative:  CASP8AP2, POLD3, G2E3, BRIP1, CKAP5, PSRC1, BLM, KIF23, CDCA2, CDCA7 
#       HELLS, AURKA, TACC3, CTCF, MCM6, CBX5, GTSE1, ANLN, ECT2, POLA1 
#       TMPO, CDC25C, NCAPD2, CDCA3, GAS2L3, KIF11, E2F8, CENPE, NDC80, CLSPN 
# Warning message:
# In PrepDR(object = object, features = features, verbose = verbose) :
#   The following 1 features requested have not been scaled (running reduction without them): MLF1IP

pdf(paste(filedir,"DimPlot_CellCycle_newPCA_afterRmDiff_res0.05_defaultK",date,".pdf",sep=""))
d1Plot <- DimPlot(FullSetMultiRes_CC, reduction="pca", group.by='Phase')  
print(d1Plot)
dev.off()

##### UPDATE Step 8: Visualization with UMAP
FullSetMultiRes_CC <- RunUMAP(object=FullSetMultiRes_CC, reduction = "pca", dims = 1:50)

#FullSetMultiRes_CC <- FindNeighbors(object=FullSetMultiRes_CC, dims = 1:50)
#FullSetMultiRes_CC <- FindClusters(object=FullSetMultiRes_CC, resolution = 0.05, n.start=100)

pdf(paste(filedir,"updatedUMAP_byClusters",date,".pdf",sep=""), width=12, height=8)
dPlot <- DimPlot(FullSetMultiRes_CC, reduction.use="umap") + xlim(-15, 15) + ylim(-15, 15)
print(dPlot)
dev.off()
pdf(paste(filedir,"updatedUMAP_clustersByType",date,".pdf",sep=""), width=12, height=8)
dPlot <- DimPlot(FullSetMultiRes_CC, reduction.use="umap", group.by='orig.ident') + xlim(-15, 15) + ylim(-15, 15)
print(dPlot)
dev.off()
pdf(paste(filedir,"updatedUMAP_clustersByType_noLegend",date,".pdf",sep=""), width=12, height=8)
dPlot <- DimPlot(FullSetMultiRes_CC, reduction.use="umap", group.by='orig.ident') + theme(legend.position = "none") + xlim(-15, 15) + ylim(-15, 15)
print(dPlot)
dev.off()

pdf(paste(filedir,"updatedUMAP_byClusters_noLims",date,".pdf",sep=""), width=12, height=8)
dPlot <- DimPlot(FullSetMultiRes_CC, reduction.use="umap")
print(dPlot)
dev.off()
pdf(paste(filedir,"updatedUMAP_clustersByType_noLims",date,".pdf",sep=""), width=12, height=8)
dPlot <- DimPlot(FullSetMultiRes_CC, reduction.use="umap", group.by='orig.ident')
print(dPlot)
dev.off()
pdf(paste(filedir,"updatedUMAP_clustersByType_noLegend_noLims",date,".pdf",sep=""), width=12, height=8)
dPlot <- DimPlot(FullSetMultiRes_CC, reduction.use="umap", group.by='orig.ident') + theme(legend.position = "none")
print(dPlot)
dev.off()

file1 <- "./analysis/rm-cell-cycle-effect/CMMLoutputDataNames_31C+7Npts_CD34only_res0.05_rmCC_20200722.csv";
file2 <- "./analysis/rm-cell-cycle-effect/CMMLoutputData_31C+7Npts_CD34only_res0.05_rmCC_20200722.csv";

cat(FullSetMultiRes_CC@meta.data$RNA_snn_res.0.05, file=file2, sep=",\n")

listNames <- FullSetMultiRes_CC@meta.data$orig.ident
write.table(data.frame(listNames),
          row.names=FALSE,
          col.names = FALSE, 
          file = file1,
          sep=",")

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
tabPerClus <- fulltab %>% group_by(cluster) %>% count()
type <- sub("\\_.*","",fulltab$UMI)
fulltab <- cbind(fulltab, type)

divout <-  paste(filedir,"CellBreakdown_PerClusterPerType_res-0.05_rmCC",date,".csv",sep="") 
A <- fulltab %>% group_by(cluster) %>% count(type)
write.csv(A, file=divout)

saveRDS(FullSetMultiRes_CC, paste(filedir,"SeuratObject_rmCC_res0.05",date,".rds",sep=""))

# ## PCA pre removal of cell cycle effect
# FullSetMultiRes <- CellCycleScoring(FullSetMultiRes, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# FullSetMultiRes <- RunPCA(FullSetMultiRes, features = c(s.genes, g2m.genes))
# FullSetMultiRes$CC.Difference <- FullSetMultiRes$S.Score - FullSetMultiRes$G2M.Score

# pdf(paste(filedir,"DimPlot_CellCycle_newPCA_preRmDiff_res0.05_defaultK",date,".pdf",sep=""))
# dPlot <- DimPlot(FullSetMultiRes, reduction="pca", group.by='Phase')  
# print(dPlot)
# dev.off()

# # saveRDS(FullSetMultiRes, paste(filedir,"SeuratObject_quant-CC-preRemoval",date,".rds",sep=""))

# # c(s.genes, g2m.genes)

# pdf(paste(filedir,"_RidgePlot_CellCycle_preRmDiff_shortList_res0.05_defaultK",date,".pdf",sep=""))
# dPlot <- RidgePlot(FullSetMultiRes, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
# print(dPlot)
# dev.off()

# pdf(paste(filedir,"_RidgePlot_CellCycle_preRmDiff_shortList2_res0.05_defaultK",date,".pdf",sep=""))
# dPlot <- RidgePlot(FullSetMultiRes, features = c("NASP", "HMGB2", "SLBP", "TUBB4B"), ncol = 2)
# print(dPlot)
# dev.off()

# pdf(paste(filedir,"_FeaturePlot_CellCycle_preRmDiff_shortList_res0.05_defaultK",date,".pdf",sep=""))
# dPlot <- FeaturePlot(FullSetMultiRes, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
# print(dPlot)
# dev.off()

# pdf(paste(filedir,"_FeaturePlot_CellCycle_preRmDiff_shortList2_res0.05_defaultK",date,".pdf",sep=""))
# dPlot <- FeaturePlot(FullSetMultiRes, features = c("NASP", "HMGB2", "SLBP", "TUBB4B"), ncol = 2)
# print(dPlot)
# dev.off()

# pdf(paste(filedir,"_VlnPlot_CellCycle_preRmDiff_shortList_res0.05_defaultK",date,".pdf",sep=""))
# dPlot <- VlnPlot(FullSetMultiRes, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2, group.by='Phase', pt.size=0)
# print(dPlot)
# dev.off()

# pdf(paste(filedir,"_VlnPlot_CellCycle_preRmDiff_shortList2_res0.05_defaultK",date,".pdf",sep=""))
# dPlot <- VlnPlot(FullSetMultiRes, features = c("NASP", "HMGB2", "SLBP", "TUBB4B"), ncol = 2, group.by='Phase', pt.size=0)
# print(dPlot)
# dev.off()

# pdf(paste(filedir,"_RidgePlot_CellCycle_preRmDiff_S-genes_res0.05_defaultK",date,".pdf",sep=""), width=25, height=100)
# dPlot <- RidgePlot(FullSetMultiRes, features = s.genes, ncol = 3)
# print(dPlot)
# dev.off()

# pdf(paste(filedir,"_VlnPlot_CellCycle_preRmDiff_S-genes_res0.05_defaultK",date,".pdf",sep=""), width=25, height=100)
# dPlot <- VlnPlot(FullSetMultiRes, features = s.genes, ncol = 2, group.by='Phase', pt.size=0)
# print(dPlot)
# dev.off()

# pdf(paste(filedir,"_RidgePlot_CellCycle_preRmDiff_G2M-genes_res0.05_defaultK",date,".pdf",sep=""), width=25, height=100)
# dPlot <- RidgePlot(FullSetMultiRes, features = g2m.genes, ncol = 3)
# print(dPlot)
# dev.off()

# pdf(paste(filedir,"_VlnPlot_CellCycle_preRmDiff_G2M-genes_res0.05_defaultK",date,".pdf",sep=""), width=25, height=100)
# dPlot <- VlnPlot(FullSetMultiRes, features = g2m.genes, ncol = 2, group.by='Phase', pt.size=0)
# print(dPlot)
# dev.off()

# pdf(paste(filedir,"_RidgePlot_CellCycle_preRmDiff_S-genes-short_res0.05_defaultK",date,".pdf",sep=""), width=25, height=100)
# dPlot <- RidgePlot(FullSetMultiRes, features = s.genes.short, ncol = 3)
# print(dPlot)
# dev.off()

# pdf(paste(filedir,"_RidgePlot_CellCycle_preRmDiff_G2M-genes-short_res0.05_defaultK",date,".pdf",sep=""), width=25, height=100)
# dPlot <- RidgePlot(FullSetMultiRes, features = g2m.genes.short, ncol = 3)
# print(dPlot)
# dev.off()

# pdf(paste(filedir,"_VlnPlot_CellCycle_preRmDiff_S-genes-short_res0.05_defaultK",date,".pdf",sep=""), width=25, height=100)
# dPlot <- VlnPlot(FullSetMultiRes, features = s.genes.short, ncol = 3, group.by='Phase', pt.size=0)
# print(dPlot)
# dev.off()

# pdf(paste(filedir,"_VlnPlot_CellCycle_preRmDiff_G2M-genes-short_res0.05_defaultK",date,".pdf",sep=""), width=25, height=100)
# dPlot <- VlnPlot(FullSetMultiRes, features = g2m.genes.short, ncol = 3, group.by='Phase', pt.size=0)
# print(dPlot)
# dev.off()


## PCA post removal cell cycle effect/difference
pdf(paste(filedir,"_DimPlot_CellCycle_newPCA_afterRmDiff_res0.05_defaultK",date,".pdf",sep=""))
dPlot <- DimPlot(FullSetMultiRes_CC, reduction="pca", group.by='Phase')  
print(dPlot)
dev.off()

pdf(paste(filedir,"_RidgePlot_CellCycle_afterRmDiff_shortList_res0.05_defaultK",date,".pdf",sep=""))
dPlot <- RidgePlot(FullSetMultiRes_CC, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
print(dPlot)
dev.off()

pdf(paste(filedir,"_RidgePlot_CellCycle_afterRmDiff_shortList2_res0.05_defaultK",date,".pdf",sep=""))
dPlot <- RidgePlot(FullSetMultiRes_CC, features = c("NASP", "HMGB2", "SLBP", "TUBB4B"), ncol = 2)
print(dPlot)
dev.off()

pdf(paste(filedir,"_FeaturePlot_CellCycle_afterRmDiff_shortList_res0.05_defaultK",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes_CC, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
print(dPlot)
dev.off()

pdf(paste(filedir,"_FeaturePlot_CellCycle_afterRmDiff_shortList2_res0.05_defaultK",date,".pdf",sep=""))
dPlot <- FeaturePlot(FullSetMultiRes_CC, features = c("NASP", "HMGB2", "SLBP", "TUBB4B"), ncol = 2)
print(dPlot)
dev.off()

pdf(paste(filedir,"_VlnPlot_CellCycle_afterRmDiff_shortList_res0.05_defaultK",date,".pdf",sep=""))
dPlot <- VlnPlot(FullSetMultiRes_CC, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2, group.by='Phase', pt.size=0)
print(dPlot)
dev.off()

pdf(paste(filedir,"_VlnPlot_CellCycle_afterRmDiff_shortList2_res0.05_defaultK",date,".pdf",sep=""))
dPlot <- VlnPlot(FullSetMultiRes_CC, features = c("NASP", "HMGB2", "SLBP", "TUBB4B"), ncol = 2, group.by='Phase', pt.size=0)
print(dPlot)
dev.off()

pdf(paste(filedir,"_RidgePlot_CellCycle_afterRmDiff_S-genes_res0.05_defaultK",date,".pdf",sep=""), width=25, height=100)
dPlot <- RidgePlot(FullSetMultiRes_CC, features = s.genes, ncol = 3)
print(dPlot)
dev.off()

pdf(paste(filedir,"_VlnPlot_CellCycle_afterRmDiff_S-genes_res0.05_defaultK",date,".pdf",sep=""), width=25, height=100)
dPlot <- VlnPlot(FullSetMultiRes_CC, features = s.genes, ncol = 2, group.by='Phase', pt.size=0)
print(dPlot)
dev.off()

pdf(paste(filedir,"_RidgePlot_CellCycle_afterRmDiff_G2M-genes_res0.05_defaultK",date,".pdf",sep=""), width=25, height=100)
dPlot <- RidgePlot(FullSetMultiRes_CC, features = g2m.genes, ncol = 3)
print(dPlot)
dev.off()

pdf(paste(filedir,"_VlnPlot_CellCycle_afterRmDiff_G2M-genes_res0.05_defaultK",date,".pdf",sep=""), width=25, height=100)
dPlot <- VlnPlot(FullSetMultiRes_CC, features = g2m.genes, ncol = 2, group.by='Phase', pt.size=0)
print(dPlot)
dev.off()

s.genes.short <- c("MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","GINS2","MCM6","CDCA7","PRIM1","UHRF1","HELLS","RFC2","RPA2","NASP","RAD51AP1","GMNN","SLBP","USP1","CLSPN","ATAD2")

g2m.genes.short <- c("HMGB2","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","CKS2","CKS1B","MKI67","TMPO","CENPF","TACC3","SMC4","CCNB2","CKAP2","ANP32E","TUBB4B","KIF20B","HN1","LBR","CTCF","CBX5")

pdf(paste(filedir,"_RidgePlot_CellCycle_afterRmDiff_S-genes-short_res0.05_defaultK",date,".pdf",sep=""), width=25, height=100)
dPlot <- RidgePlot(FullSetMultiRes_CC, features = s.genes.short, ncol = 3)
print(dPlot)
dev.off()

pdf(paste(filedir,"_RidgePlot_CellCycle_afterRmDiff_G2M-genes-short_res0.05_defaultK",date,".pdf",sep=""), width=25, height=100)
dPlot <- RidgePlot(FullSetMultiRes_CC, features = g2m.genes.short, ncol = 3)
print(dPlot)
dev.off()

pdf(paste(filedir,"_VlnPlot_CellCycle_afterRmDiff_S-genes-short_res0.05_defaultK",date,".pdf",sep=""), width=25, height=100)
dPlot <- VlnPlot(FullSetMultiRes_CC, features = s.genes.short, ncol = 3, group.by='Phase', pt.size=0)
print(dPlot)
dev.off()

pdf(paste(filedir,"_VlnPlot_CellCycle_afterRmDiff_G2M-genes-short_res0.05_defaultK",date,".pdf",sep=""), width=25, height=100)
dPlot <- VlnPlot(FullSetMultiRes_CC, features = g2m.genes.short, ncol = 3, group.by='Phase', pt.size=0)
print(dPlot)
dev.off()

#### CC corrected cluster 1 subset
clus1.subset <- subset(FullSetMultiRes_CC, subset= cluster == "1")
dim(clus1.subset)
# [1] 33694 25498

g1.clus1.subset <- subset(clus1.subset, subset= Phase == "G1")
dim(g1.clus1.subset)
# [1] 33694 10864

g2m.clus1.subset <- subset(clus1.subset, subset= Phase == "G2M")
dim(g2m.clus1.subset)
# [1] 33694  4843

s.clus1.subset <- subset(clus1.subset, subset= Phase == "S")
dim(s.clus1.subset)
# [1] 33694  9791

saveRDS(clus1.subset, paste(filedir,"SeuratObject_post-rmCC_Clus1-Only",date,".rds",sep=""))
saveRDS(g1.clus1.subset, paste(filedir,"SeuratObject_post-rmCC_Clus1_G1",date,".rds",sep=""))
saveRDS(g2m.clus1.subset, paste(filedir,"SeuratObject_post-rmCC_Clus1_G2M",date,".rds",sep=""))
saveRDS(s.clus1.subset, paste(filedir,"SeuratObject_post-rmCC_Clus1_S",date,".rds",sep=""))

### subset for looking at cell-cycle effects based on low vs high diversity groupings -- after correction; will do pre-correction after this (reloading uncorrected dataset)
## have to label Phase

clus1.subset.pre <- subset(FullSetMultiRes, subset= cluster == "1")
dim(clus1.subset.pre)
# [1] 33694 25498

FullSetMultiRes <- CellCycleScoring(FullSetMultiRes, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

g1.subset <- subset(FullSetMultiRes, subset= Phase == "G1")
dim(g1.subset)
# [1]  33694 110054

clus1.g1.subset <- subset(g1.subset, subset= cluster == "1")
dim(clus1.g1.subset)
# [1] 33694 10864

g2m.subset <- subset(FullSetMultiRes, subset= Phase == "G2M")
dim(g2m.subset)
# [1] 33694 22568

clus1.g2m.subset <- subset(g2m.subset, subset= cluster == "1")
dim(clus1.g2m.subset)
# [1] 33694  4843

s.subset <- subset(FullSetMultiRes, subset= Phase == "S")
dim(s.subset)
# [1] 33694 42092

clus1.s.subset <- subset(s.subset, subset= cluster == "1")
dim(clus1.s.subset)
# [1] 33694  9791

saveRDS(clus1.subset.pre, paste(filedir,"SeuratObject_pre-rmCC_Clus1-Only_All-Phases",date,".rds",sep=""))
saveRDS(clus1.g1.subset, paste(filedir,"SeuratObject_pre-rmCC_G1_Clus1-Only",date,".rds",sep=""))
saveRDS(clus1.g2m.subset, paste(filedir,"SeuratObject_pre-rmCC_G2M_Clus1-Only",date,".rds",sep=""))
saveRDS(clus1.s.subset, paste(filedir,"SeuratObject_pre-rmCC_S_Clus1-Only",date,".rds",sep=""))

### pre-rmCC exported cell breakdowns are fine --- don't need to re-do

# file1 <- "./analysis/CMMLoutputDataNames_31C+7Npts_CD34only_CC_res0.05_rmCC_20200710.csv";
# file2 <- "./analysis/CMMLoutputData_31C+7Npts_CD34only_CC_res0.05_rmCC_20200710.csv";

# cat(g1.subset@meta.data$RNA_snn_res.0.05, file=file2, sep=",\n")

# listNames <- g1.subset@meta.data$orig.ident
# write.table(data.frame(listNames),
#           row.names=FALSE,
#           col.names = FALSE, 
#           file = file1,
#           sep=",")

# mydat1 <- read.csv(file2)
# mydat2 <- read.csv(file1)
# fulldat <- cbind(mydat2[1],mydat1[1])
# fulltab <- as.data.table(fulldat)
# # name table columns
# names(fulltab)[1] <- paste("UMI")
# names(fulltab)[2] <- paste("cluster")
# # group data based on clusters
# group_by(fulltab, cluster)
# # create a table counting unqiue UMIs/cells per cluster
# tabPerClus <- fulltab %>% group_by(cluster) %>% count()
# type <- sub("\\_.*","",fulltab$UMI)
# fulltab <- cbind(fulltab, type)

# divout <-  paste(filedir,"_CellBreakdown_PerClusterPerType_res-0.05_pre-rmCC_G1-only",date,".csv",sep="") 
# A <- fulltab %>% group_by(cluster) %>% count(type)
# write.csv(A, file=divout)

# rm(listNames,mydat1,mydat2,fulldat,fulltab,tabPerClus,type,fulltab,divout,A)


# cat(g2m.subset@meta.data$RNA_snn_res.0.05, file=file2, sep=",\n")

# listNames <- g2m.subset@meta.data$orig.ident
# write.table(data.frame(listNames),
#           row.names=FALSE,
#           col.names = FALSE, 
#           file = file1,
#           sep=",")

# mydat1 <- read.csv(file2)
# mydat2 <- read.csv(file1)
# fulldat <- cbind(mydat2[1],mydat1[1])
# fulltab <- as.data.table(fulldat)
# # name table columns
# names(fulltab)[1] <- paste("UMI")
# names(fulltab)[2] <- paste("cluster")
# # group data based on clusters
# group_by(fulltab, cluster)
# # create a table counting unqiue UMIs/cells per cluster
# tabPerClus <- fulltab %>% group_by(cluster) %>% count()
# type <- sub("\\_.*","",fulltab$UMI)
# fulltab <- cbind(fulltab, type)

# divout <-  paste(filedir,"_CellBreakdown_PerClusterPerType_res-0.05_pre-rmCC_G2M-only",date,".csv",sep="") 
# A <- fulltab %>% group_by(cluster) %>% count(type)
# write.csv(A, file=divout)

# rm(listNames,mydat1,mydat2,fulldat,tabPerClus,type,fulltab,divout,A)


# cat(s.subset@meta.data$RNA_snn_res.0.05, file=file2, sep=",\n")

# listNames <- s.subset@meta.data$orig.ident
# write.table(data.frame(listNames),
#           row.names=FALSE,
#           col.names = FALSE, 
#           file = file1,
#           sep=",")

# mydat1 <- read.csv(file2)
# mydat2 <- read.csv(file1)
# fulldat <- cbind(mydat2[1],mydat1[1])
# fulltab <- as.data.table(fulldat)
# # name table columns
# names(fulltab)[1] <- paste("UMI")
# names(fulltab)[2] <- paste("cluster")
# # group data based on clusters
# group_by(fulltab, cluster)
# # create a table counting unqiue UMIs/cells per cluster
# tabPerClus <- fulltab %>% group_by(cluster) %>% count()
# type <- sub("\\_.*","",fulltab$UMI)
# fulltab <- cbind(fulltab, type)

# divout <-  paste(filedir,"_CellBreakdown_PerClusterPerType_res-0.05_pre-rmCC_S-only",date,".csv",sep="") 
# A <- fulltab %>% group_by(cluster) %>% count(type)
# write.csv(A, file=divout)

# rm(listNames,mydat1,mydat2,fulldat,tabPerClus,type,fulltab,divout,A)

##### 

### NEXT UP 
### GROUPING HMA TOGETHER
# highlighted UMAP
# differential gene expression

date <- "_2020-07-23"

hma1 <- WhichCells(FullSetMultiRes, expression = name =="4K001")
hma2 <- WhichCells(FullSetMultiRes, expression = name =="5E001")
hma3 <- WhichCells(FullSetMultiRes, expression = name =="4S001")
hma4 <- WhichCells(FullSetMultiRes, expression = name =="4Q001")

length(hma1)  # [1] 2304
length(hma2)  # [1] 2054
length(hma3)  # [1] 2607
length(hma4)  # [1] 4730

hma <- c(hma1,hma2,hma3,hma4)

pdf(paste(filedir,"UMAP_Clus1_res0.6_Highlighted_HMA",date,".pdf",sep=""))
dPlot <- DimPlot(clus1.subset.pre, reduction.use="umap", cells.highlight=hma) + theme(legend.position = "none")
print(dPlot)
dev.off()

clus1g1 <- WhichCells(clus1.subset, expression = Phase =="G1")
clus1g2m <- WhichCells(clus1.subset, expression = Phase =="G2M")
clus1s <- WhichCells(clus1.subset, expression = Phase =="S")

pdf(paste(filedir,"UMAP_Clus1_res0.6_Highlighted_G1",date,".pdf",sep=""))
dPlot <- DimPlot(clus1.subset.pre, reduction.use="umap", cells.highlight=clus1g1) + theme(legend.position = "none")
print(dPlot)
dev.off()

pdf(paste(filedir,"UMAP_Clus1_res0.6_Highlighted_G2M",date,".pdf",sep=""))
dPlot <- DimPlot(clus1.subset.pre, reduction.use="umap", cells.highlight=clus1g2m) + theme(legend.position = "none")
print(dPlot)
dev.off()

pdf(paste(filedir,"UMAP_Clus1_res0.6_Highlighted_S",date,".pdf",sep=""))
dPlot <- DimPlot(clus1.subset.pre, reduction.use="umap", cells.highlight=clus1s) + theme(legend.position = "none")
print(dPlot)
dev.off()

## other CMMLs
l1 <- WhichCells(FullSetMultiRes, expression = name=="SF16072200003")
l3 <- WhichCells(FullSetMultiRes, expression = name=="SF16112300029")
l5 <- WhichCells(FullSetMultiRes, expression = name=="SF14050700419")
l6 <- WhichCells(FullSetMultiRes, expression = name=="SF12042500035")
l7 <- WhichCells(FullSetMultiRes, expression = name=="LTB3966")
l8 <- WhichCells(FullSetMultiRes, expression = name=="SF100109101914")
l11 <- WhichCells(FullSetMultiRes, expression = name=="2V001")
l12 <- WhichCells(FullSetMultiRes, expression = name=="LTB4121")
h1 <- WhichCells(FullSetMultiRes, expression = name=="SF100109111451")
h2 <- WhichCells(FullSetMultiRes, expression = name=="LTB6169")
h3 <- WhichCells(FullSetMultiRes, expression = name=="SF14031800065")
h4 <- WhichCells(FullSetMultiRes, expression = name=="SF100109106293")
h5 <- WhichCells(FullSetMultiRes, expression = name=="4J003")
h6 <- WhichCells(FullSetMultiRes, expression = name=="6AC001")
h7 <- WhichCells(FullSetMultiRes, expression = name=="SF14040100158")
h8 <- WhichCells(FullSetMultiRes, expression = name=="SF16026800045")
h9 <- WhichCells(FullSetMultiRes, expression = name=="6AE001")
h10 <- WhichCells(FullSetMultiRes, expression = name=="SF16112900158")
h11 <- WhichCells(FullSetMultiRes, expression = name=="5H001")
h12 <- WhichCells(FullSetMultiRes, expression = name=="SF100109110236")
h13 <- WhichCells(FullSetMultiRes, expression = name=="SF12092600014")
h14 <- WhichCells(FullSetMultiRes, expression = name=="SF14101000049")
h15 <- WhichCells(FullSetMultiRes, expression = name=="SF12062800475")
h16 <- WhichCells(FullSetMultiRes, expression = name=="SF14072200012")
h17 <- WhichCells(FullSetMultiRes, expression = name=="SF13061200056")
h18 <- WhichCells(FullSetMultiRes, expression = name=="6AD001")
h19 <- WhichCells(FullSetMultiRes, expression = name=="SF14060200025")

treatmentnaive <- c(l1,l3,l5,l6,l7,l8,l11,l12,h1,h2,h3,h4,h6,h7,h8,h9,h10,h12,h13,h14,h15,h16,h17,h18,h19)
treatmentnaive_withHydrea <- c(treatmentnaive,h11)

n1 <- WhichCells(FullSetMultiRes, expression = name=="CD34")
n2 <- WhichCells(FullSetMultiRes, expression = name=="SettyPt1")
n3 <- WhichCells(FullSetMultiRes, expression = name=="SettyPt2")
n4 <- WhichCells(FullSetMultiRes, expression = name=="SettyPt3")
n5 <- WhichCells(FullSetMultiRes, expression = name=="HuaPt1")
n6 <- WhichCells(FullSetMultiRes, expression = name=="HuaPt2")
n7 <- WhichCells(FullSetMultiRes, expression = name=="HuaPt4")

normals <- c(n1,n2,n3,n4,n5,n6,n7)
normals_withBMT <- c(normals,h5)

### Differential Gene Expression
new.de.markers.1 <- FindMarkers(FullSetMultiRes, ident.1=hma, ident.2=treatmentnaive) ## Default by Wilcoxon Rank Sum Test
write.csv(new.de.markers.1, paste(filedir,"DiffExpGeneList_res0.05_1-HMA-vs-2-TxNaive_viaWilcoxonRankTest",date,".csv",sep=""))

new.de.markers.2 <- FindMarkers(FullSetMultiRes, ident.1=hma, ident.2=treatmentnaive_withHydrea) 
write.csv(new.de.markers.2, paste(filedir,"DiffExpGeneList_res0.05_1-HMA-vs-2-TxNaive+Hydrea_viaWilcoxonRankTest",date,".csv",sep=""))

new.de.markers.3 <- FindMarkers(FullSetMultiRes, ident.1=hma, ident.2=normals) 
write.csv(new.de.markers.3, paste(filedir,"DiffExpGeneList_res0.05_1-HMA-vs-2-Normals_viaWilcoxonRankTest",date,".csv",sep=""))

new.de.markers.4 <- FindMarkers(FullSetMultiRes, ident.1=hma, ident.2=normals_withBMT) 
write.csv(new.de.markers.4, paste(filedir,"DiffExpGeneList_res0.05_1-HMA-vs-2-Normals+BMT_viaWilcoxonRankTest",date,".csv",sep=""))



#### SUBSET DATA AND SAVE --- to be done 
hma_cells <- subset(x=FullSetMultiRes, cells=hma)

data_to_write_out_hma <- as.data.frame(as.matrix(hma_cells@assays$RNA@data))
dim(data_to_write_out_hma)
# [1] 33694 11695
fwrite(x = data_to_write_out_hma, file =paste(filedir,"CMML-scRNAseq_HMA-only_norm-counts-mtx_withHeaders",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)

data_embeddings_hma <- as.data.frame(as.matrix(Embeddings(object=hma_cells[["umap"]])))
fwrite(x = data_embeddings_hma, file =paste(filedir,"CMML-scRNAseq_HMA-only_UMAP-Embeddings_withHeaders",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)

data_clusters_hma <- as.data.frame(as.matrix(hma_cells[["RNA_snn_res.0.05"]]))
fwrite(x = data_clusters_hma, file =paste(filedir,"CMML-scRNAseq_HMA-only_Clusters_withHeaders",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)

saveRDS(hma_cells, paste(filedir,"SeuratObject_Subset-HMA-only_withCCphase",date,".rds",sep=""))

treatmentnaive_cells <- subset(x=FullSetMultiRes, cells=treatmentnaive)
dim(treatmentnaive_cells)
# [1] 33694 95564

data_to_write_out_naive1 <- as.data.frame(as.matrix(treatmentnaive_cells@assays$RNA@data[,1:50000]))
dim(data_to_write_out_naive)
# [1] 33694 50000
data_to_write_out_naive2 <- as.data.frame(as.matrix(treatmentnaive_cells@assays$RNA@data[,50001:95564]))

head(data_to_write_out_naive2)
# rows are genes
# columns are samples

data_to_write_out_naive <- cbind(data_to_write_out_naive1,data_to_write_out_naive2)

fwrite(x = data_to_write_out_naive, file =paste(filedir,"CMML-scRNAseq_Tx-naive-only_norm-counts-mtx_withHeaders",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)

data_embeddings_naive <- as.data.frame(as.matrix(Embeddings(object=treatmentnaive_cells[["umap"]])))
fwrite(x = data_embeddings_hma, file =paste(filedir,"CMML-scRNAseq_Tx-naive-only_UMAP-Embeddings_withHeaders",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)

data_clusters_naive <- as.data.frame(as.matrix(treatmentnaive_cells[["RNA_snn_res.0.05"]]))
fwrite(x = data_clusters_naive, file =paste(filedir,"CMML-scRNAseq_Tx-naive-only_Clusters_withHeaders",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)

saveRDS(treatmentnaive_cells, paste(filedir,"SeuratObject_Subset-Tx-naive-only_withCCphase",date,".rds",sep=""))

treatmentnaive_hydrea_cells <- subset(x=FullSetMultiRes, cells=treatmentnaive_withHydrea)
dim(treatmentnaive_hydrea_cells)
# [1] 33694 99660

data_to_write_out_naive_hydrea1 <- as.data.frame(as.matrix(treatmentnaive_hydrea_cells@assays$RNA@data[,1:50000]))
data_to_write_out_naive_hydrea2 <- as.data.frame(as.matrix(treatmentnaive_hydrea_cells@assays$RNA@data[,50001:99660]))
data_to_write_out_naive_hydrea <- cbind(data_to_write_out_naive_hydrea1,data_to_write_out_naive_hydrea2)

fwrite(x = data_to_write_out_naive_hydrea, file =paste(filedir,"CMML-scRNAseq_Tx-naive+hydrea-only_norm-counts-mtx_withHeaders",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)

data_embeddings_naive_hydrea <- as.data.frame(as.matrix(Embeddings(object=treatmentnaive_hydrea_cells[["umap"]])))
fwrite(x = data_embeddings_hma, file =paste(filedir,"CMML-scRNAseq_Tx-naive+hydrea-only_UMAP-Embeddings_withHeaders",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)

data_clusters_naive_hydrea <- as.data.frame(as.matrix(treatmentnaive_hydrea_cells[["RNA_snn_res.0.05"]]))
fwrite(x = data_clusters_naive_hydrea, file =paste(filedir,"CMML-scRNAseq_Tx-naive+hydrea-only_Clusters_withHeaders",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)

saveRDS(treatmentnaive_hydrea_cells, paste(filedir,"SeuratObject_Subset-Tx-naive+hydrea-only_withCCphase",date,".rds",sep=""))

normal_cells <- subset(x=FullSetMultiRes, cells=normals)
dim(normal_cells)
[1] 33694 59775

data_to_write_out_normal <- as.data.frame(as.matrix(normal_cells@assays$RNA@data))
fwrite(x = data_to_write_out_normal, file =paste(filedir,"CMML-scRNAseq_Normal-only_norm-counts-mtx_withHeaders",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)

data_embeddings_normal <- as.data.frame(as.matrix(Embeddings(object=normal_cells[["umap"]])))
fwrite(x = data_embeddings_normal, file =paste(filedir,"CMML-scRNAseq_Normal-only_UMAP-Embeddings_withHeaders",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)

data_clusters_normal <- as.data.frame(as.matrix(normal_cells[["RNA_snn_res.0.05"]]))
fwrite(x = data_clusters_normal, file =paste(filedir,"CMML-scRNAseq_Normal-only_Clusters_withHeaders",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)

saveRDS(normal_cells, paste(filedir,"SeuratObject_Subset-Normal-only_withCCphase",date,".rds",sep=""))

normal_BMT_cells <- subset(x=FullSetMultiRes, cells=normals_withBMT)
dim(normal_BMT_cells)
[1] 33694 63359

data_to_write_out_normal_BMT <- as.data.frame(as.matrix(normal_BMT_cells@assays$RNA@data))
fwrite(x = data_to_write_out_normal_BMT file =paste(filedir,"CMML-scRNAseq_Normal+BMT-only_norm-counts-mtx_withHeaders",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)

data_embeddings_normal_BMT <- as.data.frame(as.matrix(Embeddings(object=normal_BMT_cells[["umap"]])))
fwrite(x = data_embeddings_normal, file =paste(filedir,"CMML-scRNAseq_Normal+BMT-only_UMAP-Embeddings_withHeaders",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)

data_clusters_normal_BMT <- as.data.frame(as.matrix(normal_BMT_cells[["RNA_snn_res.0.05"]]))
fwrite(x = data_clusters_normal_BMT, file =paste(filedir,"CMML-scRNAseq_Normal+BMT-only_Clusters_withHeaders",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)

saveRDS(normal_BMT_cells, paste(filedir,"SeuratObject_Subset-Normal+BMT-only_withCCphase",date,".rds",sep=""))

### was able to continue on this loaded data for average expression in cluster 1 per patient for GSEA analysis
date <- "_2020-07-27"
filedir <- "./analysis/GSEA/Healthy-7pts+CMML-31pts-Cohort_"

Idents(FullSetMultiRes) <- FullSetMultiRes$cluster

hma_1_cells <- subset(x=FullSetMultiRes, cells=hma1) #4K001
dim(hma_1_cells)
# [1] 33694  2304
hma1.cluster.averages <- AverageExpression(hma_1_cells)
head(hma1.cluster.averages[["RNA"]][, 1:5])
#                    0          1 12     15         4
# TSPAN6   0.007234581 0.03016423  0 0.0000 0.0000000
# TNMD     0.000000000 0.00000000  0 0.0000 0.0000000
# DPM1     0.234642726 0.31390147  0 2.0803 0.0000000
# SCYL3    0.039560018 0.03968130  0 0.0000 0.4726791
# C1orf112 0.105491971 0.12130262  0 0.0000 0.0000000
# FGR      0.000000000 0.00219081  0 0.0000 0.0000000
hma1.avgClus1 <- cbind(rownames(hma1.cluster.averages[["RNA"]]), hma1.cluster.averages[["RNA"]][, 2])
data.hma1 <- as.data.frame(as.matrix(hma1.avgClus1))
fwrite(x = data.hma1, file =paste(filedir,"Export_Cluster-1-Avg-Exp_HMA1_4-K-001",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


hma_2_cells <- subset(x=FullSetMultiRes, cells=hma2) #5E001
dim(hma_2_cells)
# [1] 33694  2054
hma2.cluster.averages <- AverageExpression(hma_2_cells)
head(hma2.cluster.averages[["RNA"]][, 1:5])
#                   0           1 12        15         2
# TSPAN6   0.00000000 0.003783236  0 0.0000000 0.0000000
# TNMD     0.00000000 0.000000000  0 0.0000000 0.0000000
# DPM1     0.48034504 0.319121878  0 0.1128605 0.3813883
# SCYL3    0.07426206 0.053833045  0 0.0000000 0.0000000
# C1orf112 0.16497705 0.068131201  0 0.0000000 0.0000000
# FGR      0.04085871 0.013152146  0 0.0000000 0.0000000
hma2.avgClus1 <- cbind(rownames(hma2.cluster.averages[["RNA"]]), hma2.cluster.averages[["RNA"]][, 2])
data.hma2 <- as.data.frame(as.matrix(hma2.avgClus1))
fwrite(x = data.hma2, file =paste(filedir,"Export_Cluster-1-Avg-Exp_HMA2_5-E-001",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


hma_3_cells <- subset(x=FullSetMultiRes, cells=hma3) #4S001
dim(hma_3_cells)
# [1] 33694  2607
hma3.cluster.averages <- AverageExpression(hma_3_cells)
head(hma3.cluster.averages[["RNA"]][, 1:5])
#                    0           1 12       15        16
# TSPAN6   0.005600552 0.009678594  0 0.000000 0.5516568
# TNMD     0.000000000 0.000000000  0 0.000000 0.0000000
# DPM1     0.176823317 0.190782916  0 2.611307 0.1625917
# SCYL3    0.070469616 0.038021538  0 0.000000 0.0000000
# C1orf112 0.053799257 0.068724043  0 0.000000 0.2121880
# FGR      0.006226809 0.006938659  0 0.000000 0.0000000
hma3.avgClus1 <- cbind(rownames(hma3.cluster.averages[["RNA"]]), hma3.cluster.averages[["RNA"]][, 2])
data.hma3 <- as.data.frame(as.matrix(hma3.avgClus1))
fwrite(x = data.hma3, file =paste(filedir,"Export_Cluster-1-Avg-Exp_HMA3_4-S-001",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


hma_4_cells <- subset(x=FullSetMultiRes, cells=hma4) #4Q001
dim(hma_4_cells)
# [1] 33694  4730
hma4.cluster.averages <- AverageExpression(hma_4_cells)
head(hma4.cluster.averages[["RNA"]][, 1:5])
#                   0          1 12 15          3
# TSPAN6   0.00000000 0.01273291  0  0 0.00000000
# TNMD     0.00000000 0.00000000  0  0 0.00000000
# DPM1     0.27087284 0.28332196  0  0 0.25048694
# SCYL3    0.07461141 0.05310830  0  0 0.01483263
# C1orf112 0.08153525 0.05741517  0  0 0.08191321
# FGR      0.09910226 0.01556318  0  0 3.18984681
hma4.avgClus1 <- cbind(rownames(hma4.cluster.averages[["RNA"]]), hma4.cluster.averages[["RNA"]][, 2])
data.hma4 <- as.data.frame(as.matrix(hma4.avgClus1))
fwrite(x = data.hma4, file =paste(filedir,"Export_Cluster-1-Avg-Exp_HMA4_4-Q-001",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


naive_1_cells <- subset(x=FullSetMultiRes, cells=l1) #SF16072200003
dim(naive_1_cells)
# 
naive1.cluster.averages <- AverageExpression(naive_1_cells)
head(naive1.cluster.averages[["RNA"]][, 1:5])
#                   0           1          12       15          3
# TSPAN6   0.01274699 0.005123409 0.007186773 0.000000 0.00000000
# TNMD     0.00000000 0.000000000 0.000000000 0.000000 0.00000000
# DPM1     0.27506847 0.362446648 0.274120839 0.000000 0.14536166
# SCYL3    0.06521780 0.059249004 0.074251085 0.000000 0.03518174
# C1orf112 0.05539055 0.083288226 0.055915398 0.000000 0.02054351
# FGR      0.01880102 0.007414434 0.000000000 1.565435 2.91932328
naive1.avgClus1 <- cbind(rownames(naive1.cluster.averages[["RNA"]]), naive1.cluster.averages[["RNA"]][, 2])
data.naive1 <- as.data.frame(as.matrix(naive1.avgClus1))
fwrite(x = data.naive1, file =paste(filedir,"Export_Cluster-1-Avg-Exp_Naive-1_SF-160722-00003",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


naive_2_cells <- subset(x=FullSetMultiRes, cells=l3) #SF16112300029
dim(naive_2_cells)
# [1] 33694  1620
naive2.cluster.averages <- AverageExpression(naive_2_cells)
head(naive2.cluster.averages[["RNA"]][, 1:5])
#                   0          1        12 15 2
# TSPAN6   0.02198011 0.02652541 0.0000000  0 0
# TNMD     0.00000000 0.00000000 0.0000000  0 0
# DPM1     0.32831504 0.49183929 0.0000000  0 0
# SCYL3    0.11093524 0.07369301 0.8110300  0 0
# C1orf112 0.03441661 0.04222590 0.0000000  0 0
# FGR      0.02612626 0.04454616 0.7659314  0 0

naive2.avgClus1 <- cbind(rownames(naive2.cluster.averages[["RNA"]]), naive2.cluster.averages[["RNA"]][, 2])
data.naive2 <- as.data.frame(as.matrix(naive2.avgClus1))
fwrite(x = data.naive2, file =paste(filedir,"Export_Cluster-1-Avg-Exp_Naive-2_SF-161123-00029",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


naive_3_cells <- subset(x=FullSetMultiRes, cells=l5) #SF14050700419
dim(naive_3_cells)
# [1] 33694  2788
naive3.cluster.averages <- AverageExpression(naive_3_cells)
head(naive3.cluster.averages[["RNA"]][, 1:5])
#                    0           1        12 2          3
# TSPAN6   0.006084634 0.001453669 0.0000000 0 0.00000000
# TNMD     0.000000000 0.000000000 0.0000000 0 0.00000000
# DPM1     0.320065027 0.394724853 0.1675126 0 0.20622189
# SCYL3    0.094123778 0.057697021 0.0000000 0 0.05213006
# C1orf112 0.053075707 0.081295699 0.0000000 0 0.01219732
# FGR      0.034120698 0.023812846 0.5741627 0 4.21482847
naive3.avgClus1 <- cbind(rownames(naive3.cluster.averages[["RNA"]]), naive3.cluster.averages[["RNA"]][, 2])
data.naive3 <- as.data.frame(as.matrix(naive3.avgClus1))
fwrite(x = data.naive3, file =paste(filedir,"Export_Cluster-1-Avg-Exp_Naive-3_SF-140507-00419",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


naive_4_cells <- subset(x=FullSetMultiRes, cells=l6) #SF12042500035
dim(naive_4_cells)
# [1] 33694  5072
naive4.cluster.averages <- AverageExpression(naive_4_cells)
head(naive4.cluster.averages[["RNA"]][, 1:5])
#                   0           1        12 14          3
# TSPAN6   0.01019038 0.009888487 0.0000000  0 0.00000000
# TNMD     0.00000000 0.000000000 0.0000000  0 0.00000000
# DPM1     0.25069729 0.316838771 0.6108392  0 0.25899684
# SCYL3    0.08031383 0.051886840 0.0000000  0 0.04025762
# C1orf112 0.05641343 0.057660760 0.0432324  0 0.02348380
# FGR      0.04725493 0.018877666 0.4502268  0 3.21964908
naive4.avgClus1 <- cbind(rownames(naive4.cluster.averages[["RNA"]]), naive4.cluster.averages[["RNA"]][, 2])
data.naive4 <- as.data.frame(as.matrix(naive4.avgClus1))
fwrite(x = data.naive4, file =paste(filedir,"Export_Cluster-1-Avg-Exp_Naive-4_SF-120425-00035",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


naive_5_cells <- subset(x=FullSetMultiRes, cells=l7) #LTB3966
dim(naive_5_cells)
# [1] 33694  4051
naive5.cluster.averages <- AverageExpression(naive_5_cells)
head(naive5.cluster.averages[["RNA"]][, 1:5])
#                    0           1 12 2          3
# TSPAN6   0.006043723 0.009029238  0 0 0.00000000
# TNMD     0.000000000 0.000000000  0 0 0.00000000
# DPM1     0.253346563 0.308008021  0 0 0.17938831
# SCYL3    0.038155665 0.039028594  0 0 0.03002017
# C1orf112 0.062994222 0.081826966  0 0 0.00000000
# FGR      0.012256593 0.009564437  0 0 1.99864313
naive5.avgClus1 <- cbind(rownames(naive5.cluster.averages[["RNA"]]), naive5.cluster.averages[["RNA"]][, 2])
data.naive5 <- as.data.frame(as.matrix(naive5.avgClus1))
fwrite(x = data.naive5, file =paste(filedir,"Export_Cluster-1-Avg-Exp_Naive-5_LTB3966",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


naive_6_cells <- subset(x=FullSetMultiRes, cells=l8) #SF100109101914
dim(naive_6_cells)
# [1] 33694  2278
naive6.cluster.averages <- AverageExpression(naive_6_cells)
head(naive6.cluster.averages[["RNA"]][, 1:5])
#                   0           1         12        15         16
# TSPAN6   0.00000000 0.001301782 0.00000000 0.0000000 0.50676387
# TNMD     0.00000000 0.000000000 0.00000000 0.0000000 0.00000000
# DPM1     0.39697304 0.278051756 0.32363086 0.3744064 0.47828049
# SCYL3    0.04104134 0.039320632 0.05274973 0.0000000 0.08306425
# C1orf112 0.02189693 0.019821646 0.01395463 0.0000000 0.02243460
# FGR      0.02380162 0.036625387 0.01867936 0.3962642 0.00000000
naive6.avgClus1 <- cbind(rownames(naive6.cluster.averages[["RNA"]]), naive6.cluster.averages[["RNA"]][, 2])
data.naive6 <- as.data.frame(as.matrix(naive6.avgClus1))
fwrite(x = data.naive6, file =paste(filedir,"Export_Cluster-1-Avg-Exp_Naive-6_SF-100109-101914",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)

naive_7_cells <- subset(x=FullSetMultiRes, cells=l11) #2V001
dim(naive_7_cells)
# [1] 33694  3162
naive7.cluster.averages <- AverageExpression(naive_7_cells)
head(naive7.cluster.averages[["RNA"]][, 1:5])
#                   0          1         12        15        16
# TSPAN6   0.00162159 0.01380210 0.00000000 0.0000000 0.5222901
# TNMD     0.00000000 0.00000000 0.00000000 0.0000000 0.0000000
# DPM1     0.21025232 0.23543689 0.05018984 0.0000000 0.2689329
# SCYL3    0.04198603 0.07562230 0.03816852 0.1666861 0.0000000
# C1orf112 0.04994461 0.07204052 0.00000000 0.0000000 0.0000000
# FGR      0.12611827 0.02332013 0.33442797 0.1454948 0.0000000
naive7.avgClus1 <- cbind(rownames(naive7.cluster.averages[["RNA"]]), naive7.cluster.averages[["RNA"]][, 2])
data.naive7 <- as.data.frame(as.matrix(naive7.avgClus1))
fwrite(x = data.naive7, file =paste(filedir,"Export_Cluster-1-Avg-Exp_Naive-7_2-V-001",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


naive_8_cells <- subset(x=FullSetMultiRes, cells=l12) #LTB4121
dim(naive_8_cells)
# [1] 33694  4364
naive8.cluster.averages <- AverageExpression(naive_8_cells)
head(naive8.cluster.averages[["RNA"]][, 1:5])
#                    0          1        12 17 2
# TSPAN6   0.000683091 0.05640879 0.0000000  0 0
# TNMD     0.000000000 0.00000000 0.0000000  0 0
# DPM1     0.223374574 0.24594300 0.1384943  0 0
# SCYL3    0.044540427 0.04347802 0.0000000  0 0
# C1orf112 0.057143896 0.07004201 0.0000000  0 0
# FGR      0.023735222 0.01407687 0.0000000  0 0
naive8.avgClus1 <- cbind(rownames(naive8.cluster.averages[["RNA"]]), naive8.cluster.averages[["RNA"]][, 2])
data.naive8 <- as.data.frame(as.matrix(naive8.avgClus1))
fwrite(x = data.naive8, file =paste(filedir,"Export_Cluster-1-Avg-Exp_Naive-8_LTB4121",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


naive_9_cells <- subset(x=FullSetMultiRes, cells=h1) #SF100109111451
dim(naive_9_cells)
# [1] 33694  3881
naive9.cluster.averages <- AverageExpression(naive_9_cells)
head(naive9.cluster.averages[["RNA"]][, 1:5])
#                    0            1 12 15        17
# TSPAN6   0.003970496 0.0039690263  0  0 0.5712615
# TNMD     0.000000000 0.0004913507  0  0 0.0000000
# DPM1     0.165490852 0.2380959351  0  0 0.0000000
# SCYL3    0.066247579 0.0392933483  0  0 0.0000000
# C1orf112 0.059640636 0.0862615281  0  0 0.0000000
# FGR      0.044661433 0.0086839669  0  0 0.0000000
naive9.avgClus1 <- cbind(rownames(naive9.cluster.averages[["RNA"]]), naive9.cluster.averages[["RNA"]][, 2])
data.naive9 <- as.data.frame(as.matrix(naive9.avgClus1))
fwrite(x = data.naive9, file =paste(filedir,"Export_Cluster-1-Avg-Exp_Naive-9_SF-100109-111451",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


naive_10_cells <- subset(x=FullSetMultiRes, cells=h2) #LTB6169
dim(naive_10_cells)
# [1] 33694  3582
naive10.cluster.averages <- AverageExpression(naive_10_cells)
head(naive10.cluster.averages[["RNA"]][, 1:5])
#                    0          1        12       15          3
# TSPAN6   0.001486215 0.02575127 0.0000000 0.000000 0.00000000
# TNMD     0.000000000 0.00000000 0.0000000 0.000000 0.00000000
# DPM1     0.275702114 0.26031930 0.0000000 1.950458 0.27878302
# SCYL3    0.065390525 0.06649800 0.2734482 0.000000 0.09254996
# C1orf112 0.106899144 0.09285609 0.0000000 0.000000 0.06232479
# FGR      0.025617457 0.01151038 0.0000000 1.950458 2.59134839
naive10.avgClus1 <- cbind(rownames(naive10.cluster.averages[["RNA"]]), naive10.cluster.averages[["RNA"]][, 2])
data.naive10 <- as.data.frame(as.matrix(naive10.avgClus1))
fwrite(x = data.naive10, file =paste(filedir,"Export_Cluster-1-Avg-Exp_Naive-10_LTB5109",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


naive_11_cells <- subset(x=FullSetMultiRes, cells=h3) #SF14031800065
dim(naive_11_cells)
# [1] 33694  3502
naive11.cluster.averages <- AverageExpression(naive_11_cells)
head(naive11.cluster.averages[["RNA"]][, 1:5])
#                   0          1        11         12        15
# TSPAN6   0.03738067 0.01411362 0.0000000 0.00000000 0.0000000
# TNMD     0.00000000 0.00000000 0.0000000 0.00000000 0.0000000
# DPM1     0.31600711 0.62558431 0.5909118 0.24581776 0.2196065
# SCYL3    0.09027160 0.05317751 0.0000000 0.06906280 0.0000000
# C1orf112 0.04909308 0.04703194 0.0000000 0.01101615 0.0000000
# FGR      0.09743455 0.28761699 0.0000000 0.58719372 0.9390471
naive11.avgClus1 <- cbind(rownames(naive11.cluster.averages[["RNA"]]), naive11.cluster.averages[["RNA"]][, 2])
data.naive11 <- as.data.frame(as.matrix(naive11.avgClus1))
fwrite(x = data.naive11, file =paste(filedir,"Export_Cluster-1-Avg-Exp_Naive-11_SF-140318-00065",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)

naive_12_cells <- subset(x=FullSetMultiRes, cells=h4) #SF100109106293
dim(naive_12_cells)
# [1] 33694  3458
naive12.cluster.averages <- AverageExpression(naive_12_cells)
head(naive12.cluster.averages[["RNA"]][, 1:5])
#                   0          1         12 15       17
# TSPAN6   0.01142206 0.01894317 0.04203447  0 0.000000
# TNMD     0.00000000 0.00000000 0.00000000  0 0.000000
# DPM1     0.16707164 0.13564914 0.44403090  0 2.409436
# SCYL3    0.06752178 0.05371309 0.03757534  0 0.000000
# C1orf112 0.07416604 0.16863492 0.00000000  0 0.000000
# FGR      0.06936096 0.06085620 0.44080192  0 0.000000
naive12.avgClus1 <- cbind(rownames(naive12.cluster.averages[["RNA"]]), naive12.cluster.averages[["RNA"]][, 2])
data.naive12 <- as.data.frame(as.matrix(naive12.avgClus1))
fwrite(x = data.naive12, file =paste(filedir,"Export_Cluster-1-Avg-Exp_Naive-12_SF-100109-106293",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


naive_13_cells <- subset(x=FullSetMultiRes, cells=h6) #6AC001
dim(naive_13_cells)
# [1] 33694  4022
naive13.cluster.averages <- AverageExpression(naive_13_cells)
head(naive13.cluster.averages[["RNA"]][, 1:5])
#                   0           1          12        15        16
# TSPAN6   0.01838568 0.003373438 0.005265805 0.0000000 0.2561704
# TNMD     0.00000000 0.000000000 0.000000000 0.0000000 0.0000000
# DPM1     0.19337039 0.261584249 0.139219139 0.2150491 0.1115972
# SCYL3    0.04076651 0.052667196 0.005265805 0.0000000 0.0000000
# C1orf112 0.04114054 0.072373281 0.020499349 0.0000000 0.0000000
# FGR      0.04502903 0.012720796 0.133600652 0.2590003 0.0000000
naive13.avgClus1 <- cbind(rownames(naive13.cluster.averages[["RNA"]]), naive13.cluster.averages[["RNA"]][, 2])
data.naive13 <- as.data.frame(as.matrix(naive13.avgClus1))
fwrite(x = data.naive13, file =paste(filedir,"Export_Cluster-1-Avg-Exp_Naive-13_6-AC-001",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


naive_14_cells <- subset(x=FullSetMultiRes, cells=h7) #SF14040100158
dim(naive_14_cells)
# [1] 33694  4144
naive14.cluster.averages <- AverageExpression(naive_14_cells)
head(naive14.cluster.averages[["RNA"]][, 1:5])
#                   0           1           11        12 2
# TSPAN6   0.00000000 0.000000000 0.0002258927 0.0000000 0
# TNMD     0.00000000 0.000000000 0.0000000000 0.0000000 0
# DPM1     0.15071214 0.242408472 0.1332101034 0.3584147 0
# SCYL3    0.05529704 0.044975359 0.0644147198 0.0000000 0
# C1orf112 0.08234892 0.076683188 0.0396626394 0.0000000 0
# FGR      0.12550891 0.002387109 0.0298279247 0.0000000 0
naive14.avgClus1 <- cbind(rownames(naive14.cluster.averages[["RNA"]]), naive14.cluster.averages[["RNA"]][, 2])
data.naive14 <- as.data.frame(as.matrix(naive14.avgClus1))
fwrite(x = data.naive14, file =paste(filedir,"Export_Cluster-1-Avg-Exp_Naive-14_SF-140401-00158",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


naive_15_cells <- subset(x=FullSetMultiRes, cells=h8) #SF16026800045
dim(naive_15_cells)
# [1] 33694  4756
naive15.cluster.averages <- AverageExpression(naive_15_cells)
head(naive15.cluster.averages[["RNA"]][, 1:5])
#                    0            1        12       15        2
# TSPAN6   0.002500108 0.0623094433 0.0000000 0.000000 0.000000
# TNMD     0.000000000 0.0001823599 0.0000000 0.000000 0.000000
# DPM1     0.287255376 0.3929671041 0.1509888 0.000000 0.225469
# SCYL3    0.060174652 0.0338623820 0.0247480 0.000000 0.000000
# C1orf112 0.057486803 0.0770346215 0.0000000 0.000000 0.000000
# FGR      0.053622622 0.0062262373 0.0000000 6.259128 0.000000
naive15.avgClus1 <- cbind(rownames(naive15.cluster.averages[["RNA"]]), naive15.cluster.averages[["RNA"]][, 2])
data.naive15 <- as.data.frame(as.matrix(naive15.avgClus1))
fwrite(x = data.naive15, file =paste(filedir,"Export_Cluster-1-Avg-Exp_Naive-15_SF-160268-00045",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


naive_16_cells <- subset(x=FullSetMultiRes, cells=h9) #6AE001
dim(naive_16_cells)
# [1] 33694  4575
naive16.cluster.averages <- AverageExpression(naive_16_cells)
# head(naive16.cluster.averages[["RNA"]][, 1:5])
#                   0           1        12 16 2
# TSPAN6   0.00000000 0.005142797 0.0000000  0 0
# TNMD     0.00000000 0.000000000 0.0000000  0 0
# DPM1     0.47393763 0.215080697 0.0000000  0 0
# SCYL3    0.08977847 0.224502949 0.7823502  0 0
# C1orf112 0.05671092 0.007823716 0.0000000  0 0
# FGR      0.19916727 0.085150509 1.5577714  0 0
naive16.avgClus1 <- cbind(rownames(naive16.cluster.averages[["RNA"]]), naive16.cluster.averages[["RNA"]][, 2])
data.naive16 <- as.data.frame(as.matrix(naive16.avgClus1))
fwrite(x = data.naive16, file =paste(filedir,"Export_Cluster-1-Avg-Exp_Naive-16_6-AE-001",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


naive_17_cells <- subset(x=FullSetMultiRes, cells=h10) #SF16112900158
dim(naive_17_cells)
# [1] 33694  5043
naive17.cluster.averages <- AverageExpression(naive_17_cells)
head(naive17.cluster.averages[["RNA"]][, 1:5])
#                    0          1         12         15        16
# TSPAN6   0.001642434 0.02744533 0.01113680 0.00000000 0.4476037
# TNMD     0.000000000 0.00000000 0.00000000 0.00000000 0.0000000
# DPM1     0.306685743 0.26019271 0.18709127 0.16752075 0.4242809
# SCYL3    0.115085545 0.05162863 0.07404314 0.11002007 0.0000000
# C1orf112 0.069139889 0.11096421 0.02718741 0.02698137 0.0000000
# FGR      0.046225746 0.00767785 0.29604288 0.31555151 0.0000000
naive17.avgClus1 <- cbind(rownames(naive17.cluster.averages[["RNA"]]), naive17.cluster.averages[["RNA"]][, 2])
data.naive17 <- as.data.frame(as.matrix(naive17.avgClus1))
fwrite(x = data.naive17, file =paste(filedir,"Export_Cluster-17-Avg-Exp_Naive-17_SF-161129-00158",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


naive_18_cells <- subset(x=FullSetMultiRes, cells=h12) #SF100109110236
dim(naive_18_cells)
# [1] 33694  2795
naive18.cluster.averages <- AverageExpression(naive_18_cells)
head(naive18.cluster.averages[["RNA"]][, 1:5])
#                    0           1         12         15        17
# TSPAN6   0.002173536 0.009723919 0.00000000 0.00000000 0.5858101
# TNMD     0.000000000 0.000000000 0.00000000 0.00000000 0.0000000
# DPM1     0.180468041 0.242302825 0.25262843 0.16868474 0.0000000
# SCYL3    0.062193966 0.033476567 0.02256970 0.00000000 0.0000000
# C1orf112 0.086414539 0.106120535 0.05254347 0.00000000 0.1927116
# FGR      0.013715706 0.018800488 0.24564771 0.07944074 0.0000000
naive18.avgClus1 <- cbind(rownames(naive18.cluster.averages[["RNA"]]), naive18.cluster.averages[["RNA"]][, 2])
data.naive18 <- as.data.frame(as.matrix(naive18.avgClus1))
fwrite(x = data.naive18, file =paste(filedir,"Export_Cluster-1-Avg-Exp_Naive-18_SF-100109-110236",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


naive_19_cells <- subset(x=FullSetMultiRes, cells=h13) #SF12092600014
dim(naive_19_cells)
# [1] 33694  2622
naive19.cluster.averages <- AverageExpression(naive_19_cells)
head(naive19.cluster.averages[["RNA"]][, 1:5])
#                    0          1         12 15 2
# TSPAN6   0.003714636 0.02047475 0.00000000  0 0
# TNMD     0.000000000 0.00000000 0.00000000  0 0
# DPM1     0.315135547 0.37305930 0.28813512  0 0
# SCYL3    0.051811905 0.04928334 0.05922499  0 0
# C1orf112 0.071995994 0.07312312 0.00000000  0 0
# FGR      0.014934154 0.02327673 0.00000000  0 0
naive19.avgClus1 <- cbind(rownames(naive19.cluster.averages[["RNA"]]), naive19.cluster.averages[["RNA"]][, 2])
data.naive19 <- as.data.frame(as.matrix(naive19.avgClus1))
fwrite(x = data.naive19, file =paste(filedir,"Export_Cluster-1-Avg-Exp_Naive-19_SF-120926-00014",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


naive_20_cells <- subset(x=FullSetMultiRes, cells=h14) #SF14101000049
dim(naive_20_cells)
# [1] 33694  4221
naive20.cluster.averages <- AverageExpression(naive_20_cells)
head(naive20.cluster.averages[["RNA"]][, 1:5])
#                    0           1         12        15         16
# TSPAN6   0.002143853 0.006495589 0.00000000 0.0000000 1.36956996
# TNMD     0.000000000 0.000000000 0.00000000 0.0000000 0.00000000
# DPM1     0.166450738 0.202632628 0.08437789 0.1053264 0.08502895
# SCYL3    0.031107661 0.033642974 0.07504725 0.1529988 0.00000000
# C1orf112 0.066319591 0.085966989 0.03593001 0.1719424 0.00000000
# FGR      0.029127886 0.017039034 0.24470696 0.8005135 0.00000000
naive20.avgClus1 <- cbind(rownames(naive20.cluster.averages[["RNA"]]), naive20.cluster.averages[["RNA"]][, 2])
data.naive20 <- as.data.frame(as.matrix(naive20.avgClus1))
fwrite(x = data.naive20, file =paste(filedir,"Export_Cluster-1-Avg-Exp_Naive-20_SF-141010-00049",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


naive_21_cells <- subset(x=FullSetMultiRes, cells=h15) #SF12062800475
dim(naive_21_cells)
# [1] 33694  4153
naive21.cluster.averages <- AverageExpression(naive_21_cells)
head(naive21.cluster.averages[["RNA"]][, 1:5])
#                    0           1 11         12 14
# TSPAN6   0.008543777 0.004031616  0 0.00000000  0
# TNMD     0.000000000 0.000000000  0 0.00000000  0
# DPM1     0.152307817 0.180331916  0 0.09486274  0
# SCYL3    0.055325538 0.058630956  0 0.04456391  0
# C1orf112 0.077221826 0.098865195  0 0.06992865  0
# FGR      0.063161151 0.045056662  0 0.84944144  0
naive21.avgClus1 <- cbind(rownames(naive21.cluster.averages[["RNA"]]), naive21.cluster.averages[["RNA"]][, 2])
data.naive21 <- as.data.frame(as.matrix(naive21.avgClus1))
fwrite(x = data.naive21, file =paste(filedir,"Export_Cluster-1-Avg-Exp_Naive-21_SF-120628-00475",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


naive_22_cells <- subset(x=FullSetMultiRes, cells=h16) #SF14072200012
dim(naive_22_cells)
# [1] 33694  4798
naive22.cluster.averages <- AverageExpression(naive_22_cells)
head(naive22.cluster.averages[["RNA"]][, 1:5])
#                   0          1        12 15 2
# TSPAN6   0.02087303 0.04724507 0.0000000  0 0
# TNMD     0.00000000 0.00000000 0.0000000  0 0
# DPM1     0.26778448 0.24353178 0.3482951  0 0
# SCYL3    0.04578396 0.05512252 0.0000000  0 0
# C1orf112 0.06767744 0.11140108 0.0000000  0 0
# FGR      0.07860754 0.08370985 0.6349632  0 0
naive22.avgClus1 <- cbind(rownames(naive22.cluster.averages[["RNA"]]), naive22.cluster.averages[["RNA"]][, 2])
data.naive22 <- as.data.frame(as.matrix(naive22.avgClus1))
fwrite(x = data.naive22, file =paste(filedir,"Export_Cluster-1-Avg-Exp_Naive-22_SF-140722-00012",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


naive_23_cells <- subset(x=FullSetMultiRes, cells=h17) #SF13061200056
dim(naive_23_cells)
# [1] 33694  3794
naive23.cluster.averages <- AverageExpression(naive_23_cells)
head(naive23.cluster.averages[["RNA"]][, 1:5])
#                   0          1         12 13        15
# TSPAN6   0.00347483 0.01450063 0.00000000  0 0.0000000
# TNMD     0.00000000 0.00000000 0.00000000  0 0.0000000
# DPM1     0.17799338 0.23167791 0.10671060  0 0.0000000
# SCYL3    0.05502777 0.06491095 0.10513912  0 0.0000000
# C1orf112 0.07192317 0.10638702 0.08060405  0 0.0000000
# FGR      0.01775478 0.01597328 0.13606869  0 0.4409171
naive23.avgClus1 <- cbind(rownames(naive23.cluster.averages[["RNA"]]), naive23.cluster.averages[["RNA"]][, 2])
data.naive23 <- as.data.frame(as.matrix(naive23.avgClus1))
fwrite(x = data.naive23, file =paste(filedir,"Export_Cluster-1-Avg-Exp_Naive-23_SF-130612-00056",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


naive_24_cells <- subset(x=FullSetMultiRes, cells=h18) #6AD001
dim(naive_24_cells)
# [1] 33694  7108
naive24.cluster.averages <- AverageExpression(naive_24_cells)
head(naive24.cluster.averages[["RNA"]][, 1:5])
#                    0          1         12          14        15
# TSPAN6   0.003340236 0.01983156 0.00000000 0.007861264 0.0000000
# TNMD     0.000000000 0.00000000 0.00000000 0.000000000 0.0000000
# DPM1     0.230228634 0.27907379 0.24273622 0.346899893 0.3863916
# SCYL3    0.055995222 0.05365615 0.03372711 0.265423585 0.0000000
# C1orf112 0.056126181 0.06286921 0.01708328 0.164160406 0.0000000
# FGR      0.025258424 0.01403932 0.17678124 1.889521749 0.0000000
naive24.avgClus1 <- cbind(rownames(naive24.cluster.averages[["RNA"]]), naive24.cluster.averages[["RNA"]][, 2])
data.naive24 <- as.data.frame(as.matrix(naive24.avgClus1))
fwrite(x = data.naive24, file =paste(filedir,"Export_Cluster-1-Avg-Exp_Naive-24_SF-6-AD-001",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


naive_25_cells <- subset(x=FullSetMultiRes, cells=h19) #SF14060200025
dim(naive_25_cells)
# [1] 33694  2586
naive25.cluster.averages <- AverageExpression(naive_25_cells)
head(naive25.cluster.averages[["RNA"]][, 1:5])
#                   0           1 12          13        2
# TSPAN6   0.00000000 0.007976226  0 0.001613415 0.000000
# TNMD     0.00000000 0.000000000  0 0.000000000 0.000000
# DPM1     0.22784850 0.255607558  0 0.202602038 0.000000
# SCYL3    0.05505397 0.098137348  0 0.080058596 0.000000
# C1orf112 0.04430458 0.058474798  0 0.041689381 0.000000
# FGR      0.02341062 0.009974824  0 0.003486599 2.027575
naive25.avgClus1 <- cbind(rownames(naive25.cluster.averages[["RNA"]]), naive25.cluster.averages[["RNA"]][, 2])
data.naive25 <- as.data.frame(as.matrix(naive25.avgClus1))
fwrite(x = data.naive25, file =paste(filedir,"Export_Cluster-1-Avg-Exp_Naive-25_SF-140602-00025",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


hydrea_cells <- subset(x=FullSetMultiRes, cells=h11) #5H001
dim(hydrea_cells)
# [1] 33694  4096
hydrea.cluster.averages <- AverageExpression(hydrea_cells)
head(hydrea.cluster.averages[["RNA"]][, 1:5])
#                     0          1 12 15          17
# TSPAN6   0.0001664432 0.01079884  0  0 0.328791137
# TNMD     0.0000000000 0.00000000  0  0 0.003906087
# DPM1     0.2947256374 0.26704108  0  0 0.401824205
# SCYL3    0.0360492916 0.03867915  0  0 0.101517118
# C1orf112 0.0683596364 0.07677144  0  0 0.018444923
# FGR      0.0499702480 0.01585489  0  0 0.057410905
hydrea.avgClus1 <- cbind(rownames(hydrea.cluster.averages[["RNA"]]), hydrea.cluster.averages[["RNA"]][, 2])
data.hydrea <- as.data.frame(as.matrix(hydrea.avgClus1))
fwrite(x = data.hydrea, file =paste(filedir,"Export_Cluster-1-Avg-Exp_Hydrea_5-H-001",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


normal_1_cells <- subset(x=FullSetMultiRes, cells=n1) #CD34
dim(normal_1_cells)
# [1] 33694  9133
normal1.cluster.averages <- AverageExpression(normal_1_cells)
head(normal1.cluster.averages[["RNA"]][, 1:5])
#                   0         1        12          15 2
# TSPAN6   0.00000000 0.0000000 0.0000000 0.000000000 0
# TNMD     0.00000000 0.0000000 0.0000000 0.000000000 0
# DPM1     0.60828413 0.9126117 0.3398502 0.237526527 0
# SCYL3    0.00000000 0.0000000 0.2332437 0.009378313 0
# C1orf112 0.07024939 0.1926411 0.0000000 0.026006599 0
# FGR      0.00000000 0.0000000 0.7291585 0.433323897 0
normal1.avgClus1 <- cbind(rownames(normal1.cluster.averages[["RNA"]]), normal1.cluster.averages[["RNA"]][, 2])
data.normal1 <- as.data.frame(as.matrix(normal1.avgClus1))
fwrite(x = data.normal1, file =paste(filedir,"Export_Cluster-1-Avg-Exp_Normal-1_Zheng-Pt-CD34",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


normal_2_cells <- subset(x=FullSetMultiRes, cells=n2) #SettyPt1
dim(normal_2_cells)
# [1] 33694 14556
normal2.cluster.averages <- AverageExpression(normal_2_cells)
head(normal2.cluster.averages[["RNA"]][, 1:5])
#                   1        15       17 2        3
# TSPAN6   0.00000000 0.0000000 0.000000 0 0.000000
# TNMD     0.00000000 0.0000000 0.000000 0 0.000000
# DPM1     0.23930036 1.5556473 4.464286 0 1.081577
# SCYL3    0.02049988 0.0000000 0.000000 0 0.000000
# C1orf112 0.03080952 0.0000000 0.000000 0 0.000000
# FGR      0.01205337 0.4516508 0.000000 0 1.153622
normal2.avgClus1 <- cbind(rownames(normal2.cluster.averages[["RNA"]]), normal2.cluster.averages[["RNA"]][, 1])
data.normal2 <- as.data.frame(as.matrix(normal2.avgClus1))
fwrite(x = data.normal2, file =paste(filedir,"Export_Cluster-1-Avg-Exp_Normal-2_Setty-Pt-1",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


normal_3_cells <- subset(x=FullSetMultiRes, cells=n3) #SettyPt2
dim(normal_3_cells)
# [1] 33694  5751
normal3.cluster.averages <- AverageExpression(normal_3_cells)
head(normal3.cluster.averages[["RNA"]][, 1:5])
#                 0           1       12 15       17
# TSPAN6   0.000000 0.006409599  0.00000  0 2.575992
# TNMD     0.000000 0.000000000  0.00000  0 0.000000
# DPM1     0.920556 0.874430602  0.00000  0 0.000000
# SCYL3    0.000000 0.157720247  0.00000  0 0.000000
# C1orf112 0.000000 0.110253957  0.00000  0 0.000000
# FGR      0.000000 0.007725348 12.18027  0 0.000000
normal3.avgClus1 <- cbind(rownames(normal3.cluster.averages[["RNA"]]), normal3.cluster.averages[["RNA"]][, 2])
data.normal3 <- as.data.frame(as.matrix(normal3.avgClus1))
fwrite(x = data.normal3, file =paste(filedir,"Export_Cluster-1-Avg-Exp_Normal-3_Setty-Pt-2",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


normal_4_cells <- subset(x=FullSetMultiRes, cells=n4) #SettyPt3
dim(normal_4_cells)
# [1] 33694  4391
normal4.cluster.averages <- AverageExpression(normal_4_cells)
head(normal4.cluster.averages[["RNA"]][, 1:5])
#          0          1 15 2        3
# TSPAN6   0 0.00000000  0 0 0.000000
# TNMD     0 0.00000000  0 0 0.000000
# DPM1     0 0.29656560  0 0 0.000000
# SCYL3    0 0.22433077  0 0 0.000000
# C1orf112 0 0.19699951  0 0 0.000000
# FGR      0 0.02221492  0 0 2.423655
normal4.avgClus1 <- cbind(rownames(normal4.cluster.averages[["RNA"]]), normal4.cluster.averages[["RNA"]][, 2])
data.normal4 <- as.data.frame(as.matrix(normal4.avgClus1))
fwrite(x = data.normal4, file =paste(filedir,"Export_Cluster-1-Avg-Exp_Normal-4_Setty-Pt-3",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


normal_5_cells <- subset(x=FullSetMultiRes, cells=n5) #HuaPt1
dim(normal_5_cells)
# [1] 33694 10134
normal5.cluster.averages <- AverageExpression(normal_5_cells)
head(normal5.cluster.averages[["RNA"]][, 1:5])
#                  0          1          10 15           2
# TSPAN6   0.0000000 0.00000000 0.013911923  0 0.010392695
# TNMD     0.0000000 0.00000000 0.000000000  0 0.000000000
# DPM1     1.0780459 0.85059788 1.083499907  0 0.710521911
# SCYL3    0.1404095 0.05864827 0.050890306  0 0.050799799
# C1orf112 0.1317131 0.13954329 0.074070533  0 0.036935076
# FGR      0.0000000 0.02261834 0.006205873  0 0.003990728
normal5.avgClus1 <- cbind(rownames(normal5.cluster.averages[["RNA"]]), normal5.cluster.averages[["RNA"]][, 2])
data.normal5 <- as.data.frame(as.matrix(normal5.avgClus1))
fwrite(x = data.normal5, file =paste(filedir,"Export_Cluster-1-Avg-Exp_Normal-5_Hua-Pt-1",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


normal_6_cells <- subset(x=FullSetMultiRes, cells=n6) #HuaPt2
dim(normal_6_cells)
# [1] 33694 10902
normal6.cluster.averages <- AverageExpression(normal_6_cells)
head(normal6.cluster.averages[["RNA"]][, 1:5])
#                 0          1          10        15 17
# TSPAN6   0.000000 0.00000000 0.009394837 0.0000000  0
# TNMD     0.000000 0.00000000 0.000000000 0.0000000  0
# DPM1     1.257729 0.84074001 1.099467932 1.9802161  0
# SCYL3    0.156128 0.05861424 0.053793631 0.0000000  0
# C1orf112 0.000000 0.10982653 0.061912496 0.0000000  0
# FGR      0.000000 0.04823089 0.004208627 0.7916403  0
normal6.avgClus1 <- cbind(rownames(normal6.cluster.averages[["RNA"]]), normal6.cluster.averages[["RNA"]][, 2])
data.normal6 <- as.data.frame(as.matrix(normal6.avgClus1))
fwrite(x = data.normal6, file =paste(filedir,"Export_Cluster-1-Avg-Exp_Normal-6_Hua-Pt-2",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


normal_7_cells <- subset(x=FullSetMultiRes, cells=n7) #HuaPt4
dim(normal_7_cells)
# [1] 33694  4908
normal7.cluster.averages <- AverageExpression(normal_7_cells)
head(normal7.cluster.averages[["RNA"]][, 1:5])
#                  0          1 10           2        4
# TSPAN6   0.0000000 0.06460032  0 0.010845651 0.000000
# TNMD     0.0000000 0.00000000  0 0.000000000 0.000000
# DPM1     2.5839516 1.73731625  0 1.154713049 1.223691
# SCYL3    0.0000000 0.00000000  0 0.063116667 0.000000
# C1orf112 0.3131655 0.00000000  0 0.135814026 0.000000
# FGR      0.0000000 0.00000000  0 0.005145397 0.000000
normal7.avgClus1 <- cbind(rownames(normal7.cluster.averages[["RNA"]]), normal7.cluster.averages[["RNA"]][, 2])
data.normal7 <- as.data.frame(as.matrix(normal7.avgClus1))
fwrite(x = data.normal7, file =paste(filedir,"Export_Cluster-1-Avg-Exp_Normal-7_Hua-Pt-4",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)


BMT_cells <- subset(x=FullSetMultiRes, cells=h5) #4J003
dim(BMT_cells)
# [1] 33694  3584
BMT.cluster.averages <- AverageExpression(BMT_cells)
head(BMT.cluster.averages[["RNA"]][, 1:5])
#                    0           1       12         2          3
# TSPAN6   0.004096494 0.006062676 0.000000 0.0000000 0.06866248
# TNMD     0.000000000 0.000000000 0.000000 0.0000000 0.05489227
# DPM1     0.284685034 0.345852328 0.000000 0.1053474 0.34973367
# SCYL3    0.041859724 0.033978647 1.584284 0.0000000 0.05814583
# C1orf112 0.079465293 0.131484527 0.000000 0.0000000 0.00000000
# FGR      0.021675605 0.006329807 0.000000 0.0000000 0.82314362
BMT.avgClus1 <- cbind(rownames(BMT.cluster.averages[["RNA"]]), BMT.cluster.averages[["RNA"]][, 2])
data.BMT <- as.data.frame(as.matrix(BMT.avgClus1))
fwrite(x = data.BMT, file =paste(filedir,"Export_Cluster-1-Avg-Exp_Normal-8_BMT_4-J-003",date,".csv",sep=""), row.names=TRUE, col.names=TRUE)

#<-------------------------- HERE

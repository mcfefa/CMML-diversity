library(Seurat)
library(dplyr)

rm(list=ls())

#Read in seurat with scores, set directory
dir <- '/Users/Brian/scCode/wu_cellType_data/'
allSeurat <- readRDS(paste0(dir, 'allSeurat_39+8_postStandardPipeline_withHarmony_withSeuratGeneScores_05-20-2021.rds'))

#Get the Wu Score and find cell type fractions for all samples
wuScore <- data.frame("Sample" = allSeurat$orig.ident, "HSC" = allSeurat$wu_HSC, 
                      "MEP" = allSeurat$wu_MEP, "GMP" = allSeurat$wu_GMP, 
                      "ProB" = allSeurat$wu_ProB, "ETP" = allSeurat$wu_ETP)

#Get the GMP score and cluster 2 info in a df, write to csv
clus2_tru_false <- allSeurat$clusterResolution_0.05 == 2
wuGMP_clus2 <- data.frame("Sample" = allSeurat$orig.ident, "GMP" = allSeurat$wu_GMP, "Clus" = clus2_tru_false)
out_gmp_list <- split(wuGMP_clus2, wuGMP_clus2$Clus)
nonClus2.df <- out_gmp_list[["FALSE"]]
clus2.df <- out_gmp_list[["TRUE"]]
write.csv(nonClus2.df, paste0(dir, "nonClus2_wuGMP_score_05-20-2021.csv"))
write.csv(clus2.df, paste0(dir, "clus2_wuGMP_score_05-20-2021.csv"))

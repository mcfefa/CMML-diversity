library(Seurat)
library(Matrix)
library(stringr)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(RColorBrewer)
library(gplots)
library(mitoClone)

rm(list = ls())


dir <- '/Users/4472241/scCode/varTrix_IWG_matchedCorrect/'
list_bam <- c("CMML_1", "CMML_2", "CMML_3", "CMML_9", "CMML_10", "CMML_11", 
          "CMML_12", "CMML_13", "CMML_14", "CMML_16", "CMML_19", "CMML_20",  "CMML_21", 
          "CMML_22", "CMML_23", "CMML_24", "CMML_25", "CMML_26", "CMML_27",  "CMML_29", 
          "CMML_30", "CMML_31", "CMML_32", "CMML_33", "CMML_35", "CMML_36",  "CMML_37", "CMML_39")

allSeurat <- readRDS('/Users/4472241/scCode/redoEntirePipeline_05-05-2021/standard_pipeline_withHarmony/seuratObj_05-11-2021_postStandardPipeline_withHarmony_allSamples_39+8.rds')

for (i in 1:length(list_bam)){
  countCheck <- 0
  
  #Set dir and get file names (only file in dir, use list.files)
  sample <- list_bam[i]
  sampleDir <- paste0(dir, sample, "/")
  snv_matrix_file <- list.files(sampleDir, pattern = ".mtx", full.names = T)
  snps_file <- list.files(sampleDir, pattern = ".txt", full.names = T)
  
  # Read in the sparse genotype matrix
  snv_matrix <- readMM(snv_matrix_file)
  
  # convert the matrix to a dataframe
  snv_matrix <- as.data.frame(as.matrix(t(snv_matrix)))
  
  #read in the cell barcodes output by Cell Ranger
  barcodes <- read.table(paste0(sampleDir, sample, ".tsv"), header = F)
  
  # Construct the final table of SNPs
  snps <- read.table(snps_file, header = F)
  
  #Get row and columns names for snv_matrix
  row.names(snv_matrix) <- barcodes$V1
  colnames(snv_matrix) <- snps$V1
  
  #Convert from 0, 1, 2, 3 to what they actually mean ("No Call", "ref", etc.)
  for (k in 1:dim(snv_matrix)[2]){
    # No reads detected
    snv_matrix[,k] <- str_replace(as.character(snv_matrix[,k]), "0", "No Call")
    # Only ref detected
    snv_matrix[,k] <- str_replace(as.character(snv_matrix[,k]), "1", "ref/ref")
    # Only alt detected
    snv_matrix[,k] <- str_replace(as.character(snv_matrix[,k]), "2", "alt/alt")
    # Both alleles detected
    snv_matrix[,k] <- str_replace(as.character(snv_matrix[,k]), "3", "alt/ref")
  }
  
  #Filter the ones where we don't see mutation in scRNA
  cut.df <- snv_matrix
  for (k in 1:dim(snv_matrix)[2]){
    if("alt/alt" %in% snv_matrix[,k] | "alt/ref" %in% snv_matrix[,k]){
      #print(paste0("Column ", k, ", (", colnames(snv_matrix)[k],  ") has a mutation in scRNA"))
      cut.df[,k] <- snv_matrix[,k]
    }else{
      cut.df[,k] <- 0
    }
    
  }
  
  #Check how many mutations we have compared to how many we expected
  cut2.df <- cut.df[,cut.df[1,] != 0]
  
  #Calculate the VAF for the mutation (count alt/ref as mutated?)
  vaf <- data.frame((colSums(cut2.df == "alt/alt")+colSums(cut2.df == "alt/ref"))/
                      (colSums(cut2.df == "ref/ref")+colSums(cut2.df == "alt/alt")+colSums(cut2.df == "alt/ref")))
  
  write.csv(vaf, paste0(sampleDir, 'vaf_scRNA.csv'))
  
  #Now plot in FDL and analyze mitoclone clones by nuclear mutation
  cols <- c("alt/alt" = "red", "alt/ref" = "blue", 
            "a ref/ref" = "black", 'No Call' = "grey86")
  dir.create(file.path(sampleDir, 'FDL_plots'), showWarnings = FALSE)
  outdir <- paste0(sampleDir, 'FDL_plots/')
  
  cmml_number <- gsub("CMML_", "", list_bam[i])
  
  #Mitoclone load dir
  mitoDir <- '/Users/4472241/scCode/mitoclone/mitoclone_RDS_defaultParams_MINREADS3_MINCELL10_MINFRAC0.05_MINRelative.Other0.99/'
  #Read in the data from the mitoclone analysis
  P_sample <- readRDS(paste0(mitoDir, 'P', cmml_number, '_results_mitoclone.rds'))
  
  #Cluster metaclones
  #CMML_sample_tree <- quick_cluster(P_sample)
  CMML_sample_tree <- clusterMetaclones(P_sample, min.lik = 1)
  
  #Get clonal info and barcodes from mutaCluster output
  cmml_clonalInfo <- data.frame(CMML_sample_tree@mainClone)
  #clone <- data.frame("CB"=row.names(cmml_clonalInfo), "Clone" = NA)
  clone <- apply(CMML_sample_tree@mainClone,1,which.max)
  #for (j in 1:length(cmml_clonalInfo[,1])){
  #  clone[j, "Clone"] <- which.max(cmml_clonalInfo[j,])[[1]]
  #}
  
  #Get barcodes from mitoclone and keep only the actual barcode
  barcodesFromMito <- row.names(cmml_clonalInfo)
  barcodesOnlyFromMito <- sub(".*CB_", "", barcodesFromMito)
  barcodesOnlyMito <- gsub("-1.*","",barcodesOnlyFromMito)
  clone.df <- data.frame("CB" = barcodesOnlyMito, "Clone" = clone)
  
  for (k in 1:dim(cut2.df)[2]){
    # convert so ordering plots grey in background
    cut2.df[,k] <- str_replace(as.character(cut2.df[,k]), "ref/ref", "a ref/ref")
  }
  
  for(k in 1:length(cut2.df[1,])){
    df <- data.frame("Mutated" = sort(cut2.df[,k], decreasing = T))
    row.names(df) <- row.names(cut2.df)
      
    #Configure into palantir gene expression FDL
    paldir <- '/Users/4472241/scCode/runPalantir/outputNew/'
    palantir_embeddings <- read.csv(paste0(paldir, 'cmml', cmml_number,'fdl.csv'))
    #cmml_cellNames <- read.csv(paste0(dir, 'cmml', cmml_number, 'CellNames.csv'))
    barcodes <- palantir_embeddings$X
    barcodes <- sub(".*_", "", barcodes)
    palantir_embeddings$X <- barcodes
    colnames(palantir_embeddings) <- c("CB", "X", "Y")
    
    barcodesFromMutations <- row.names(df)
    barcodesOnlyFromMutations <- gsub("-1.*","",barcodesFromMutations)
    df$CB <- barcodesOnlyFromMutations
    df$Clone <- clone.df$Clone
    
    #Match palantir and mutated barcodes
    combine.df <- data.frame("CB" = NA, "X" = NA, "Y" = NA, "Mutated" = NA, "Clone" = NA, "CheckCB" = NA)
    for (j in 1:length(df$CB)){
      if (df$CB[j] %in% barcodes){
        combine.df[j,"Mutated"] <- df$Mutated[j]
        combine.df[j,"CB"] <- df$CB[j]
        combine.df[j, "Clone"] <- df$Clone[j]
        combine.df[j, "CheckCB"] <- df$CB[j]
        index <- which(palantir_embeddings$CB %in% df$CB[j])
        #For the sample that was combined, only take the barcodes for the first one
        if (cmml_number == "25"){
          index <- index[index < 4537]
        }
        if (length(index)>0){
          combine.df[j,"X"] <- palantir_embeddings[index, "X"]
          combine.df[j,"Y"] <- palantir_embeddings[index,"Y"]
        }
      }
    }
    
    #Remove nas, convert clone to factor in combine.df and plot
    combine.df <- combine.df[!is.na(combine.df$CB),]
    combine.df$Mutated <- as.factor(combine.df$Mutated)
    
    combine.df_removeNoCall <- combine.df[!combine.df$Mutated %in% "No Call",]
    removeNoCall.split <- split(combine.df_removeNoCall, combine.df_removeNoCall$Clone)
    if (length(removeNoCall.split) > 0){
      countCheck <- countCheck + 1
      summary.df <- data.frame("Mutation" = NA, "Clone" = NA, "alt.alt" = NA, 
                               "ref.ref" = NA, "alt.ref" = NA)
      for (j in 1:length(removeNoCall.split)){
        summary.df[j,"alt.alt"] <- sum(as.character(removeNoCall.split[[j]]$Mutated) %in% "alt/alt")
        summary.df[j,"ref.ref"] <- sum(as.character(removeNoCall.split[[j]]$Mutated) %in% "a ref/ref")
        summary.df[j,"alt.ref"] <- sum(as.character(removeNoCall.split[[j]]$Mutated) %in% "alt/ref")
        summary.df[j,"Clone"] <- names(removeNoCall.split)[j]
      }
      summary.df[,"Mutation"] <- colnames(cut2.df)[k]
      
      if (countCheck == 1 ){
        summary.all <- summary.df
      }else{
        summary.all <- rbind(summary.all, summary.df)
      }
    }
    
    
    #pdf(paste0(outdir, 'palantirFDL_mutationColored_CMML_', cmml_number, '_cut2dfColumn_',k,'.pdf'))
    #print(ggplot(combine.df, aes(x=X, y=Y, color = Mutated))+ geom_point(size = .2) +
    #        theme(panel.background = element_rect(fill = "white"))+
    #        scale_color_manual(values=cols)+
    #        ggtitle(paste0(colnames(cut2.df)[k])))
    #dev.off()
    
    
  }
  
  write.csv(summary.all, paste0(sampleDir, "CMML_", cmml_number, "mitoclone_nuclear_comparison.csv"))
}



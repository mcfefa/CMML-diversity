# Run mitoclone at higher resolution with samples from same patient combined
library(mitoClone)
library(deepSNV)
library(GenomicRanges)
library(pheatmap)
library(Seurat)
library(ggplot2)
library(dplyr)
library("reticulate")
setwd('~')
Sys.setenv(RETICULATE_PYTHON = "./anaconda3/envs/mitoclone/bin/python", force = T)
knitr::opts_chunk$set(echo = TRUE)
py_config()

rm(list = ls())

# Set dir to save to and from which we pull countTables
dir <- '/Users/Brian/scCode/mitoclone/'
dirCountTables <- '/Users/Brian/scCode/mitoclone/allSampleCountTables/'
saveDir <- paste0(dir, 'sequential_params_from_02022022_run_03262022/')

# Get mutationCallsFromBlacklist as single core fn
source(paste0(dir, 'mutationCallsFromBlacklist_singleCoreFn.R'))

# Import seurat info
umap_embeddings <- readRDS('~/scCode/mitoclone/umap_embeddings_pre_mito_02022022.rds')

#Run each sample separately
ptm <- proc.time()

# Note which samples are from the same patients
sequential_list <- list("A" = c(1, 32, 39), "B" = c(13, 35), "C" = c(25, 33), "D" = c(27, 37),
                        "E" = c(28, 38), "F" = c(15, 34), "G" = c(12, 36), "H" = c(26, 3))

for (i in 1:length(sequential_list)){

  # Get current patient "name"
  currentPatient <- paste0("patient_", names(sequential_list)[i])
  
  currentSample <- paste0("CMML_", sequential_list[[i]][1])
  
  #Read in the countTables for this patient only (must do this one at a time for memory)
  patient_countTables <- readRDS(paste0(dirCountTables, currentSample, "_countTables_01112022.rds"))
  
  # Read in the next sample(s) count tables and combine into one object "patient_countTables"
  for (j in 2:length(sequential_list[[i]])) {
    currentSample <- paste0("CMML_", sequential_list[[i]][j])
    patient_countTables <- append(patient_countTables, readRDS(paste0(dirCountTables, currentSample, "_countTables_01112022.rds")))
  }
  
  # Run mutationCallsFromBlacklist to filter sites and assign as mutated or not in each cell
  lim.cov <- 20 # Number of reads to classify a cell as covered (default = 20)
  min.af <- 0.1 # Number of reads that need to be mutant for a single cell to be classified as a mutant (default = .2)
  min.af.universal <- min.af # Number of reads that need to be mutant for other cells to be classified as a mutant for universal calc (default = min.af = .2)
  min.num.samples.frac <- 0.01 # minimum fraction of mutant cells required for a variant to be kept. (default = 0.01)
  universal.var.cells.frac <- 0.9 # remove all variants classified as mutant in at least this fraction of cells (default = 0.95)
  max.var.na <- 0.5 # variants must have less than this % NA's (default = 0.5)
  max.cell.na <- 0.75 # cells must have less than this % NA's (default = 0.75)
  CMML <- mutationCallsFromBlacklist(patient_countTables, 
                                     min.af = min.af,
                                     min.af.universal = min.af.universal,
                                     lim.cov = lim.cov,
                                     min.num.samples = min.num.samples.frac * length(patient_countTables), 
                                     universal.var.cells = universal.var.cells.frac * length(patient_countTables), 
                                     binarize = 0.1,
                                     max.var.na = max.var.na, 
                                     max.cell.na = max.cell.na)
  
  # Run gurobi to cluster clones
  object <- CMML
  flag <- TRUE # Catch error
  CMML <- tryCatch(
    expr = {
      muta_cluster(object, cores=8, time = 1000, tempfolder = paste0(getwd(),"/CMML_temp"), force_recalc = T)
    },
    error = function(e){ 
      flag <<- FALSE
    }
  )
  if (!flag) next # Catch cases where no sites are included (should not happen with combined run)
  
  # Cluster and plot
  CMML <- clusterMetaclones(CMML, min.lik =1)
  dev.off()
  pdf(paste0(saveDir, currentPatient, '_plotClones.pdf'))
  print(plotClones(CMML))
  dev.off()
  
  #Get clonal info from mutaCluster output
  cmml_clonalInfo <- data.frame(CMML@cell2clone)
  clone <- data.frame("CB" = row.names(cmml_clonalInfo), "Clone" = NA)
  for (j in 1:length(cmml_clonalInfo[,1])){
    clone[j, "Clone"] <- which.max(cmml_clonalInfo[j,])[[1]]
  }
  
  # Format barcodes
  barcodesFromMito <- clone[,1]
  barcodesOnlyFromMito <- sub("CB_", "", barcodesFromMito)
  barcodesOnlyMito <- gsub("-1.*","",barcodesOnlyFromMito)
  clone$CB <- barcodesOnlyMito
  clone$Sample <- substr(clone$CB, 1, 7)
  
  #Match seurat and mitoclone barcodes
  combine.df <- data.frame("CB" = NA, "X" = NA, "Y" = NA, "Clone" = NA, "Clus" = NA, "Sample" = NA)
  barcodes <- row.names(umap_embeddings)
  for (j in 1:length(clone$CB)){
    if (clone$CB[j] %in% barcodes){
      combine.df[j,"Clone"] <- clone$Clone[j]
      combine.df[j,"CB"] <- clone$CB[j]
      index <- which(row.names(umap_embeddings) %in% clone$CB[j])
      combine.df[j,"X"] <- umap_embeddings[index, "X"]
      combine.df[j,"Y"] <- umap_embeddings[index,"Y"]
      combine.df[j,"Palantir_X"] <- umap_embeddings[index, "palantir_embedding_X"]
      combine.df[j,"Palantir_Y"] <- umap_embeddings[index,"palantir_embedding_Y"]
      combine.df[j, "Clus"] <- as.numeric(as.character(umap_embeddings[index, "Clus"]))
      combine.df[j, "cellType"] <- umap_embeddings[index, "cellType"]
      combine.df[j, "Sample"] <- clone$Sample[j]
    }
  }
  
  # Add binary mutational status to combine.df
  binary_mut_status <- mutate_all(as.data.frame(CMML@ternary)[1:dim(combine.df)[1],], function(x) as.integer(x))
  row.names(binary_mut_status) <- barcodesOnlyMito[1:dim(combine.df)[1]]
  combine_muts <- cbind(combine.df, binary_mut_status)
  
  # Remove na's (cells that were filtered by our QC)
  combine_muts <- combine_muts[!is.na(combine_muts$CB),]
  
  # Plot clone information in palantir embeddings
  pdf(paste0(saveDir, currentPatient, '_palantir_embeddings_mitoclone.pdf'), width = 8, height = 6)
  print(ggplot(combine_muts) + geom_point(aes(x = Palantir_X, y = Palantir_Y, color = as.factor(Clone)), size = .4))
  dev.off()
  
  # Save 
  mito_object_dir <- '~/scCode/mitoclone/sequential_params_from_02022022_run_03262022/'
  saveRDS(CMML, paste0(mito_object_dir, currentPatient, '_clustered_mitoclone_object_03262022.rds'))
  saveRDS(combine_muts, paste0(mito_object_dir, currentPatient, '_combine_muts_df.rds'))
  rm(CMML)
}
                                  
proc.time() - ptm

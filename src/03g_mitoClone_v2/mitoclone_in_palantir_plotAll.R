if (paste0(mito_object_dir, '/cmml_', i, '_clustered_mitoclone_object_02022022.rds') %in% mito_objects) {
    CMML <- readRDS(paste0(mito_object_dir, 'cmml_', i, '_clustered_mitoclone_object_02022022.rds'))
  } else {
    next
  }
  
  #Get clonal info from mutaCluster output
  cmml_clonalInfo <- data.frame(CMML@cell2clone)
  clone <- data.frame("CB" = row.names(cmml_clonalInfo), "Clone" = NA)
  for (j in 1:length(cmml_clonalInfo[,1])){
    clone[j, "Clone"] <- which.max(cmml_clonalInfo[j,])[[1]]
  }
  
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
  
  #pdf(paste0(saveDir, 'cmml', i, '_palantir_embeddings_mitoclone.pdf'), width = 6, height = 4.5)
  p1 <- ggplot(combine_muts) + geom_jitter(aes(x = Palantir_X, y = Palantir_Y, color = as.factor(Clone)), size = .6) +
          theme_bw()+ theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(face = "bold", size = 15),
                            axis.ticks = element_blank(), axis.text = element_blank(),
                            axis.title = element_blank())+ggtitle(paste0(currentSample))+
    guides(color="none")

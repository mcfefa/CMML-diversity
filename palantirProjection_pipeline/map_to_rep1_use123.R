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

#*********************REP SENSITIVE*******************************************
filedir = "/Users/4472241/scCode/palantirEmbeddings/mapToRep1/"

#Create the dataframe to store poor integration data
storeIntegration_df <- data.frame("Sample"=character(), "Cluster"=integer(), "ClusterComposition"=double())

#***********Load the rep3 data from csv
all_content <- readLines("/Users/4472241/scCode/rep3Info/csvRawFilteredDataRep3.csv")
dataRep3 <- read.csv(textConnection(all_content), header = FALSE, stringsAsFactors = FALSE)
geneNamesRep3 <- read.csv("/Users/4472241/scCode/rep3Info/var.csv")
cellNamesExtraRep3 <- read.csv("/Users/4472241/scCode/rep3Info/obs.csv")
cellNamesRep3.df <- data.frame(cellNamesExtraRep3)
cellNamesRep3 <- data.frame(cellNamesRep3.df[,1])
row.names(dataRep3) <- cellNamesRep3[,1]
colnames(dataRep3) <- geneNamesRep3[,1]
dataRep3 <- t(dataRep3)

#Create Seurat Object
rep3 <- CreateSeuratObject(dataRep3)
rep3[["Identity"]] <- "rep3"
rep3[["Integration"]] <- "Reference"

#Import Rep3 Embeddings
rep3Embeddings <- read.table('/Users/4472241/scCode/rep3Info/rep3Embeddings.csv', sep = "" , header = F ,
                             na.strings ="", stringsAsFactors= F)

#***********Load the rep2 data from csv
all_content <- readLines("/Users/4472241/scCode/rep2Info/csvRawFilteredDataRep2.csv")
dataRep2 <- read.csv(textConnection(all_content), header = FALSE, stringsAsFactors = FALSE)
geneNamesRep2 <- read.csv("/Users/4472241/scCode/rep2Info/var.csv")
cellNamesExtraRep2 <- read.csv("/Users/4472241/scCode/rep2Info/obs.csv")
cellNamesRep2.df <- data.frame(cellNamesExtraRep2)
cellNamesRep2 <- data.frame(cellNamesRep2.df[,1])
row.names(dataRep2) <- cellNamesRep2[,1]
colnames(dataRep2) <- geneNamesRep2[,1]
dataRep2 <- t(dataRep2)

#Create Seurat Object
rep2 <- CreateSeuratObject(dataRep2)
rep2[["Identity"]] <- "rep2"
rep2[["Integration"]] <- "Reference"

#Import Rep2 Embeddings
rep2Embeddings <- read.table('/Users/4472241/scCode/rep2Info/rep2Embeddings.csv', sep = "," , header = F ,
                             na.strings ="", stringsAsFactors= F)

#***********Load the rep1 data from csv
all_content <- readLines("/Users/4472241/scCode/rep1Info/csvRawFilteredDataRep1.csv")
dataRep1 <- read.csv(textConnection(all_content), header = FALSE, stringsAsFactors = FALSE)
geneNamesRep1 <- read.csv("/Users/4472241/scCode/rep1Info/var.csv")
cellNamesExtraRep1 <- read.csv("/Users/4472241/scCode/rep1Info/obs.csv")
cellNamesRep1.df <- data.frame(cellNamesExtraRep1)
cellNamesRep1 <- data.frame(cellNamesRep1.df[,1])
row.names(dataRep1) <- cellNamesRep1[,1]
colnames(dataRep1) <- geneNamesRep1[,1]
dataRep1 <- t(dataRep1)

#Create Seurat Object
rep1 <- CreateSeuratObject(dataRep1)
rep1[["Identity"]] <- "rep1"
rep1[["Integration"]] <- "Reference"

#Import Rep1 Embeddings
rep1Embeddings <- read.table('/Users/4472241/scCode/rep1Info/rep1Embeddings.csv', sep = "," , header = F ,
                             na.strings ="", stringsAsFactors= F)

#Merge reps 1 2 and 3
reps123 <- merge(rep1, y = c(rep2,rep3))

#Normalize, scale, and run pca
reps123 <- NormalizeData(reps123)
reps123 <- FindVariableFeatures(reps123, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(reps123)
reps123 <- ScaleData(reps123, features = VariableFeatures(reps123))
reps123 <- RunPCA(reps123, features = VariableFeatures(object = reps123))

#Harmonize and run umap, plotting before and after
pdf(paste0(filedir,'pcaBeforeHarmony.pdf'), width = 8, height = 6)
print(DimPlot(reps123, reduction = 'pca', group.by = "Identity"))
dev.off()
#Harmony
reps123 <- RunHarmony(reps123, group.by.vars = c("Identity"), theta = 1, max.iter.harmony = 20)
pdf(paste0(filedir,'pcaAfterHarmony.pdf'), width = 8, height = 6)
print(DimPlot(reps123, reduction = 'harmony', group.by = "Identity"))
dev.off()

#Run UMAP and then plot
reps123 <- RunUMAP(reps123, reduction = 'harmony', dims = 1:10)
pdf(paste0(filedir,'umapAfterHarmony.pdf'), width = 8, height = 6)
print(DimPlot(reps123, reduction = "umap", group.by = "Identity"))
dev.off()

#Map from reps 1 and 2 to 3
##Match Neighbors
#Find harmony embeddings of the query (cmml) and reference (rep3)
#*******************************REP SENSITIVE**********************************
dimBreak <- dim(rep1)[[2]]
harmonyEmbeddings_query <- reps123@reductions[["harmony"]]@cell.embeddings[(1+dimBreak):dim(reps123)[[2]],]
harmonyEmbeddings_ref <- reps123@reductions[["harmony"]]@cell.embeddings[1:dimBreak,]

#Find nearest neighbor of each query in reference set
repsQuery_to_Ref_embedIndex <- queryKNN(harmonyEmbeddings_ref, harmonyEmbeddings_query, k=1, BNPARAM=KmknnParam())

#Create function to copy the embeddings
copy2Query <- function(refEmbeddings, queryIndex)
{
  lengthQuery <- length(queryIndex)
  output <- matrix(data=NA,nrow=lengthQuery,ncol=2)
  for (i in 1:lengthQuery){
    output[i,1] <- refEmbeddings[queryIndex[i],1]
    output[i,2] <- refEmbeddings[queryIndex[i],2]
  }
  return(output)
}

#Run function to map query to palantir embeddings and convert to dataframe
#***********************************REP SENSITIVE******************************
RefEmbeddings <- rep1Embeddings
query_palantir_embeddings <- copy2Query(RefEmbeddings, repsQuery_to_Ref_embedIndex[["index"]])
dfrepsQuery <- data.frame("X" = query_palantir_embeddings[,1], "Y" = query_palantir_embeddings[,2])
#****************************REP SENSITIVE*************************************
dfrepRef <- data.frame("X" = rep1Embeddings[,1], "Y" = rep1Embeddings[,2])
#***************************REP SENSITIVE**************************************
dfreps123 <- rbind(dfrepRef, dfrepsQuery)

#Plot reps 1 and 2 in palantir rep 3 embedding space and visualize density
pdf(paste0(filedir,'2PalantirRepsinThirdembeddingsBinHex.pdf'),width=8, height=5)
print(ggplot(dfrepsQuery,aes(x=X,y=Y)) + stat_binhex(aes(fill=log10(..density..)),binwidth = 1) + 
        scale_fill_gradientn(colours=c("white","cyan", "blue", "black"),name = "Fraction in Hexagon",na.value=NA))
dev.off()

#Add to seurat object metadata
reps123@meta.data[["Palantir Embeddings"]] <- dfreps123

#Import the branch probabilities and differentiation potential of all three reps
###Rep 3
branchProbs_rep3 <- read.csv('/Users/4472241/scCode/rep3Info/rep3branchProbs.csv', sep = ' ', header = F)
diffPot_rep3 <- read.csv('/Users/4472241/scCode/rep3Info/rep3DiffPot.csv', header = F)
cellStates_rep3 <- read.csv('/Users/4472241/scCode/rep3Info/uns/palantir_branch_probs_cell_types.csv', header = F)
row.names(branchProbs_rep3) <- cellNamesRep3[,1]
row.names(diffPot_rep3) <- cellNamesRep3[,1]
colnames(branchProbs_rep3) <- cellStates_rep3[,1]
colnames(diffPot_rep3) <- "Diff. Potential"

###Rep2
branchProbs_rep2 <- read.csv('/Users/4472241/scCode/rep2Info/rep2branchProbs.csv', sep = ' ', header = F)
diffPot_rep2 <- read.csv('/Users/4472241/scCode/rep2Info/rep2DiffPot.csv', header = F)
cellStates_rep2 <- read.csv('/Users/4472241/scCode/rep2Info/rep2BranchProbsCellTypes.csv', header = F)
row.names(branchProbs_rep2) <- cellNamesRep2[,1]
row.names(diffPot_rep2) <- cellNamesRep2[,1]
colnames(branchProbs_rep2) <- cellStates_rep2[,1]
colnames(diffPot_rep2) <- "Diff. Potential"

###Rep 1
branchProbs_rep1 <- read.csv('/Users/4472241/scCode/rep1Info/rep1branchProbs.csv', sep = ' ', header = F)
diffPot_rep1 <- read.csv('/Users/4472241/scCode/rep1Info/rep1DiffPot.csv', header = F)
cellStates_rep1 <- read.csv('/Users/4472241/scCode/rep1Info/rep1BranchProbsCellTypes.csv', header = F)
row.names(branchProbs_rep1) <- cellNamesRep1[,1]
row.names(diffPot_rep1) <- cellNamesRep1[,1]
colnames(branchProbs_rep1) <- cellStates_rep1[,1]
colnames(diffPot_rep1) <- "Diff. Potential"

#set reference branchProb and diffPot dataframes
branchProbs_reps123 <- rbind(branchProbs_rep1, branchProbs_rep2, branchProbs_rep3)
diffPot_reps123 <- rbind(diffPot_rep1, diffPot_rep2, diffPot_rep3)

#Save to cut down time when re-running
#saveRDS(reps123, paste0(filedir,'reps123Integrated.rds'))
#saveRDS(branchProbs_reps123, paste0(filedir,'branchProbsAll3.rds'))
#saveRDS(diffPot_reps123, paste0(filedir, 'diffPotential.rds'))
reps123 <- readRDS(paste0(filedir,'reps123Integrated.rds'))
branchProbs_reps123 <- readRDS(paste0(filedir,'branchProbsAll3.rds'))
diffPot_reps123 <- readRDS(paste0(filedir, 'diffPotential.rds'))

rdsFileList <- c('CMML1_LTB3966_2020-10-15', 'CMML2_LTB4121_2020-10-15', 'CMML3_LTB5109_2020-10-15', 'CMML4_4-J-003_postBMT_2020-10-15',
                 'CMML5_4-K-001_HMA_2020-10-15', 'CMML6_4-Q-001_HMA_2020-10-15', 'CMML7_5-E-001_HMA_2020-10-15',
                 'CMML8_5-H-001_2020-10-15', 'CMML9_SF-100109-106293_2020-10-15', 'CMML10_SF-100109-111451_2020-10-15',
                 'CMML11_SF-100109-110236_2020-10-15', 'CMML12_SF-140401-00158_2020-10-15', 'CMML13_SF-140602-00025_2020-10-15',
                 'CMML14_SF-120628-00475_2020-10-15', 'CMML15_SF-140722-00012_2020-10-15', 'CMML16_SF-130612-00056_2020-10-15',
                 'CMML17_4-S-001_HMA_2020-10-15','CMML18_2-V-001_2020-10-15', 'CMML19_SF-141010-00049_2020-10-15',
                 'CMML20_SF-161129-00158_2020-10-15', 'CMML21_6-AE-001_2020-10-15', 'CMML22_6-AC-001_2020-10-15',
                 'CMML24_SF-100109-101914_2020-10-15','CMML25_SF-120425-00035_2020-10-15', 'CMML23+26_6-AD-001_2020-10-15',
                 'CMML27_SF-120926-00014_2020-10-15','CMML28_SF-140318-00065_2020-10-15', 'CMML29_SF-140507-00419_2020-10-15',
                 'CMML30_SF-160268-00045_2020-10-15', 'CMML31_SF-160722-00003_2020-10-15', 'CMML32_SF-161123-00029_2020-10-15',
                 'CMML33_SF-130328-00016_2020-10-15', 'CMML34_SF-141104-00108_2020-10-15', 'CMML35_SF-141114-00033_2020-10-15',
                 'CMML36_SF-140925-00135_2020-10-15', 'CMML37_SF-140613-00036_2020-10-15', 'CMML38_SF-140804-00065_2020-10-15',
                 'CMML39_SF-150102-00008_2020-10-15', 'CMML40_SF-130709-00171_2020-10-15')

for (i in 1:length(rdsFileList)){
  #Import the file
  cmml <- readRDS(paste0("/Users/4472241/scCode/",rdsFileList[i],".rds"))
  
  #Find mito percentage in each cell
  mito.genes <- grep(pattern = "^MT-", cmml@assays$RNA@counts@Dimnames[[1]], value = TRUE)
  percent.mito <- Matrix::colSums(x=GetAssayData(object=cmml, slot='counts')[mito.genes, ])/Matrix::colSums(x=GetAssayData(object=cmml, slot='counts'))
  cmml[['percent.mito']] <- percent.mito
  
  #Set filter min/max
  mean(cmml@meta.data$nFeature_RNA, na.rm=TRUE) + 2*sd(cmml@meta.data$nFeature_RNA, na.rm=TRUE) # [1] 3326.408
  nFeatUpperN <- mean(cmml@meta.data$nFeature_RNA, na.rm=TRUE) + 2*sd(cmml@meta.data$nFeature_RNA, na.rm=TRUE) 
  nFeatLowerN <- 200 #Unnecessary because we set min.genes to 315 anyway in a few lines
  #To match setty, set mito to 0.2 (instead of 0.25 which we use in our CMML pipeline)
  perMitoUpperN <- 0.20
  
  #filter out those that don't make the cut
  cmml <- subset(x=cmml, subset=nFeature_RNA > nFeatLowerN & percent.mito < perMitoUpperN)
  
  #Get cmml data
  cmml.data <- cmml@assays[["RNA"]]@data
  cmml.data <- data.frame(cmml.data)
  #To match setty paper, remove cells with less than 1000 molecules
  cmml.data <- cmml.data[,colSums(cmml.data)>1000]
  #To match setty paper, remove "low-complexity" cells mapping to less than 315 genes (Why 315???)
  cmml <- CreateSeuratObject(cmml.data, min.genes = 315)
  cmml[["Identity"]] <- paste0("cmml",i)
  cmml[["Integration"]] <- "Query"
  
  # Match Genes from reps123 to cmml
  a <- cmml@assays$RNA@counts@Dimnames[[1]] 
  b <- reps123@assays$RNA@counts@Dimnames[[1]]  
  c <- intersect(a,b) 
  cmml <- cmml[c,]
  reps123 <- reps123[c,]
  
  #Normalize before merging
  cmml <- NormalizeData(cmml)
  
  #Combine seurat objects
  cmmlreps123 <- merge(cmml, reps123, merge.data = TRUE)
  
  #scale together
  cmmlreps123 <- FindVariableFeatures(cmmlreps123)
  cmmlreps123 <- ScaleData(cmmlreps123, features = VariableFeatures(cmmlreps123))
  cmmlreps123 <- RunPCA(cmmlreps123, features = VariableFeatures(cmmlreps123))
  
  pdf(paste0(filedir,"PCA_preHarmony_cmml",i,".pdf"), width=8, height=5)
  print(DimPlot(cmmlreps123, reduction = "pca", group.by = "Integration", order= c("Query", 'Reference')))
  dev.off()
  
  cmmlreps123 <- RunHarmony(cmmlreps123, reduction = "pca", group.by.vars = "Integration", 
                            theta = 2, max.iter.harmony = 20)
  
  pdf(paste0(filedir,"PCA_postHarmony_cmml",i,".pdf"), width=8, height=5)
  print(DimPlot(cmmlreps123, reduction = "harmony", group.by = "Integration", order = c("Query", "Reference")))
  dev.off()
  
  #Compute SNN
  cmmlreps123 <- FindNeighbors(cmmlreps123, reduction = "harmony", dims = 1:30)
  
  #Cluster
  cmmlreps123 <- FindClusters(cmmlreps123, resolution = 0.1)
  
  #convert cluster info to df
  cluster_df <- data.frame(cmmlreps123@meta.data[["seurat_clusters"]])
  
  #Split the df into the cmml sample and the reps123
  lengthcmml <- dim(cmml)[[2]]
  cluster_cmml <- cluster_df[1:lengthcmml,]
  cluster_reps123 <- cluster_df[1+lengthcmml:dim(cluster_df)[[1]],]
  
  #Find the sum of counts in each df
  cmml_cluster_counts <- table(cluster_cmml)
  reps123_cluster_counts <- table(cluster_reps123)
  counter <- 0
  #Record the percentage in each cluster
  for (j in 1:length(table(cluster_df))){
    
    sum <- table(cluster_df)[[j]]
    fractionCmml <- cmml_cluster_counts[[j]]/sum
    
    length_df <- dim(storeIntegration_df)[[1]]
    storeIntegration_df[[length_df+1,1]] <- rdsFileList[[i]]
    storeIntegration_df[[length_df+1,2]] <- j-1
    storeIntegration_df[[length_df+1,3]] <- fractionCmml
  }
  
  #UMAP
  cmmlreps123 <- RunUMAP(cmmlreps123, dims = 1:30, reduction = 'harmony')
  
  #Save UMAP plot
  pdf(paste0(filedir,"UMAP_postHarmony_cmml",i,".pdf"),width=8, height=5)
  p1 <- DimPlot(cmmlreps123, reduction = "umap", group.by = 'Identity', order = c(paste0("cmml",i),"rep1","rep2", "rep3"))
  p2 <- DimPlot(cmmlreps123, reduction = "umap", group.by = 'seurat_clusters')
  p3 <- DimPlot(cmmlreps123, reduction = 'umap', group.by = "Integration", order = c("Query", "Reference"))
  print(p1)
  print(p2)
  print(p3)
  dev.off()
  
  ######Match Neighbors
  #Create function to copy the diffPotential and branchProbs
  copyBranchDiffPot <- function(refDiffPot, refBranchProb, knnData)
  {
    
    knn_index <- knnData[["index"]]
    knn_distance <- knnData[["distance"]]
    
    lengthQuery <- dim(knn_index)[[1]]
    
    output <- data.frame("DiffPot" = double(), "pDC" = double(), "CLP" = double(),
                         "Mono" = double(), "Ery" = double(), "Mega" = double(), "cDC" = double())
    
    
    for (j in 1:lengthQuery){
      
      distance_sum <- sum(1/(knn_distance[j,]))
      
      #Weight by inverse distance
      weightedDiffPot <- sum(1/(knn_distance[j,])*refDiffPot[knn_index[j,],])/distance_sum
      weightedBranchProbs <- colSums(1/(knn_distance[j,])*refBranchProb[knn_index[j,],])/distance_sum
      
      #Assign to each column
      output[j,"DiffPot"] <- weightedDiffPot
      output[j,"Mono"] <- weightedBranchProbs[[1]]
      output[j,"Ery"] <- weightedBranchProbs[[2]]
      output[j,"pDC"] <- weightedBranchProbs[[3]]
      output[j,"Mega"] <- weightedBranchProbs[[4]]
      output[j,"CLP"] <- weightedBranchProbs[[5]]
      output[j,"cDC"] <- weightedBranchProbs[[6]]
    }
    return(output)
  }
  
  #Find harmony embeddings of the query (cmml) and reference (reps123)
  harmonyEmbeddings_Query <- cmmlreps123@reductions[["harmony"]]@cell.embeddings[1:(dim(cmml)[[2]]),]
  harmonyEmbeddings_Reference <- cmmlreps123@reductions[["harmony"]]@cell.embeddings[(dim(cmml)[[2]]+1):dim(cmmlreps123)[[2]],]
  
  #Find nearest neighbor of each query in reference set
  query2ref_embedIndex <- queryKNN(harmonyEmbeddings_Reference, harmonyEmbeddings_Query, k=1, BNPARAM=KmknnParam())
  #Run function (defined previously) to map query to palantir embeddings and convert to dataframe
  query_palantir_embeddings <- copy2Query(dfreps123, query2ref_embedIndex[["index"]])
  dfcmml <- data.frame(query_palantir_embeddings)
  
  #Find k nearest neighbors for diff pot and branch prob calculation
  query2ref_diffPot <- queryKNN(harmonyEmbeddings_Reference, harmonyEmbeddings_Query, k=1, BNPARAM=KmknnParam())
  #Run separate function to average k neighbors to determine branch prob. and diff pot.
  cmml_branchDiffProb_df <- copyBranchDiffPot(diffPot_reps123, branchProbs_reps123, query2ref_diffPot)
  
  #Plot cmml in palantir embedding space and visualize density
  pdf(paste0(filedir,"cmml",i,'embeddingsBinHex.pdf'),width=8, height=5)
  print(ggplot(dfcmml,aes(x=X1,y=X2)) + stat_binhex(aes(fill=log10(..density..)),binwidth = 1) + 
          scale_fill_gradientn(colours=c("white","cyan", "blue", "black"),name = "Log Fraction in Hexagon",na.value=NA)+
          theme(panel.background = element_rect(fill = "gray")))
  dev.off()
  
  #Separate cmml and rep3 data
  seurat.list <- SplitObject(cmmlreps123, split.by = "Identity")
  
  #Take only the cmml
  cmmlOnly <- seurat.list[[1]]
  
  #Plot and save the diffPotential
  pdf(paste0(filedir,"cmml",i,'_diffPotential.pdf'),width=8, height=5)
  print(ggplot(dfcmml,aes(x=X1,y=X2, color = cmml_branchDiffProb_df[['DiffPot']])) + 
          geom_point(size = 0.7)+scale_color_gradient2(midpoint=0.5,low="black",mid = "purple",high="yellow" , name = "Diff. Pot.")+
          theme(panel.background = element_rect(fill = "gray")))
  dev.off()
  
  
  #Make plots for branch prob
  p1 <- ggplot(dfcmml,aes(x=X1,y=X2, color = cmml_branchDiffProb_df[['Ery']])) + 
    geom_point(size = 0.7)+scale_color_gradient2(midpoint=0.5,limits = c(0, 1),low="white", mid = "cyan",high="black",name = "Ery Prob")+
    theme(panel.background = element_rect(fill = "gray"))
  p2 <- ggplot(dfcmml,aes(x=X1,y=X2, color = cmml_branchDiffProb_df[['Mega']])) + 
    geom_point(size = 0.7)+scale_color_gradient2(midpoint=0.5,limits = c(0, 1),low="white", mid = "cyan",high="black",name = "Mega Prob")+
    theme(panel.background = element_rect(fill = "gray"))
  p3 <- ggplot(dfcmml,aes(x=X1,y=X2, color = cmml_branchDiffProb_df[['Mono']])) + 
    geom_point(size = 0.7)+scale_color_gradient2(midpoint=0.5,limits = c(0, 1),low="white", mid = "cyan",high="black",name = "Mono Prob")+
    theme(panel.background = element_rect(fill = "gray"))
  p4 <- ggplot(dfcmml,aes(x=X1,y=X2, color = cmml_branchDiffProb_df[['pDC']])) + 
    geom_point(size = 0.7)+scale_color_gradient2(midpoint=0.5,limits = c(0, 1),low="white", mid = "cyan",high="black",name = "pDC Prob")+
    theme(panel.background = element_rect(fill = "gray"))
  p5 <- ggplot(dfcmml,aes(x=X1,y=X2, color = cmml_branchDiffProb_df[['cDC']])) + 
    geom_point(size = 0.7)+scale_color_gradient2(midpoint=0.5,limits = c(0, 1),low="white", mid = "cyan",high="black",name = "cDC Prob")+
    theme(panel.background = element_rect(fill = "gray"))
  p6 <- ggplot(dfcmml,aes(x=X1,y=X2, color = cmml_branchDiffProb_df[['CLP']])) + 
    geom_point(size = 0.7)+scale_color_gradient2(midpoint=0.5,limits = c(0, 1),low="white", mid = "cyan",high="black",name = "CLP Prob")+
    theme(panel.background = element_rect(fill = "gray"))
  
  #Save plots for branch Prob
  pdf(paste0(filedir,"cmml",i,'_branchProbs.pdf'),width=8, height=5)
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  print(p5)
  print(p6)
  dev.off()
  
  
  #Save normalized data, gene names, cell names, and palantir embeddings
  write.csv(cmmlOnly@assays[["RNA"]]@data, paste0(filedir,'cmml', i,"normalizedData.csv"))
  write.csv(cmmlOnly@assays[["RNA"]]@counts@Dimnames[[1]], paste0(filedir,'cmml', i,'GeneNames.csv'))
  write.csv(cmmlOnly@assays[["RNA"]]@counts@Dimnames[[2]], paste0(filedir,'cmml', i,'CellNames.csv'))
  write.csv(cmmlOnly@reductions[['harmony']]@cell.embeddings, paste0(filedir,'cmml', i,'harmonyEmbeddings.csv'))
  write.csv(cmmlOnly@reductions[['umap']]@cell.embeddings, paste0(filedir,'cmml', i,'UMAPEmbeddings.csv'))
  write.csv(dfcmml, paste0(filedir,'cmml', i,'PalantirEmbeddings.csv'))
  write.csv(cmml_branchDiffProb_df, paste0(filedir,'cmml',i,'branchProbs_diffPot.csv'))
  #saveRDS(cmmlOnly, paste0(filedir,'cmml',i,'.rds'))
  
  rm(cmml, cmmlOnly, cmml_branchDiffProb_df, dfcmml)
  
}

dfrep1 <- data.frame(rep1Embeddings)
pdf(paste0(filedir,"mainRepEmbeddingsBinHex.pdf"), width=8, height=5)
print(ggplot(dfrep1,aes(x=V1,y=V2)) + stat_binhex(aes(fill=log10(..density..)),binwidth = 1) + 
        scale_fill_gradientn(colours=c("white","cyan", "blue", "black"),name = "Fraction in Hexagon",na.value=NA))
dev.off()

#Export the dataframe containing cluster info for each sample/cluster to csv
write.csv(storeIntegration_df, paste0(filedir,'clusterPercentages_maptorep1.csv'))

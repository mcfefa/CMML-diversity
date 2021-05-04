library(Seurat)
library(Matrix)
library(stringr)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(RColorBrewer)
library(gplots)

rm(list = ls())

dir <- '/Users/4472241/scCode/varTrix_CMML_20/'
# Read in the sparse genotype matrix
snv_matrix <- readMM(paste0(dir, "outputMatrix.mtx"))

# convert the matrix to a dataframe
snv_matrix <- as.data.frame(as.matrix(t(snv_matrix)))

#read in the cell barcodes output by Cell Ranger
barcodes <- read.table(paste0(dir, "CMML_20.tsv"), header = F)

# read in SNV loci
# Should be constructed a single column. For example

# chr1:1234-1235
# chr2:2345-2346

# Construct the final table to add to the Seurat object
snps <- read.table(paste0(dir, "SNV.loci.withID.txt"), header = F)

row.names(snv_matrix) <- barcodes$V1

colnames(snv_matrix) <- snps$V1

x <- colSums(snv_matrix>2)

saveRDS(snv_matrix, paste0(dir, 'mutationMatrix_CMML_20_singleCellRes.rds'))

for (i in 1:dim(snv_matrix)[2]){
  # No reads detected
  snv_matrix[,i] <- str_replace(as.character(snv_matrix[,i]), "0", "No Call")
  # Only ref detected
  snv_matrix[,i] <- str_replace(as.character(snv_matrix[,i]), "1", "ref/ref")
  # Only alt detected
  snv_matrix[,i] <- str_replace(as.character(snv_matrix[,i]), "2", "alt/alt")
  # Both alleles detected
  snv_matrix[,i] <- str_replace(as.character(snv_matrix[,i]), "3", "alt/ref")
}

#Filter the ones where we don't see mutation in scRNA
cut.df <- snv_matrix
for (i in 1:dim(snv_matrix)[2]){
  if("alt/alt" %in% snv_matrix[,i] | "alt/ref" %in% snv_matrix[,i]){
    print(paste0("Column ", i, " has a mutation in scRNA"))
    cut.df[,i] <- snv_matrix[,i]
  }else{
    cut.df[,i] <- 0
  }
  
}

#Check how many mutations we have compared to how many we expected
cut2.df <- cut.df[,cut.df[1,] != 0]

#Calculate VAF
vaf <- data.frame((colSums(cut2.df == "alt/alt")+colSums(cut2.df == "alt/ref"))/
                    (colSums(cut2.df == "ref/ref")+colSums(cut2.df == "alt/alt")+colSums(cut2.df == "alt/ref")))

write.csv(vaf, '/Users/4472241/scCode/vaf_scRNA_CMML_20.csv')


#plot using palantir FDLs
dir <- '/Users/4472241/scCode/mitoclone/mitoclone_RDS_defaultParams_MINREADS3_MINCELL10_MINFRAC0.05_MINRelative.Other0.99/'
outdir <- '/Users/4472241/scCode/varTrix_CMML_20/plotMutationsInFDL/'

cols <- c("alt/alt" = "red", "alt/ref" = "blue", 
          "aref/ref" = "black", 'No Call' = "grey86")

for (i in 1:dim(cut2.df)[2]){
  # convert so ordering plots grey in background
  cut2.df[,i] <- str_replace(as.character(cut2.df[,i]), "ref/ref", "aref/ref")
}

for(k in 1:length(cut2.df[1,])){
  tet2 <- data.frame("Mutated" = sort(cut2.df[,k], decreasing = T))
  row.names(tet2) <- row.names(cut2.df)
  i <- 20
  
  #Configure into palantir gene expression FDL
  paldir <- '/Users/4472241/scCode/runPalantir/outputNew/'
  palantir_embeddings <- read.csv(paste0(paldir, 'cmml', i,'fdl.csv'))
  #cmml_cellNames <- read.csv(paste0(dir, 'cmml', i, 'CellNames.csv'))
  barcodes <- palantir_embeddings$X
  barcodes <- sub(".*_", "", barcodes)
  palantir_embeddings$X <- barcodes
  colnames(palantir_embeddings) <- c("CB", "X", "Y")
  
  barcodesFromMutations <- row.names(tet2)
  barcodesOnlyFromMutations <- gsub("-1.*","",barcodesFromMutations)
  tet2$CB <- barcodesOnlyFromMutations
  
  #Match palantir and mitoclone barcodes
  combine.df <- data.frame("CB" = NA, "X" = NA, "Y" = NA, "Mutation" = NA)
  for (j in 1:length(tet2$CB)){
    if (tet2$CB[j] %in% barcodes){
      combine.df[j,"Mutated"] <- tet2$Mutated[j]
      combine.df[j,"CB"] <- tet2$CB[j]
      index <- which(palantir_embeddings$CB %in% tet2$CB[j])
      combine.df[j,"X"] <- palantir_embeddings[index, "X"]
      combine.df[j,"Y"] <- palantir_embeddings[index,"Y"]
    }
  }
  
  #Remove nas, convert clone to factor in combine.df and plot
  combine.df <- combine.df[!is.na(combine.df$CB),]
  combine.df$Mutated <- as.factor(combine.df$Mutated)
  
  pdf(paste0(outdir, 'palantirFDL_mutationColored_CMML_', i, '_cut2dfColumn_',k,'.pdf'))
  print(ggplot(combine.df, aes(x=X, y=Y, color = Mutated))+ geom_point(size = .2) +
          theme(panel.background = element_rect(fill = "white"))+
          scale_color_manual(values=cols)+
          ggtitle(paste0(colnames(cut2.df)[k])))
  dev.off()
  
  
}

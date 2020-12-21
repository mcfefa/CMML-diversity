library(dplyr)
library(harmony)
library(ggplot2)
library(patchwork)
library(BiocNeighbors)
library(KSgeneral)
library(Matching)
library(uwot)
library(igraph)

rm(list=ls())

#This is the directory where I have files stored for the diff. potential of each sample
filedir <- '/Users/4472241/scCode/palantirEmbeddings/mapToRep1/'

#These samples are either HMA or some other Tx 
#We include HMA as a separate category, but exclude the Rux and Chemo samples
hmaList <- c(17,6,5,7,32,33,35,37,39)
excludeList <- c(5, 34, 36, 38)

#Import the diff potential for each individual and add it to the metadata to construct a violin plot
for (i in 1:46){
  
  #Create sampleID (used to ID file while reading csv) and assign condition (HMA, normal, or cmml)
  if (i < 8){
    sampleID <- paste0('normal', i)
    condition <- "normal"
  } else{
    sampleID <- paste0('cmml', i-7)
    condition <- "cmml"
    if ((i-7) %in% hmaList){
      condition <- 'HMA'
    }
    if ((i-7) %in% excludeList){
      condition <- 'Other'
    }
  }
  
  #Read the csv containing diff. potential info for each cell in the sample
  diffPotential <- read.csv(paste0(filedir, sampleID, 'branchProbs_diffPot.csv'))
  diffPot <- diffPotential[,'DiffPot']
  diffPot <- data.frame("DiffPotential" = diffPot, "Treatment" = condition)
  
  #Combine all the cells of all the samples into one large dataframe with condition labeled
  if (i == 1){
    diffPotAll <- diffPot
  } else{
    diffPotAll <- rbind(diffPotAll, diffPot)
  }
  rm(diffPot)
}

#Exclude those that aren't hma, cmml, or normal (rux and chemo)
diffPotAll <- diffPotAll[diffPotAll$Treatment!='Other',]

#Save as pdf
pdf('/Users/4472241/scCode/Paper_Figures/mainPaper/diffPotentialVlnPlotsByCondition_MatchUMAPColors.pdf', width = 10, height = 10)

#Violin plot (keep colors used for sample UMAP to assign CMML vs HMA vs Normal
ggplot(diffPotAll, aes(x = DiffPotential, y=Treatment, fill=Treatment, group = Treatment)) + geom_violin() + 
  guides(fill=FALSE)+theme(panel.background = element_rect(fill = "white")) + 
  scale_fill_manual(values=c("#B5D8E4", "#D43130", "#BEBEBD"))
  
dev.off()

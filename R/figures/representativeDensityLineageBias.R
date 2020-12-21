library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
library(patchwork)
library(BiocNeighbors)

rm(list = ls())

#Make representative density plot for over/under branch expression
#One expressing Ery Bias and one expressing Mono Bias

inputfiledir <- "/Users/4472241/scCode/palantirEmbeddings/mapToRep1/"
outputfiledir <- "/Users/4472241/scCode/Paper_Figures/"


#First make Erythrocyte bias plot, for this we use sample 11
cmml11Embeddings <- read.csv(paste0(inputfiledir, 'cmml11PalantirEmbeddings.csv'))
write.csv(cmml11Embeddings, paste0(outputfiledir,'dataMain/cmml11PalantirEmbeddings.csv'))

pdf(paste0(outputfiledir,"/mainPaper/cmml11embeddingsBinHex.pdf"),width=8, height=5)
print(ggplot(cmml11Embeddings,aes(x=X1,y=X2)) + stat_binhex(aes(fill=log10(..density..)),binwidth = 1) + 
        scale_fill_gradient2(midpoint = -2.8,low="red",mid="white", high="blue",name = "Log Fraction in Hexagon",na.value=NA,
                             limits = c(-3.5,-1.5))+
        theme(panel.background = element_rect(fill = "white")))
dev.off()

cmml25Embeddings <- read.csv(paste0(inputfiledir, 'cmml25PalantirEmbeddings.csv'))
write.csv(cmml25Embeddings, paste0(outputfiledir,'dataMain/cmml25PalantirEmbeddings.csv'))

pdf(paste0(outputfiledir,"/mainPaper/cmml24embeddingsBinHex.pdf"),width=8, height=5)
print(ggplot(cmml25Embeddings,aes(x=X1,y=X2)) + stat_binhex(aes(fill=log10(..density..)),binwidth = 1) + 
        scale_fill_gradient2(midpoint = -2.8,low="red",mid="white", high="blue",name = "Log Fraction in Hexagon",na.value=NA,
                             limits = c(-3.5,-1.5))+
        theme(panel.background = element_rect(fill = "white")))
dev.off()

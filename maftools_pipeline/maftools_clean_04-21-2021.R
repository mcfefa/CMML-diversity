
library(maftools)

rm(list = ls())

#Define functions for color mapping the continuous variables
colorMap <- function(vec, min=-3.4, max=3.4){
  mid = (max+min)/2
  output <- c()
  for (i in 1: length(vec)){
    r <- (vec[i]-mid)/(min-mid)
    b <- (vec[i]-mid)/(max-mid)
    g <- 0
    if (b>0){
      output[i] <- rgb(0, g, 1, b)
    }else{
      output[i] <- rgb(1, g, 0, r)
    }
  }
  print(output)
}

colorMapBlast <- function(vec, min=0, max=20){
  output <- c()
  for (i in 1: length(vec)){
    r <- 0
    b <- (vec[i])/(max)
    g <- 0
    output[i] <- rgb(0, g, 1, b)
  }
  print(output)
}

colorMapToSample <- function(vec, min=0, max=36){
  output <- c()
  for (i in 1: length(vec)){
    r <- 0
    b <- (vec[i])/(max)
    g <- 0
    output[i] <- rgb(0, g, 1, b)
  }
  print(output)
}

#path to  MAF file
laml.maf = maf = '/Users/4472241/scCode/MAF_sampleOnly_IWGifBoth_no4J.txt'
#clinical information containing lots of information
laml.clin = '/Users/4472241/scCode/clinical_data_forMAF.txt'

laml = read.maf(maf = laml.maf, clinicalData = laml.clin, removeDuplicatedVariants = FALSE)

#Shows sample summry.
#getSampleSummary(laml)
#Shows gene summary.
#getGeneSummary(laml)
#shows clinical data associated with samples
#getClinicalData(laml)
#Shows all fields in MAF
#getFields(laml)
#Writes maf summary to an output file with basename laml.
#write.mafSummary(maf = laml, basename = 'laml')
#plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)


laml@clinical.data[,Age := as.numeric(as.character(Age))]
laml@clinical.data[,Blast_Percent := as.numeric(as.character(Blast_Percent))]
laml@clinical.data[,Palantir_HSC := as.numeric(as.character(Palantir_HSC))]
laml@clinical.data[,Eppert_HSC := as.numeric(as.character(Eppert_HSC))]
laml@clinical.data[,Time_From_Diagnosis_to_Sample := as.numeric(as.character(Time_From_Diagnosis_to_Sample))]

Blast_Percent_Color <- colorMapBlast(laml@clinical.data[,Blast_Percent], min = 0, max = 20)
names(Blast_Percent_Color) <- laml@clinical.data[,Blast_Percent]
Eppert_HSC_Color <- colorMap(laml@clinical.data[,Eppert_HSC])
names(Eppert_HSC_Color) <- laml@clinical.data[,Eppert_HSC]
ToSample_Color <- colorMapToSample(laml@clinical.data[,Time_From_Diagnosis_to_Sample])
names(ToSample_Color) <- laml@clinical.data[,Time_From_Diagnosis_to_Sample]

colsAnn = list(
  'Blast_Percent' = Blast_Percent_Color,
  'Eppert_HSC' = Eppert_HSC_Color,
  'Time_From_Diagnosis_to_Sample' = ToSample_Color,
  'MDACC' = c('1' = 'green', '2' = 'yellow', '3' = 'orange', '4' = 'red', 'Post-BMT' = 'black',
              'AML' = 'pink'),
  'Mayo' = c('1' = 'green', '2' = 'yellow', '3' = 'red', 'AML' = 'pink', 'Post-BMT' = 'black'),
  'FAB' = c('1' = 'deepskyblue', '2' = 'deeppink', 'AML' = 'pink', 'Post-BMT' = 'black'),
  'WHO' = c('0' = 'green', '1' = 'yellow', '2' = 'red', 'AML' = 'pink', 'Post-BMT' = 'black'),
  'Tx_Naive' = c('yes' = 'black', 'no' = 'red')
)

clin_data = data.table::fread(input = laml.clin)
clin_data_sorted = clin_data[order(Eppert_HSC)]


#Oncoplot with clinical data
pdf(file = "/Users/4472241/scCode/oncoPlots/oncoPlot_continuousWithTimeToSample_onlyTimeOfSample_orderByWHO.pdf", width = 12, height = 10, paper = "special", bg = "white")
oncoplot(maf = laml, clinicalFeatures = c(#"WHO",
  #"MDACC", "Mayo", "FAB", "Tx_Naive"),
  "Blast_Percent", "Eppert_HSC", "Time_From_Diagnosis_to_Sample"),
  sampleOrder = clin_data_sorted$Tumor_Sample_Barcode, top = 20,
  annotationColor = colsAnn,
  showTumorSampleBarcodes = TRUE, removeNonMutated = FALSE,
  anno_height = 15/4, barcode_mar = 8)
dev.off()
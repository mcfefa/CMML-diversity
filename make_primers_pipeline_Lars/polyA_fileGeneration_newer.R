#script to tell us which sample to start at and how many cb's until next sample

library(biomaRt)
rm(list=ls())

#Read in arguments from bash
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
iteration <- as.numeric(args[2])

#Read csv
df <- read.csv(input_file)

#Set base directory
baseDir <- '/share/lab_altrock/MeghanCluster/BrianCluster/make_primers/'

#Split into list
df.split <- split(df, df$Tumor_Sample_Barcode)
df <- df.split[[iteration]]

#Get the sample name in "CMML_#" format
sampleNumber <- unique(df$Sample.Number)

#Get gene name list
geneVec <- unique(df$Hugo_Symbol)

#Get position of start and end for given genes
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
bed <- getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'strand'),
      filters=c('hgnc_symbol'),
      values=geneVec,
      mart=ensembl)

#Write this out to a bed file
if (!dir.exists(paste0(baseDir, sampleNumber))){
    dir.create(paste0(baseDir, sampleNumber))
    bed_outfile <- paste0(baseDir, sampleNumber, '/', sampleNumber,'.bed')
    write.table(bed, bed_outfile, row.names = F, sep = "\t", col.names = F)

    #Make csv with the variant info
    to_make_csv <- df
    to_make_csv$CHROM <- gsub("chr", "", df$Chromosome)
    to_make_csv$POS <- df$Start_Position
    to_make_csv$symbol <- df$Hugo_Symbol
    csv_outfile <- paste0(baseDir, sampleNumber, '/', sampleNumber,'_variants.csv')
    write.csv(to_make_csv, csv_outfile)
}


#Get the path to the individual RDS file
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
raw.number <- as.numeric(gsub("CMML_", "", sampleNumber))
index <- raw.number
path_to_RDS <- paste0('/share/lab_altrock/MeghanCluster/BrianCluster/individual-CMML-rds/', 
               rdsFileList[index], '.rds')


#Output back to bash file
cat(c(as.character(sampleNumber), bed_outfile, csv_outfile, path_to_RDS)) 

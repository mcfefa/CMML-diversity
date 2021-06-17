library(mitoClone)
library(deepSNV)
library(GenomicRanges)
library(devtools)
library(pheatmap)
library("reticulate")
#Sys.setenv(RETICULATE_PYTHON = "/Users/4472241/anaconda3/envs/mitoclone/bin/python", force = T)
knitr::opts_chunk$set(echo = TRUE)
py_config()

# GOAL: Run Mitoclone individually, then manually filter mutation sites

#Directions to run on cluster:
#screen
#request job
#conda activate mitocloneV2
#open R
#Run this file

rm(list=ls())

#Open and merge all of the previously computed
#dirFirstTen <- '/share/lab_altrock/MeghanCluster/BrianCluster_PA/BAM/mitocloneOutput_firstTen/'
#dirSecondTen <- '/share/lab_altrock/MeghanCluster/BrianCluster_PA/BAM/mitocloneOutput_secondTen/'
#dirThirdTen <- '/share/lab_altrock/MeghanCluster/BrianCluster_PA/BAM/mitocloneOutput_thirdTen/'
#dirFourthTen <- '/share/lab_altrock/MeghanCluster/BrianCluster_PA/BAM/mitocloneOutput_fourthTen_defaultParams_MINREADS3_MINCELL10_MINFRAC0.05_MINRelative.Other0.99/'

#countTables_firstTen <- readRDS(paste0(dirFirstTen, 'countTablesTen.rds'))
#countTables_secondTen <- readRDS(paste0(dirSecondTen, 'countTablesSecondTen.rds'))
#countTables_thirdTen <- readRDS(paste0(dirThirdTen, 'countTablesThirdTen.rds'))
#countTables_fourthTen <- readRDS(paste0(dirFourthTen, 'countTablesFourthTen.rds'))

#Merge the countTables and save the larger RDS file in both directories (altrock and padron)
#countTables_allSamples <- c(countTables_firstTen, countTables_secondTen, countTables_thirdTen, countTables_fourthTen)
#countTables <- countTables_allSamples

#This is how we would generate count tables from individual bam files, but count tables makes this unnecessary
#dir <- '/share/lab_padron/Meghan/scRNAseq/BrianCluster/BAM/mitocloneOutput_firstTen/' #For cluster
#bamdir <- '/share/lab_altrock/MeghanCluster/BrianCluster/BAM/firstTenSamples_RustOutput_UseThisOne' #For cluster
#bamdir <- '/Users/4472241/scCode/mitoclone/cmml_bamFiles/rustOutput' #For local
#dir <- '/Users/4472241/scCode/mitoclone' #For local
#Get bam files, exclude bai files (which are in same dir), but make sure each bam has a bai
#bamfiles <- list.files(bamdir, pattern = "bam", full.names = T)
#baifiles <- list.files(bamdir, pattern = 'bai', full.names = T)
#bamfiles <- bamfiles[!(bamfiles %in% baifiles)]
#baifilesCut <- gsub(".bai\\.*","",baifiles)
#bamfiles <- bamfiles[bamfiles %in% baifilesCut]
#Generate count table (list of matrices with one matrix for each bam file, so one for each cell)
#countTables <- baseCountsFromBamList( bamfiles,  sites = "MT:1-16569", ncores = 12)
#Then save as rds, because this^ will take a while
#saveRDS(countTables, paste0(dir, 'countTablesTen.rds'))

#Assign names to each matrix (cell)
#bamdir <- '/share/lab_altrock/MeghanCluster/BrianCluster_PA/BAM/firstTenSamples_RustOutput_UseThisOne'
#bamfiles <- list.files(bamdir, pattern = "bam", full.names = T)
#baifiles <- list.files(bamdir, pattern = 'bai', full.names = T)
#bamfiles <- bamfiles[!(bamfiles %in% baifiles)]
#baifilesCut <- gsub(".bai\\.*","",baifiles)
#bamfiles <- bamfiles[bamfiles %in% baifilesCut]
#namesFrom_bamFiles_firstTen <- sub(".*UseThisOne/", "", bamfiles)

#bamdir <- '/share/lab_altrock/MeghanCluster/BrianCluster_PA/BAM/secondTenSamples_RustOutput'
#bamfiles <- list.files(bamdir, pattern = "bam", full.names = T)
#baifiles <- list.files(bamdir, pattern = 'bai', full.names = T)
#bamfiles <- bamfiles[!(bamfiles %in% baifiles)]
#baifilesCut <- gsub(".bai\\.*","",baifiles)
#bamfiles <- bamfiles[bamfiles %in% baifilesCut]
#namesFrom_bamFiles_secondTen <- sub(".*RustOutput/", "", bamfiles)

#bamdir <- '/share/lab_altrock/MeghanCluster/BrianCluster_PA/BAM/thirdTenSamples_RustOutput'
#bamfiles <- list.files(bamdir, pattern = "bam", full.names = T)
#baifiles <- list.files(bamdir, pattern = 'bai', full.names = T)
#bamfiles <- bamfiles[!(bamfiles %in% baifiles)]
#baifilesCut <- gsub(".bai\\.*","",baifiles)
#bamfiles <- bamfiles[bamfiles %in% baifilesCut]
#namesFrom_bamFiles_thirdTen <- sub(".*RustOutput/", "", bamfiles)

#bamdir <- '/share/lab_altrock/MeghanCluster/BrianCluster_PA/BAM/fourthTenSamples_RustOutput'
#bamfiles <- list.files(bamdir, pattern = "bam", full.names = T)
#baifiles <- list.files(bamdir, pattern = 'bai', full.names = T)
#bamfiles <- bamfiles[!(bamfiles %in% baifiles)]
#baifilesCut <- gsub(".bai\\.*","",baifiles)
#bamfiles <- bamfiles[bamfiles %in% baifilesCut]
#namesFrom_bamFiles_fourthTen <- sub(".*RustOutput/", "", bamfiles)

#namesFrom_bamFiles <- c(namesFrom_bamFiles_firstTen, namesFrom_bamFiles_secondTen, namesFrom_bamFiles_thirdTen, namesFrom_bamFiles_fourthTen)
#names(countTables) <- namesFrom_bamFiles

#saveRDS(countTables, '/share/lab_altrock/MeghanCluster/BrianCluster_PA/countTablesAllForMitoclone_06-17-2021.rds')
#In terminal (scp faster than saveRDS): scp /share/lab_altrock/MeghanCluster/BrianCluster_PA/countTablesAllForMitoclone_06-17-2021.rds /share/lab_padron/Meghan/scRNAseq/BrianCluster/
#saveRDS(countTables, '/share/lab_padron/Meghan/scRNAseq/BrianCluster/countTablesAllForMitoclone_06-17-2021.rds')

#Once we have merged all the count tables, we can just load it from wherever it is located
dir <- '/share/lab_altrock/MeghanCluster/BrianCluster_PA/'
countTables <- readRDS(paste0(dir, "countTablesAllForMitoclone_06-17-2021.rds"))

#Assign names to each patient and run mutationCallsFromCohort
namesFrom_bamFiles <- names(countTables)
patient <- sub("_CB.*", "", namesFrom_bamFiles)

countTables.patientList <- split(countTables, patient)

#Try running all patients at once with default params...see what happens (expecting memory problem)
#result <- mutationCallsFromCohort(countTables, patient) #Will not work: memory error

#Instead run each individually, creating a blacklist which will have to be manually edited (remove repeats in non-sequential samples)
for (i in 1:39){

    CMML.count <- countTables.patientList[[i]]
    sample <- gsub("_CB.*", "", names(countTables.patientList[[i]])[1])
    CMML <- mutationCallsFromBlacklist(CMML.count, min.af=0.1, min.num.samples=5, universal.var.cells = 0.8 * length(CMML.count), binarize = 0.1)

    #CMML.Meta <- data.frame(row.names = rownames(CMML_1@N), Clone = gsub("_.*","",gsub("Donor1_","",rownames(CMML_1@N))))

    #Keep track of who is what site (to filter later)
    sites <- data.frame("Site" = row.names(as.data.frame(CMML@cluster)), "Sample" = sample)
    if (i ==1){
        sitesAll <- sites
        CMML_All <- list(CMML)
        patientOrdered <- sample
    }else{
        sitesAll <- rbind(sitesAll, sites)
        CMML_All <- append(CMML_All, list(CMML))
        patientOrdered <- c(patientOrdered, sample)
    }

}

names(CMML_All) <- patientOrdered
#Write to csv for analyzing which sites to keep and which to get rid of
write.csv(sitesAll, paste0(dir, 'siteListForMakingBlacklist_06-17-2021.csv'))
saveRDS(CMML_All, paste0(dir, 'mitocloneOutputList_unfilteredSites_runIndividually_06-17-2021.rds'))

########## ACTION REQUIRED #######################################################################################################
#Manually filter all of the mutations by looking through the csv of sites we created above
##################################################################################################################################

#Re-run checking against a blacklist
blacklist <- read.csv(paste0(dir, "blacklistMitoclone_06-17-2021.csv"))

#Check for errors by running first example
CMML.count <- countTables.patientList[[i]]
sample <- gsub("_CB.*", "", names(countTables.patientList[[i]])[1])
CMML <- mutationCallsFromBlacklist(CMML.count, min.af=0.1, min.num.samples=5, universal.var.cells = 0.8 * length(CMML.count), 
    binarize = 0.1, blacklists.use = blacklist)

#Loop to do everything over again, except this time we will also cluster
for (i in 1:39){

    CMML.count <- countTables.patientList[[i]]
    sample <- gsub("_CB.*", "", names(countTables.patientList[[i]])[1])
    CMML <- mutationCallsFromBlacklist(CMML.count, min.af=0.1, min.num.samples=5, universal.var.cells = 0.8 * length(CMML.count), 
        binarize = 0.1, blacklists.use = blacklist)

    #CMML.Meta <- data.frame(row.names = rownames(CMML_1@N), Clone = gsub("_.*","",gsub("Donor1_","",rownames(CMML_1@N))))

    #Make a list with all of the outputs
    sites <- data.frame("Site" = row.names(as.data.frame(CMML@cluster)), "Sample" = sample)
    if (i ==1){
        CMML_All <- list(CMML)
    }else{
        CMML_All <- append(CMML_All, list(CMML))
    }

    CMML <- muta_cluster(CMML, cores = 12, tempfolder = paste0(getwd(),"/", sample, "_temp"), force_recalc = T)
    saveRDS(CMML, paste0(dir, '/mitocloneOutputAll_06-17-2021/', sample, '_results_postFilter_postCluster_mitoclone_06-17-2021.rds'))

}

#Check that we have the expected output
print(colnames(CMML@M))

#Save the filtered RDS file with all of the samples (before clustering, may not be worth saving)
#saveRDS(CMML_All, paste0(dir, 'mitocloneOutputList_filteredSites_runIndividually_06-17-2021.rds'))




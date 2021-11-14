library(mitoClone)
library(deepSNV)
library(GenomicRanges)
library(pheatmap)
library(ggplot2)
library("reticulate")
Sys.setenv(RETICULATE_PYTHON = "./anaconda3/envs/mitoclone/bin/python", force = T)
knitr::opts_chunk$set(echo = TRUE)
py_config()

#Once we have merged all the count tables, we can just load it from wherever it is located
dir <- './mitocloneDir/'
countTables <- readRDS(paste0(dir, "countTablesAllForMitoclone.rds"))

# Edit the mutationCallsFromBlacklist function to remove multiple cores, which
# caused errors on my local computer
mutationCallsFromBlacklist <- function(BaseCounts,lim.cov=20, min.af=0.2, 
                                       min.num.samples=0.01*length(BaseCounts), 
                                       min.af.universal =min.af, 
                                       universal.var.cells=0.95*length(BaseCounts), 
                                       blacklists.use = blacklists, max.var.na = 0.5, 
                                       max.cell.na = 0.95, ...) {
  varaf <- parallel::mclapply(BaseCounts,function(x){
    ## focus on A,G,C,T
    x <- x[,1:4]
    ## find cell that have less than 100 cov over agct at a given pos
    zeroes <- rowSums(x) < lim.cov
    ## af calc
    #x.af <- x/rowSums(x)
    x.af <- x / (x+apply(x,1,max))
    x.af <- reshape2::melt(x.af)
    colnames(x.af) <- c('pos','nt','af')
    ## remove reference af's
    x.af <- x.af[!(mito.dna[x.af$pos] == x.af$nt),]
    ## remove N site
    x.af <- x.af[!(mito.dna[x.af$pos] == 'N'),]
    x.af$name <- paste0(x.af$pos,' ',mito.dna[x.af$pos],'>',x.af$nt)
    ## find dominant NT
    x.af$af[x.af$pos %in% which(zeroes)] <- NA
    x <- x.af$af
    names(x) <- x.af$name
    return(x)
  }, mc.cores=1) ####################### remove parallelism here ###############
  varaf <- do.call(cbind, varaf)
  ## you could allow for only sites with coverage! currently you filter at a rate of 10% cells dropping out max
  ##varaf <- varaf[rowSums(is.na(varaf))/length(mc.out) < max.fraction.na,]
  varaf <- varaf[rowSums(varaf > min.af,na.rm=TRUE) >= min.num.samples,]
  
  is.names <- sapply(blacklists.use, function(x) typeof(x) == "character")
  #part 2 - filter based on the blacklist
  if(sum(is.names) > 0){
    removal.names.list <- unique(unlist(blacklists.use[is.names]))
    varaf <- varaf[!row.names(varaf) %in% removal.names.list,]
  }
  if(sum(!is.names) > 0){
    removal.ranges.list <- unique(unlist(GenomicRanges::GRangesList(blacklists.use[!is.names])))
    varaf <- varaf[-c(S4Vectors::queryHits(GenomicRanges::findOverlaps(mut2gr(row.names(varaf)),removal.ranges.list))),]
  }
  #if(drop.empty){
  varaf <- varaf[rowSums(varaf,na.rm=T) > 0,] #colSums(varaf,na.rm=T) > 0
  #}
  
  varaf <- varaf[!rowSums(varaf >= min.af.universal,na.rm=TRUE) >= universal.var.cells,]
  ## vars must have less than X % NA's
  varaf <- varaf[rowSums(is.na(varaf)) < max.var.na*NCOL(varaf),]
  ## cells must have less than X % NA's
  varaf <- varaf[,colSums(is.na(varaf)) < max.cell.na*NROW(varaf)]
  
  MN <- pullcounts.vars(BaseCounts, rownames(varaf), colnames(varaf))
  mutationCallsFromMatrix(t(MN$M), t(MN$N), ...)
}

#Assign names to each patient and run mutationCallsFromCohort/Blacklist
namesFrom_bamFiles <- names(countTables)
patient <- sub("_CB.*", "", namesFrom_bamFiles)

countTables.patientList <- split(countTables, patient)
alphaVec <- c("A", "B", "C", "D", "E", "F", "G", "H")
sampleList <- list(c(1, 32, 39), c(13, 35), c(25, 33), c(27, 37), c(28, 38), c(15, 34), c(12, 36), c(26, 3))
names(sampleList) <- alphaVec

#Run each patient's sequential samples together
for (i in 1:8){
  
  #Get the right samples for each patient
  currentPatient <- alphaVec[i]
  samples <- sampleList[[currentPatient]]
  
  #Combine samples to one count table for each patient (3 samples for A, 2 for others)
  CMML.count1 <- countTables.patientList[[paste0("CMML_", samples[1])]]
  CMML.count2 <- countTables.patientList[[paste0("CMML_", samples[2])]]
  if (currentPatient == "A"){
    CMML.count3 <- countTables.patientList[[paste0("CMML_", samples[3])]]
    combinedPatient.count <- c(CMML.count1, CMML.count2, CMML.count3)
  }else{
    combinedPatient.count <- c(CMML.count1, CMML.count2)
  }
  
  #CMML.count <- countTables.patientList[[paste0("CMML_", samples[1])]]
  samples <- unique(gsub("_CB.*", "", names(combinedPatient.count)))
  
  # Run mutationCallsFromBlacklist to filter sites and assign as mutated or not in each cell
  min.af <- 0.2
  min.num.samples.frac <- 0.01
  universal.var.cells.frac <- 0.95
  max.var.na <- 0.5
  max.cell.na <- 0.75
  CMML <- mutationCallsFromBlacklist(combinedPatient.count, min.af=min.af, min.num.samples=min.num.samples.frac * length(combinedPatient.count), 
                                     universal.var.cells = universal.var.cells.frac * length(combinedPatient.count), binarize = 0.1,
                                     max.var.na = max.var.na, max.cell.na = max.cell.na)
  
  # Keep track of which patients have mutations at what site (to filter later)
  # Make list of output for each patient
  sites <- data.frame("Site" = row.names(as.data.frame(CMML@cluster)), "Sample" = samples[1])
  if (i ==1){
    sitesAll <- sites
    CMML_All <- list(CMML)
    #patientOrdered <- sample
  }else{
    sitesAll <- rbind(sitesAll, sites)
    CMML_All <- append(CMML_All, list(CMML))
    #patientOrdered <- c(patientOrdered, sample)
  }
  
}

names(CMML_All) <- alphaVec
#Write to csv for analyzing which sites to keep and which to get rid of
write.csv(sitesAll, paste0(dir, 'siteListForMakingBlacklist_sequentials_probablyDontNeed_07-22-2021.csv'))
# Save output for use in next steps (after filtering sites)
saveRDS(CMML_All, paste0(dir, 'mitocloneOutputList_unfilteredSites_sequentials_07-22-2021.rds'))

########## ACTION REQUIRED #######################################################################################################
#Manually filter all of the mutations by looking through the csv of sites we created above
##################################################################################################################################

#Before clustering, set sites in blacklist to false
dir <- '/Users/Brian/scCode/mitoclone/'
blacklist_csv <- read.csv(paste0(dir, "blacklistMitoclone_06-29-2021.csv"))
blacklist <- blacklist_csv$Blacklist[1:31]

CMML_All <- readRDS(paste0(dir, 'mitocloneOutputList_unfilteredSites_sequentials_07-22-2021.rds'))

allSeurat <- readRDS('/Users/Brian/Downloads/seuratObj_05-11-2021_postStandardPipeline_withHarmony_allSamples_39+8.rds')

#Configure into UMAP space
umap_embeddings <- data.frame(allSeurat@reductions[['umap']]@cell.embeddings)
cluster_info <- allSeurat$clusterResolution_0.05
#cmml_cellNames <- read.csv(paste0(dir, 'cmml', i, 'CellNames.csv'))
barcodes <- row.names(umap_embeddings)
barcodes <- sub(".*_", "", barcodes)
umap_embeddings$CB <- barcodes
umap_embeddings$Clus <- cluster_info
colnames(umap_embeddings) <- c( "X", "Y", "CB", "Clus")

#Import the csv matching cmml# to deID
number2ID <- read.csv(paste0(dir, 'numbering_to_deID.csv'))
number2ID$DeID <- gsub("LTB5109", "LTB6169", number2ID$DeID)

#Get the seurat sample name only
sampleNameOnly <- gsub("_.*", "", row.names(umap_embeddings))
sampleNameOnly <- gsub(".*-1", "HuaPt4", sampleNameOnly)
#sampleNameOnly <- gsub("LTB6169", "LTB5109", sampleNameOnly)
unique(sampleNameOnly)
#convert this to cmml number
cmmlNumberOnly <- c()
for (i in 1: length(sampleNameOnly)){
  if(sampleNameOnly[i] %in% number2ID$DeID){
    index <- which(number2ID$DeID %in% sampleNameOnly[i])
    cmmlNumberOnly[i] <- number2ID$Brian.Numbering[index]
  }else{
    cmmlNumberOnly[i] <- sampleNameOnly[i]
  }
}

for (i in 1:length(unique(sampleNameOnly))){
  if (! is.na(unique(cmmlNumberOnly)[i])){
    row.names(umap_embeddings) <- gsub(unique(sampleNameOnly)[i], unique(cmmlNumberOnly)[i], 
                                       row.names(umap_embeddings))
  }
}

row.names(umap_embeddings) <- paste0("CMML_", row.names(umap_embeddings))
barcodes <- row.names(umap_embeddings)
for (i in 1:length(CMML_All)){
  
  object <- CMML_All[[i]]
  names(CMML_All)[i]
  sites <- object@cluster
  
  # Remove sites which are in blacklist
  sites.remove <- which(names(sites) %in% blacklist)
  object@cluster[sites.remove] <- F
  
  #Cluster and plot, skipping samples that don't have any sites (which will error)
  tryCatch({
    CMML <- muta_cluster(object, cores=4, tempfolder = paste0(getwd(),"/CMML_temp"), force_recalc = T)
    #CMML <- quick_cluster(object)
    if (sum(CMML@cluster) > 0){
      CMML <- clusterMetaclones(CMML, min.lik =1)
      dev.off()
      pdf(paste0(dir, names(CMML_All)[i], "plotClones.pdf"))
      print(plotClones(CMML))
      dev.off()
      
      saveRDS(CMML, paste0(dir, 'clustered_mito_output_07-05-2021/', names(CMML_All)[i], '.rds'))
      
      #Get clonal info from mutaCluster output
      cmml_clonalInfo <- data.frame(CMML@cell2clone)
      clone <- data.frame("CB"=row.names(cmml_clonalInfo), "Clone" = NA)
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
      for (j in 1:length(clone$CB)){
        if (clone$CB[j] %in% barcodes){
          combine.df[j,"Clone"] <- clone$Clone[j]
          combine.df[j,"CB"] <- clone$CB[j]
          index <- which(row.names(umap_embeddings) %in% clone$CB[j])
          combine.df[j,"X"] <- umap_embeddings[index, "X"]
          combine.df[j,"Y"] <- umap_embeddings[index,"Y"]
          combine.df[j, "Clus"] <- umap_embeddings[index, "Clus"]
          combine.df[j, "Sample"] <- clone$Sample[j]
        }
      }
      
      #Remove nas, convert clone to factor in combine.df and plot
      combine.df <- combine.df[!is.na(combine.df$CB),]
      combine.df$Clone <- as.factor(combine.df$Clone)
      
      pdf(paste0(dir, 'mitocloneInUMAP/UMAP_mitocloneColored_CMML_', names(CMML_All)[i], '.pdf'))
      print(ggplot()+ geom_point(data = umap_embeddings, aes(x = X, y = Y), size = 0.05, color = "grey")+
              geom_point(data = combine.df, aes(x=X, y=Y, color = Clone), size = .2, alpha = 0.8) +
              theme(panel.background = element_rect(fill = "white"))+
              scale_color_manual(values=c("blue", "green", "red", "pink", 'orange', 'yellow',
                                          'purple', 'gray', 'cyan')))
      dev.off()
      
      #Make table for combined output and save it
      two_way_table <- with(combine.df, table(Clone, Clus))
      colnames(two_way_table) <- as.character(as.numeric(colnames(two_way_table))-1)
      write.csv(two_way_table, paste0(dir, 'cluster_clone/two_way_table_', names(CMML_All)[i], '.csv'))
      
      #Make table for each sample too
      split.df <- split(combine.df, combine.df$Sample)
      for (i in 1:length(split.df)){
        #Make table for combined output and save it
        two_way_table <- with(split.df[[i]], table(Clone, Clus))
        clone_frac <- rowSums(two_way_table)
        colnames(two_way_table) <- as.character(as.numeric(colnames(two_way_table))-1)
        write.csv(two_way_table, paste0(dir, 'cluster_clone_sequentialCombined/two_way_table_', names(split.df)[i], '.csv'))
        write.csv(clone_frac, paste0(dir, 'cluster_clone_sequentialCombined/clone_frac_', names(split.df)[i], '.csv'))
      }
    }
  }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
  rm(CMML, object)
}

#Check that we have the expected output
print(colnames(CMML@M))

CMML <- clusterMetaclones(CMML, min.lik =1)
plotClones(CMML)

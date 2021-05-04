library(mitoClone)
library(deepSNV)
library(GenomicRanges)
library(devtools)
library(pheatmap)
library("reticulate")
py_config()

rm(list=ls())

# GOAL: Run Mitoclone

#Directions to run on cluster:
#screen
#request job
#conda activate mitocloneV2
#open R
#Run this file

#Define custom "pullcounts.vars" function "my.pullcounts.vars"
my.pullcounts.vars <- function(mc.out,vars, cells=NULL, shift=0){
  pos <- as.integer(gsub(" *[NACGT].+","",vars)) + shift
  ref <- gsub("\\d+ *([NACGT])>(.+)","\\1",vars)
  alt <- gsub("\\d+ *([NACGT])>(.+)","\\2",vars)
  N <- sapply(mc.out, function(cell) {
    mapply(function(p,x) cell[p,x], pos, ref)
  })
  M <- sapply(mc.out, function(cell) {
    mapply(function(p,x) cell[p,x], pos, alt)
  })
  if(!is.matrix(M)) {
    M <- matrix(M, ncol = length(M),dimnames = list(vars, names(M)))
    N <- matrix(N, ncol = length(N),dimnames = list(vars,names(N)))
  }
  rownames(M) <- vars -> rownames(N)
  if (is.null(cells)) return(list(M = M, N = N)) else return(list(M = M[,cells], N = N[,cells]))
}

#Set dir
dir <- '/share/lab_altrock/MeghanCluster/BrianCluster/BAM/mitocloneOutput_thirdTen_defaultParams_MINREADS3_MINCELL10_MINFRAC0.05_MINRelative.Other0.99/' #For cluster
setwd(dir)
bamdir <- '/share/lab_altrock/MeghanCluster/BrianCluster/BAM/thirdTenSamples_RustOutput' #For cluster
#bamdir <- '/Users/4472241/scCode/mitoclone/cmml_bamFiles/rustOutput' #For local
#dir <- '/Users/4472241/scCode/mitoclone' #For local

#Get bam files, exclude bai files (which are in same dir)
bamfiles <- list.files(bamdir, pattern = "bam", full.names = T)
baifiles <- list.files(bamdir, pattern = 'bai', full.names = T)
bamfiles <- bamfiles[!(bamfiles %in% baifiles)]
baifilesCut <- gsub(".bai\\.*","",baifiles)
bamfiles <- bamfiles[bamfiles %in% baifilesCut]

#Generate count table (list of matrices with one matrix for each bam file, so one for each cell)
#countTables <- baseCountsFromBamList( bamfiles,  sites = "MT:1-16569", ncores = 12)

#Then save as rds, because this^ will take a while
#saveRDS(countTables, paste0(dir, '/countTablesThirdTen.rds'))
countTables <- readRDS('/share/lab_altrock/MeghanCluster/BrianCluster/BAM/mitocloneOutput_thirdTen/countTablesThirdTen.rds')

#Assign names to each matrix (cell)
namesFrom_bamFiles <- sub(".*RustOutput/", "", bamfiles)
names(countTables) <- namesFrom_bamFiles

#Assign names to each patient and run mutationCallsFromCohort
patient <- sub("_CB.*", "", namesFrom_bamFiles)

#Paste function from github (wasn't working in called function: mutationCallsFromCohort)
#result <- mutationCallsFromCohort(countTables, patient)

################### INSIDE FUNCTION ##############################
BaseCounts <- countTables
#patient
MINREADS <- 3
MINCELL <- 10
MINFRAC <- 0.05
MINCELLS.PATIENT <- 10
MINRELATIVE.PATIENT <- 0.01
MINRELATIVE.OTHER <-  0.99
nuc.count.per.position.array <- array(data = 0,
                                    dim = c(length(BaseCounts),
                                            nrow(BaseCounts[[1]]),
                                            ncol(BaseCounts[[1]])),
                                    dimnames = list(names(BaseCounts),
                                                    paste0("X", 1:nrow(BaseCounts[[1]])),
                                                    colnames(BaseCounts[[1]])
                                    ))
  

for (i in 1:length(BaseCounts)) nuc.count.per.position.array[i,,] <- BaseCounts[[i]]

#determine the overall reference
sum.overall <- apply(nuc.count.per.position.array, c(2,3), sum)
reference <- colnames(sum.overall)[apply(sum.overall, 1, which.max)]

mt.reads.per.cell <- apply(nuc.count.per.position.array, 1, sum)

#turn the array into a binary array of VARIANTS
variant_calls <- lapply(1:length(reference), function(pos) {
    #which variants exist at this site?
    support <- apply(nuc.count.per.position.array[,pos,] > MINREADS,2,sum )
    support <- support[!names(support) %in% c(reference[pos], "N")]
    use <- names(support)[support > MINCELL]

    if (length(use) == 0) NULL else {
      out <- matrix(data =NA, ncol = length(use), nrow = nrow(nuc.count.per.position.array[,pos,]),
                    dimnames = list(rownames(nuc.count.per.position.array[,pos,]), paste0(pos,reference[pos],">",use)))
      for (i in 1:length(use)) {
        pos_sum <- apply(nuc.count.per.position.array[,pos,],1,sum)
        condition_mut <- nuc.count.per.position.array[,pos,use[i]] > MINREADS &  nuc.count.per.position.array[,pos,use[i]] > MINFRAC * pos_sum
        condition_ref <- nuc.count.per.position.array[,pos,reference[pos]] > MINREADS &  nuc.count.per.position.array[,pos,reference[pos]] > MINFRAC * pos_sum

        out[, i] <- ifelse(condition_mut,
                           ifelse(condition_ref, "BOTH", "MUT"),
                           ifelse(condition_ref, "WT", "DROP")
        )
      }
      out
    }

})

variant_calls <- do.call(cbind, variant_calls)

#now check:
#how often does each variant exist per patient?
varcount.bypatient <- sapply(unique(patient), function(pa) {
    apply(variant_calls[patient == pa, ] ,2, function(x) sum(x %in% c("BOTH","MUT")))
})


patient.count <-  as.vector(table(patient)[colnames(varcount.bypatient)])
names(patient.count) <- colnames(varcount.bypatient)
varcount.bypatient.fraction <- t(t(varcount.bypatient) / patient.count)

#throw out anything with less than
#a) 10 cells and 1% support in any patient
filter <- apply(varcount.bypatient, 1, max) > MINCELLS.PATIENT & apply(varcount.bypatient, 1, function(x) max(x) / patient.count[which.max(x)] ) > MINRELATIVE.PATIENT

#b) support of more than 1000 cells or 99% the level in a second patient (i.e. exclude none...for now)
patientfilter <- filter & apply(varcount.bypatient.fraction, 1, function(x) sum(x > MINRELATIVE.OTHER*max(x)) ) == 1 &
    apply(varcount.bypatient, 1, function(x) sum(x >= 1000) ) == 1

  #what are the variants that occur in multiple patients abundantly?
multipatient <- filter & apply(varcount.bypatient.fraction, 1, function(x) sum(x > MINRELATIVE.OTHER*max(x)) ) > 1 &
    apply(varcount.bypatient, 1, function(x) sum(x >= MINCELLS.PATIENT) ) > 1 & !grepl(">-",rownames(varcount.bypatient))

mutation.bypatient <- colnames(varcount.bypatient)[apply(varcount.bypatient[patientfilter,],1,which.max)]

variant_calls_selected <- variant_calls[,patientfilter]

  #now, prepare return values.
mutation.bypatient <- mutation.bypatient[!grepl("->",colnames(variant_calls_selected))]
variant_calls_selected <- variant_calls_selected[,!grepl("->",colnames(variant_calls_selected))]

out <- lapply(unique(patient), function(pa) {
    #a, retrieve matrices of allele counts for patient specific variants
    if (sum(mutation.bypatient == pa) == 0) return(NULL)
    MN <- my.pullcounts.vars(BaseCounts[patient == pa], colnames(variant_calls_selected)[mutation.bypatient == pa])
    #b, create mutationCalls object
    o <- mutationCallsFromMatrix(t(MN$M), t(MN$N))
})

names(out) <- unique(patient)

#discuss this part with ben - not consistent with his blacklist.
out$blacklist <- rownames(varcount.bypatient[multipatient,])
out$blacklist <- gsub("(\\d+)(.+)","\\1 \\2", out$blacklist)
result <- out

##################### OUTSIDE Function ###################################

saveRDS(result, paste0(dir, 'result_04-26-2021.rds'))
#result <- readRDS(paste0(dir, 'result_04-26-2021.rds'))

#Check first 3
print(colnames(result$CMML_21@M))
print(colnames(result$CMML_22@M))
print(colnames(result$CMML_23@M))
print(colnames(result$CMML_24@M))
print(colnames(result$CMML_25@M))
print(colnames(result$CMML_26@M))
print(colnames(result$CMML_27@M))
print(colnames(result$CMML_28@M))
print(colnames(result$CMML_29@M))
print(colnames(result$CMML_30@M))

#Assign patient 21 through 30
P21 <- result$CMML_21
P22 <- result$CMML_22
P23 <- result$CMML_23
P24 <- result$CMML_24
P24@cluster[which((gsub(".*\\.\\.", "", colnames(result$CMML_24@M)) %in% ""))] <- F
P25 <- result$CMML_25
P26 <- result$CMML_26
P27 <- result$CMML_27
P28 <- result$CMML_28
P29 <- result$CMML_29
P30 <- result$CMML_30
P30@cluster[which((gsub(".*\\.\\.", "", colnames(result$CMML_30@M)) %in% ""))] <- F

#Construct most likely phylogenetic tree and save
P21 <- muta_cluster(P21, cores = 12, tempfolder = paste0(getwd(),"/CMML_21_temp"), force_recalc = T)
saveRDS(P21, paste0(dir, 'P21_results_mitoclone.rds'))
P22 <- muta_cluster(P22, cores = 12, tempfolder = paste0(getwd(),"/CMML_22_temp"), force_recalc = T)
saveRDS(P22, paste0(dir, 'P22_results_mitoclone.rds'))
P23 <- muta_cluster(P23, cores = 12, tempfolder = paste0(getwd(),"/CMML_23_temp"), force_recalc = T)
saveRDS(P23, paste0(dir, 'P23_results_mitoclone.rds'))
P24 <- muta_cluster(P24, cores = 12, tempfolder = paste0(getwd(),"/CMML_24_temp"), force_recalc = T)
saveRDS(P24, paste0(dir, 'P24_results_mitoclone.rds'))
P25 <- muta_cluster(P25, cores = 12, tempfolder = paste0(getwd(),"/CMML_25_temp"), force_recalc = T)
saveRDS(P25, paste0(dir, 'P25_results_mitoclone.rds'))
P26 <- muta_cluster(P26, cores = 12, tempfolder = paste0(getwd(),"/CMML_26_temp"), force_recalc = T)
saveRDS(P26, paste0(dir, 'P26_results_mitoclone.rds'))
P27 <- muta_cluster(P27, cores = 12, tempfolder = paste0(getwd(),"/CMML_27_temp"), force_recalc = T)
saveRDS(P27, paste0(dir, 'P27_results_mitoclone.rds'))
P28 <- muta_cluster(P28, cores = 12, tempfolder = paste0(getwd(),"/CMML_28_temp"), force_recalc = T)
saveRDS(P28, paste0(dir, 'P28_results_mitoclone.rds'))
P29 <- muta_cluster(P29, cores = 12, tempfolder = paste0(getwd(),"/CMML_29_temp"), force_recalc = T)
saveRDS(P29, paste0(dir, 'P29_results_mitoclone.rds'))
P30 <- muta_cluster(P30, cores = 12, tempfolder = paste0(getwd(),"/CMML_30_temp"), force_recalc = T)
saveRDS(P30, paste0(dir, 'P30_results_mitoclone.rds'))
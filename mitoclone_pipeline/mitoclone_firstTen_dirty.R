library(mitoClone)
library(deepSNV)
library(GenomicRanges)
library(devtools)
library(pheatmap)
library("reticulate")
Sys.setenv(RETICULATE_PYTHON = "/Users/4472241/anaconda3/envs/mitoclone/bin/python", force = T)
knitr::opts_chunk$set(echo = TRUE)
py_config()

rm(list=ls())

# GOAL: Run Mitoclone

#Directions to run on cluster:
#screen
#request job
#conda activate mitocloneV2
#open R
#Run this file

#Set dir
dir <- '/share/lab_altrock/MeghanCluster/BrianCluster/BAM/mitocloneOutput_firstTen/' #For cluster
bamdir <- '/share/lab_altrock/MeghanCluster/BrianCluster/BAM/firstTenSamples_RustOutput_UseThisOne' #For cluster
#bamdir <- '/Users/4472241/scCode/mitoclone/cmml_bamFiles/rustOutput' #For local
#dir <- '/Users/4472241/scCode/mitoclone' #For local

#Get bam files, exclude bai files (which are in same dir)
bamfiles <- list.files(bamdir, pattern = "bam", full.names = T)
baifiles <- list.files(bamdir, pattern = 'bai', full.names = T)
bamfiles <- bamfiles[!(bamfiles %in% baifiles)]
baifilesCut <- gsub(".bai\\.*","",baifiles)
bamfiles <- bamfiles[bamfiles %in% baifilesCut]

#Generate count table (list of matrices with one matrix for each bam file, so one for each cell)
countTables <- baseCountsFromBamList( bamfiles,  sites = "MT:1-16569", ncores = 12)
#Then save as rds, because this^ will take a while
saveRDS(countTables, paste0(dir, 'countTablesTen.rds'))

#countTables <- readRDS(paste0(dir, 'countTablesTen.rds'))

#Assign names to each matrix (cell)
namesFrom_bamFiles <- sub(".*UseThisOne/", "", bamfiles)
names(countTables) <- namesFrom_bamFiles

#Assign names to each patient and run mutationCallsFromCohort
patient <- sub("_CB.*", "", namesFrom_bamFiles)

#Paste function from github (wasn't working in called function: mutationCallsFromCohort)
#result <- mutationCallsFromCohort(countTables, patient)

################### INSIDE FUNCTION ##############################
BaseCounts <- countTables
#patient
MINREADS <- 5
MINCELL <- 20
MINFRAC <- 0.1
MINCELLS.PATIENT <- 10
MINRELATIVE.PATIENT <- 0.01
MINRELATIVE.OTHER <-  0.1
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

#b) support of more than 10 cells or 10% the level in a second patient
patientfilter <- filter & apply(varcount.bypatient.fraction, 1, function(x) sum(x > MINRELATIVE.OTHER*max(x)) ) == 1 &
    apply(varcount.bypatient, 1, function(x) sum(x >= MINCELLS.PATIENT) ) == 1

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
    MN <- pullcounts.vars(BaseCounts[patient == pa], colnames(variant_calls_selected)[mutation.bypatient == pa])
    #b, create mutationCalls object
    o <- mutationCallsFromMatrix(t(MN$M), t(MN$N))
})

#Check why the above lapply isn't working
for(count in 9:10){
    if (sum(mutation.bypatient == unique(patient)[count]) == 0){
        print("NULL")
        return(NULL)
    }
    print(count)
    MN <- pullcounts.vars(BaseCounts[patient == unique(patient)[count]], colnames(variant_calls_selected)[mutation.bypatient == unique(patient)[count]])
    #o <- mutationCallsFromMatrix(t(MN$M), t(MN$N))
}

#################### INSIDE 'pullcounts.vars' FUNCTION #######################
count <- 9
mc.out <- BaseCounts[patient == unique(patient)[count]]
vars <- colnames(variant_calls_selected)[mutation.bypatient == unique(patient)[count]]
cells <- NULL
shift <- 0
pos <- as.integer(gsub(" *[ACGT].+","",vars)) + shift
ref <- gsub("\\d+ *([ACGT])>(.+)","\\1",vars)
alt <- gsub("\\d+ *([ACGT])>(.+)","\\2",vars)
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

names(out) <- unique(patient)

#discuss this part with ben - not consistent with his blacklist.
out$blacklist <- rownames(varcount.bypatient[multipatient,])
out$blacklist <- gsub("(\\d+)(.+)","\\1 \\2", out$blacklist)
result <- out

##################### OUTSIDE Function ###################################

saveRDS(result, paste0(dir, 'result_04-21-2021.rds'))


print(colnames(result$CMML_1@M))
print(colnames(result$CMML_2@M))
print(colnames(result$CMML_3@M))
print(colnames(result$CMML_9@M))

#Assign patient 1 through 10
P1 <- result$CMML_1
P2 <- result$CMML_2
P3 <- result$CMML_3
P4 <- result$CMML_4
P5 <- result$CMML_5
P6 <- result$CMML_6
P7 <- result$CMML_7
P8 <- result$CMML_8
P9 <- result$CMML_9
P10 <- result$CMML_10


#Construct most likely phylogenetic tree
#P1 <- muta_cluster(P1, cores=4, tempfolder = paste0(getwd(),"/CMML_1_debug"), 
#                   python_env = 'mitoclone')
P1 <- muta_cluster(P1, cores = 12, tempfolder = paste0(getwd(),"/CMML_1_temp"), force_recalc = T)
P2 <- muta_cluster(P2, cores = 12, tempfolder = paste0(getwd(),"/CMML_2_temp"), force_recalc = T)
P3 <- muta_cluster(P3, cores = 12, tempfolder = paste0(getwd(),"/CMML_3_temp"), force_recalc = T)
P4 <- muta_cluster(P4, cores = 12, tempfolder = paste0(getwd(),"/CMML_4_temp"), force_recalc = T)
P5 <- muta_cluster(P5, cores = 12, tempfolder = paste0(getwd(),"/CMML_5_temp"), force_recalc = T)
P6 <- muta_cluster(P6, cores = 12, tempfolder = paste0(getwd(),"/CMML_6_temp"), force_recalc = T)
P7 <- muta_cluster(P7, cores = 12, tempfolder = paste0(getwd(),"/CMML_7_temp"), force_recalc = T)
P8 <- muta_cluster(P8, cores = 12, tempfolder = paste0(getwd(),"/CMML_8_temp"), force_recalc = T)
P9 <- muta_cluster(P9, cores = 12, tempfolder = paste0(getwd(),"/CMML_9_temp"), force_recalc = T)
P10 <- muta_cluster(P10, cores = 12, tempfolder = paste0(getwd(),"/CMML_10_temp"), force_recalc = T)

usedata <- P1@ternary[,P1@cluster]

rownames(usedata)
plotTree(P1, file = "P1.ps")
plotTree(P1, file = "/Users/4472241/scCode/CMML_1_tree.ps")
m2c <- getMut2Clone(P1)

CMML_1_tree <- clusterMetaclones(P1, min.lik = 1)
CMML_2_tree <- clusterMetaclones(P2, min.lik = 1)
CMML_3_tree <- clusterMetaclones(P3, min.lik = 1)
CMML_9_tree <- clusterMetaclones(P9, min.lik = 1)
plotClones(CMML_1_tree)
plotClones(CMML_2_tree)
plotClones(CMML_3_tree)
plotClones(CMML_9_tree)

saveRDS(P1, './P1_results_mitoclone.rds')

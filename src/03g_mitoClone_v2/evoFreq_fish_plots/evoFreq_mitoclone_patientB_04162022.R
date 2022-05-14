# Evo freq for mitoclone
library(EvoFreq)
library(gridExtra)
library(colormap)
library(viridis)

rm(list = ls())

# NOTE: This has the correct phylogeny for the clones present according to mitoclone for patient A

# Set dir to save to and from which we pull countTables
dir <- '/Users/Brian/scCode/mitoclone/'
saveDir <- paste0(dir, 'evoFreq_04042022_newParams/')
importDir <- paste0(dir, 'sequential_clus_clone_03272022_newParams/')

#######   CHANGE PATIENT #############
currentSample <- "patientB"
patient <- read.csv(paste0(importDir, currentSample, '_cluster_and_clone_info.csv'))


# Clone
patient.table <- table(patient$Sample, patient$Clone)
patient.mat <- matrix(patient.table, ncol = ncol(patient.table), dimnames = dimnames(patient.table))
patient.wide <- as.data.frame(t(100*patient.mat/rowSums(patient.mat)))

#######  CHANGE: SET INITIAL COMPOSITION (SET INITIAL TO 0 IF IT STEMS FROM ANOTHER CLONE) ##########
patient.wide <- cbind(c(.01, rep(0, dim(patient.wide)[1]-1)), patient.wide)

colnames(patient.wide) <- c("", "Sample 1", "Sample 2")

###### NUMBER CLONES FROM 1 TO N. NOTE: THIS WILL AFFECT THE PARENTS YOU SET
patient.wide$clones <- c(1:dim(patient.wide)[1])
########  CHANGE PARENTS TO FIT MITOCLONE PHYLOGENY (MAKE SURE CLONES CONSISTENT WITH ABOVE NUMBERING) #########
patient.wide$parents <- c(0, rep(1, dim(patient.wide)[1] - 1))



# Then get the frequency data. (Use ?get_evofreq for options)
freq_frame <- get_evofreq(patient.wide[,seq(1,dim(patient.wide)[2]-2)], patient.wide$clones, 
                          patient.wide$parents, clone_cmap = "YIGnBu")

source("~/scCode/evoFreq_fixColor.R")
freq_frame <- evoFreq_fixColor_clone(freq_frame)

# Create the plot (shown on the left below)
evo_freq_p <- plot_evofreq(freq_frame) + ylab("Relative Pop. Size") + ggtitle("Clonal Dynamics") + xlab(NULL) + 
  theme_bw() + theme(axis.text.x = element_text(face = "bold", family = "Arial", size = 12), 
                     axis.title.y = element_text(face = "bold", family = "Arial", size = 12),
                     plot.title = element_text(face = "bold", family = "Arial", size = 14))
print(evo_freq_p)




# Do for cluster now
patient.table <- table(patient$Sample, patient$Clus)
patient.mat <- matrix(patient.table, ncol = ncol(patient.table), dimnames = dimnames(patient.table))
patient.wide <- as.data.frame(t(100*patient.mat/rowSums(patient.mat)))

####### CHANGE: THIS WILL DEPEND ON CLUSTREE (IF SPAWNING FROM OTHER CLUSTER, SET TO 0) ########
patient.wide <- cbind(rep(.01, dim(patient.wide)[1]), patient.wide)
patient.wide[5, 2] <- 0.001

colnames(patient.wide) <- c("", "Sample 1", "Sample 2")


####### CHANGE: THIS WILL DEPEND ON THE CLUSTREE #########
patient.wide$parents <- rep(0, dim(patient.wide)[1])
patient.wide$clones <- colnames(patient.mat)


# Then get the frequency data. (Use ?get_evofreq for options)
freq_frame <- get_evofreq(patient.wide[,seq(1,dim(patient.wide)[2]-2)], patient.wide$clones, 
                          patient.wide$parents, clone_cmap = "jet")

freq_frame <- evoFreq_fixColor_clus(freq_frame)

# Create the plot (shown on the left below)
clus_freq_p <- plot_evofreq(freq_frame) + ylab("Relative Pop. Size") + ggtitle("Cluster Dynamics") + xlab(NULL) +
  theme_bw() + theme(axis.text.x = element_text(face = "bold", family = "Arial", size = 12), 
                     axis.title.y = element_text(face = "bold", family = "Arial", size = 12),
                     plot.title = element_text(face = "bold", family = "Arial", size = 14))
print(clus_freq_p)

#pdf(paste0(saveDir, "patient_D_clone_and_cluster.pdf"))

####### TAKE SCREENSHOT FOR NOW...PDF WASN'T WORKING FOR ME
grid.arrange(evo_freq_p, clus_freq_p, ncol = 1)

table(evo_freq_p[["data"]][["clone_id"]], evo_freq_p[["data"]][["plot_color"]])

table(clus_freq_p[["data"]][["clone_id"]], clus_freq_p[["data"]][["plot_color"]])

#dev.off()

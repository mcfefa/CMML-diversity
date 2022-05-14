# Evo freq for mitoclone
library(EvoFreq)
library(gridExtra)
library(colormap)
library(viridis)

rm(list = ls())

# NOTE: This has the correct phylogeny for the clones present according to mitoclone for patient H

# Set dir to save to and from which we pull countTables
dir <- '/Users/Brian/scCode/mitoclone/'
saveDir <- paste0(dir, 'evoFreq_04042022_newParams/')
importDir <- paste0(dir, 'sequential_clus_clone_03272022_newParams/')

#######   SET PATIENT #############
currentSample <- "patientH"
H <- read.csv(paste0(importDir, currentSample, '_cluster_and_clone_info.csv'))


# Clone
H.table <- table(H$Sample, H$Clone)
H.mat <- matrix(H.table, ncol = ncol(H.table), dimnames = dimnames(H.table))
H.wide <- as.data.frame(t(100*H.mat/rowSums(H.mat)))

#######  SET INITIAL COMPOSITION (set initial to 0 if a clone stems from another clone) ##########
H.wide <- cbind(c(.01,rep(0, dim(H.wide)[1]-1)), H.wide)

colnames(H.wide) <- c("", "Sample 1", "Sample 2")

###### NUMBER CLONES FROM 1 TO N. NOTE: THIS WILL AFFECT THE PARENTS YOU SET
H.wide$clones <- c(1:dim(H.wide)[1])
########  CHANGE PARENTS TO FIT MITOCLONE PHYLOGENY (MAKE SURE CLONES CONSISTENT WITH ABOVE NUMBERING) #########
H.wide$parents <- c(0, rep(1, dim(H.wide)[1] - 1))

row.names(H.wide) <- c(1, 2)
H.wide[2,2] <- 1.01

# Then get the frequency data. (Use ?get_evofreq for options)
freq_frame <- get_evofreq(H.wide[,seq(1,dim(H.wide)[2]-2)], H.wide$clones, 
                          H.wide$parents, clone_cmap = "YIGnBu")

# Manually set to single clone (was giving error previously)
freq_frame$draw_order <- rev(freq_frame$draw_order)
freq_frame$clone_id <- 1

# Run function to set color scheme
source("~/scCode/evoFreq_fixColor.R")
freq_frame <- evoFreq_fixColor_clone(freq_frame)

# Create the plot
evo_freq_p <- plot_evofreq(freq_frame) + ylab("Relative Pop. Size") + ggtitle("Clonal Dynamics") + xlab(NULL) + 
  theme_bw() + theme(axis.text.x = element_text(face = "bold", family = "Arial", size = 12), 
                     axis.title.y = element_text(face = "bold", family = "Arial", size = 12),
                     plot.title = element_text(face = "bold", family = "Arial", size = 14))
print(evo_freq_p)




# Do for cluster now
H.table <- table(H$Sample, H$Clus)
H.mat <- matrix(H.table, ncol = ncol(H.table), dimnames = dimnames(H.table))
H.wide <- as.data.frame(t(100*H.mat/rowSums(H.mat)))

####### Set initial frequency ########
H.wide <- cbind(rep(.01, dim(H.wide)[1]), H.wide)

colnames(H.wide) <- c("", "Sample 1", "Sample 2")

####### Set parents (here all arise separately) #########
H.wide$parents <- rep(0, dim(H.wide)[1])
H.wide$clones <- colnames(H.mat)


# Then get the frequency data. (Use ?get_evofreq for options)
freq_frame <- get_evofreq(H.wide[,seq(1,dim(H.wide)[2]-2)], H.wide$clones, 
                          H.wide$parents, clone_cmap = "jet")

# Set color scheme
freq_frame <- evoFreq_fixColor_clus(freq_frame)

# Create the plot
clus_freq_p <- plot_evofreq(freq_frame) + ylab("Relative Pop. Size") + ggtitle("Cluster Dynamics") + xlab(NULL) +
  theme_bw() + theme(axis.text.x = element_text(face = "bold", family = "Arial", size = 12), 
                     axis.title.y = element_text(face = "bold", family = "Arial", size = 12),
                     plot.title = element_text(face = "bold", family = "Arial", size = 14))
print(clus_freq_p)

#pdf(paste0(saveDir, "patient_D_clone_and_cluster.pdf"))

####### TAKE SCREENSHOT FOR NOW...PDF WASN'T WORKING FOR ME
grid.arrange(evo_freq_p, clus_freq_p, ncol = 1)

#dev.off()

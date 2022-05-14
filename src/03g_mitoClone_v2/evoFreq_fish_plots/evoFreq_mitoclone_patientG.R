# Evo freq for mitoclone
library(EvoFreq)
library(gridExtra)
library(colormap)
library(viridis)
library(RColorBrewer)

rm(list = ls())

# NOTE: This has the correct phylogeny for the clones present according to mitoclone for patient G

# Set dir to save to and from which we pull countTables
dir <- '/Users/Brian/scCode/mitoclone/'
saveDir <- paste0(dir, 'evoFreq_04042022_newParams/')
importDir <- paste0(dir, 'sequential_clus_clone_03272022_newParams/')


#######   SET PATIENT #############
currentSample <- "patientG"
G <- read.csv(paste0(importDir, currentSample, '_cluster_and_clone_info.csv'))


# Clone
G.table <- table(G$Sample, G$Clone)
G.mat <- matrix(G.table, ncol = ncol(G.table), dimnames = dimnames(G.table))
G.wide <- as.data.frame(t(100*G.mat/rowSums(G.mat)))

#######  CHANGE: SET INITIAL COMPOSITION (SET INITIAL TO 0 IF IT STEMS FROM ANOTHER CLONE) ##########
G.wide <- cbind(c(.01,rep(0, dim(G.wide)[1]-1)), G.wide)

colnames(G.wide) <- c("", "Sample 1", "Sample 2")

###### NUMBER CLONES FROM 1 TO N. NOTE: THIS WILL AFFECT THE PARENTS YOU SET
G.wide$clones <- c(1:dim(G.wide)[1])
########  CHANGE PARENTS TO FIT MITOCLONE PHYLOGENY (MAKE SURE CLONES CONSISTENT WITH ABOVE NUMBERING) #########
G.wide$parents <- c(0, rep(1, dim(G.wide)[1] - 1))



# Then get the frequency data. (Use ?get_evofreq for options)
freq_frame <- get_evofreq(G.wide[,seq(1,dim(G.wide)[2]-2)], G.wide$clones, 
                          G.wide$parents, clone_cmap = "YIGnBu")

# Run function to set custom color scheme
source("~/scCode/evoFreq_fixColor.R")
freq_frame <- evoFreq_fixColor_clone(freq_frame)

# Create the plot
evo_freq_p <- plot_evofreq(freq_frame) + ylab("Relative Pop. Size") + ggtitle("Clonal Dynamics") + xlab(NULL) + 
  theme_bw() + theme(axis.text.x = element_text(face = "bold", family = "Arial", size = 12), 
                     axis.title.y = element_text(face = "bold", family = "Arial", size = 12),
                     plot.title = element_text(face = "bold", family = "Arial", size = 14))
print(evo_freq_p)

table(evo_freq_p$data$clone_id, evo_freq_p$data$plot_color)




# Do for cluster now
G.table <- table(G$Sample, G$Clus)
G.mat <- matrix(G.table, ncol = ncol(G.table), dimnames = dimnames(G.table))
G.wide <- as.data.frame(t(100*G.mat/rowSums(G.mat)))

####### Set initial frequency ########
G.wide <- cbind(rep(.01, dim(G.wide)[1]), G.wide)

colnames(G.wide) <- c("", "Sample 1", "Sample 2")

####### Set parents (here all arise separately) #########
G.wide$parents <- rep(0, dim(G.wide)[1])
G.wide$clones <- colnames(G.mat)


# Then get the frequency data. (Use ?get_evofreq for options)
freq_frame <- get_evofreq(G.wide[,seq(1,dim(G.wide)[2]-2)], G.wide$clones, 
                          G.wide$parents, clone_cmap = "jet")

freq_frame <- evoFreq_fixColor_clus(freq_frame)

# Create the plot
clus_freq_p <- plot_evofreq(freq_frame) + ylab("Relative Pop. Size") + ggtitle("Cluster Dynamics") + xlab(NULL) +
  theme_bw() + theme(axis.text.x = element_text(face = "bold", family = "Arial", size = 12), 
                     axis.title.y = element_text(face = "bold", family = "Arial", size = 12),
                     plot.title = element_text(face = "bold", family = "Arial", size = 14))
print(clus_freq_p)

#pdf(paste0(saveDir, "patient_G_clone_and_cluster.pdf"))

table(clus_freq_p$data$clone_id, clus_freq_p$data$plot_color)

####### TAKE SCREENSHOT FOR NOW...PDF WASN'T WORKING FOR ME
grid.arrange(evo_freq_p, clus_freq_p, ncol = 1)

#dev.off()

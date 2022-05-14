# Evo freq for mitoclone
library(EvoFreq)
library(gridExtra)
library(colormap)
library(viridis)

rm(list = ls())

# NOTE: This has the correct phylogeny for the clones present according to mitoclone for patient E

# Set dir to save to and from which we pull countTables
dir <- '/Users/Brian/scCode/mitoclone/'
saveDir <- paste0(dir, 'evoFreq_04042022_newParams/')
importDir <- paste0(dir, 'sequential_clus_clone_03272022_newParams/')

#######   SET PATIENT #############
currentSample <- "patientE"
E <- read.csv(paste0(importDir, currentSample, '_cluster_and_clone_info.csv'))


# Clone
E.table <- table(E$Sample, E$Clone)
E.mat <- matrix(E.table, ncol = ncol(E.table), dimnames = dimnames(E.table))
E.wide <- as.data.frame(t(100*E.mat/rowSums(E.mat)))

#######  SET INITIAL COMPOSITION (set initial to 0 if a clone stems from another clone) ##########
E.wide <- cbind(c(.01,rep(0, dim(E.wide)[1]-1)), E.wide)

colnames(E.wide) <- c("", "Sample 1", "Sample 2")

###### NUMBER CLONES FROM 1 TO N. NOTE: THIS WILL AFFECT THE PARENTS YOU SET
E.wide$clones <- c(1:dim(E.wide)[1])
########  CHANGE PARENTS TO FIT MITOCLONE PHYLOGENY (MAKE SURE CLONES CONSISTENT WITH ABOVE NUMBERING) #########
E.wide$parents <- c(0, rep(1, dim(E.wide)[1] - 1))



# Then get the frequency data. (Use ?get_evofreq for options)
freq_frame <- get_evofreq(E.wide[,seq(1,dim(E.wide)[2]-2)], E.wide$clones, 
                          E.wide$parents, clone_cmap = "YIGnBu")

# Run function to set custom color scheme
source("~/scCode/evoFreq_fixColor.R")
freq_frame <- evoFreq_fixColor_clone(freq_frame)

# Create the plot
evo_freq_p <- plot_evofreq(freq_frame) + ylab("Relative Pop. Size") + ggtitle("Clonal Dynamics") + xlab(NULL) + 
  theme_bw() + theme(axis.text.x = element_text(face = "bold", family = "Arial", size = 12), 
                     axis.title.y = element_text(face = "bold", family = "Arial", size = 12),
                     plot.title = element_text(face = "bold", family = "Arial", size = 14))
print(evo_freq_p)

#table(evo_freq_p$data$clone_id, evo_freq_p$data$plot_color)


# Do for cluster now
E.table <- table(E$Sample, E$Clus)
E.mat <- matrix(E.table, ncol = ncol(E.table), dimnames = dimnames(E.table))
E.wide <- as.data.frame(t(100*E.mat/rowSums(E.mat)))

####### Set initial frequency ########
E.wide <- cbind(rep(.01, dim(E.wide)[1]), E.wide)

colnames(E.wide) <- c("", "Sample 1", "Sample 2")

####### Set parents (here all arise separately) #########
E.wide$parents <- rep(0, dim(E.wide)[1])
E.wide$clones <- colnames(E.mat)


# Then get the frequency data. (Use ?get_evofreq for options)
freq_frame <- get_evofreq(E.wide[,seq(1,dim(E.wide)[2]-2)], E.wide$clones, 
                          E.wide$parents, clone_cmap = "jet")

# Set color scheme
freq_frame <- evoFreq_fixColor_clus(freq_frame)

# Create the plot
clus_freq_p <- plot_evofreq(freq_frame) + ylab("Relative Pop. Size") + ggtitle("Cluster Dynamics") + xlab(NULL) +
  theme_bw() + theme(axis.text.x = element_text(face = "bold", family = "Arial", size = 12), 
                     axis.title.y = element_text(face = "bold", family = "Arial", size = 12),
                     plot.title = element_text(face = "bold", family = "Arial", size = 14))

print(clus_freq_p)

#pdf(paste0(saveDir, "patient_E_clone_and_cluster.pdf"))

####### TAKE SCREENSHOT FOR NOW...PDF WASN'T WORKING FOR ME
grid.arrange(evo_freq_p, clus_freq_p, ncol = 1)


table(clus_freq_p[["data"]][["clone_id"]], clus_freq_p[["data"]][["plot_color"]])

#dev.off()


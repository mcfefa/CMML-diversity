# Evo freq for mitoclone
library(EvoFreq)
library(gridExtra)
library(colormap)
library(viridis)

rm(list = ls())

# NOTE: This has the correct phylogeny for the clones present according to mitoclone for patient F

# Set dir to save to and from which we pull countTables
dir <- '/Users/Brian/scCode/mitoclone/'
saveDir <- paste0(dir, 'evoFreq_04042022_newParams/')
importDir <- paste0(dir, 'sequential_clus_clone_03272022_newParams/')

#######   SET PATIENT #############
currentSample <- "patientF"
F <- read.csv(paste0(importDir, currentSample, '_cluster_and_clone_info.csv'))


# Clone
F.table <- table(F$Sample, F$Clone)
F.mat <- matrix(F.table, ncol = ncol(F.table), dimnames = dimnames(F.table))
F.wide <- as.data.frame(t(100*F.mat/rowSums(F.mat)))

#######  SET INITIAL COMPOSITION (set initial to 0 if a clone stems from another clone) ##########
F.wide <- cbind(c(.01,rep(0, dim(F.wide)[1]-1)), F.wide)

colnames(F.wide) <- c("", "Sample 1", "Sample 2")

###### NUMBER CLONES FROM 1 TO N. NOTE: THIS WILL AFFECT THE PARENTS YOU SET
F.wide$clones <- c(1:dim(F.wide)[1])
########  CHANGE PARENTS TO FIT MITOCLONE PHYLOGENY (MAKE SURE CLONES CONSISTENT WITH ABOVE NUMBERING) #########
F.wide$parents <- c(0, rep(1, dim(F.wide)[1] - 1))



# Then get the frequency data. (Use ?get_evofreq for options)
freq_frame <- get_evofreq(F.wide[,seq(1,dim(F.wide)[2]-2)], F.wide$clones, 
                          F.wide$parents, clone_cmap = "YIGnBu")

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
F.table <- table(F$Sample, F$Clus)
F.mat <- matrix(F.table, ncol = ncol(F.table), dimnames = dimnames(F.table))
F.wide <- as.data.frame(t(100*F.mat/rowSums(F.mat)))

####### Set initial frequency ########
F.wide <- cbind(rep(.01, dim(F.wide)[1]), F.wide)

colnames(F.wide) <- c("", "Sample 1", "Sample 2")

####### Set parents (here all arise separately) #########
F.wide$parents <- rep(0, dim(F.wide)[1])
F.wide$clones <- colnames(F.mat)


# Then get the frequency data. (Use ?get_evofreq for options)
freq_frame <- get_evofreq(F.wide[,seq(1,dim(F.wide)[2]-2)], F.wide$clones, 
                          F.wide$parents, clone_cmap = "jet")

# Set color scheme
freq_frame <- evoFreq_fixColor_clus(freq_frame)

# Create the plot
clus_freq_p <- plot_evofreq(freq_frame) + ylab("Relative Pop. Size") + ggtitle("Cluster Dynamics") + xlab(NULL) +
  theme_bw() + theme(axis.text.x = element_text(face = "bold", family = "Arial", size = 12), 
                     axis.title.y = element_text(face = "bold", family = "Arial", size = 12),
                     plot.title = element_text(face = "bold", family = "Arial", size = 14))

print(clus_freq_p)

#pdf(paste0(saveDir, "patient_F_clone_and_cluster.pdf"))

####### TAKE SCREENSHOT FOR NOW...PDF WASN'T WORKING FOR ME
grid.arrange(evo_freq_p, clus_freq_p, ncol = 1)


table(clus_freq_p[["data"]][["clone_id"]], clus_freq_p[["data"]][["plot_color"]])

#dev.off()

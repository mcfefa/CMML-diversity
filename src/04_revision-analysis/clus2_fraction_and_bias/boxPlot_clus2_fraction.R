library(Seurat)

rm(list = ls())

allSeurat <- readRDS("~/scCode/allSeurat_withsingleR_withBias_03092022.rds")

split_cluster_by_ind <- split(allSeurat$clusterResolution_0.05, paste0(allSeurat$bias, allSeurat$orig.ident))

clusFraction <- lapply(split_cluster_by_ind, function(x) {sum(x == "2")/sum(x != "-1")})

mep <- grep("MEP", names(clusFraction))
mono <- grep("Mono", names(clusFraction))
normal_like <- grep("Normal-like", names(clusFraction))
normal <- c(1:length(clusFraction))[-c(mep, mono, normal_like)]

df <- data.frame("clus2_fraction" = unlist(clusFraction), "Bias" = rep(NA, length(clusFraction)))

df$Bias[mep] <- "MEP"
df$Bias[mono] <- "Mono"
df$Bias[normal_like] <- "Normal_like"
df$Bias[normal] <- "Normal"

pdf("~/scCode/clus2_and_bias_boxplot_05112022.pdf")
ggplot(df, aes(factor(Bias, levels = c("MEP", "Mono", "Normal_like", "Normal")), y = clus2_fraction)) + 
  geom_boxplot(aes(fill = factor(Bias, levels = c("MEP", "Mono", "Normal_like", "Normal"))))+
  geom_point() + theme_bw() + labs(fill = "Bias")+
  xlab(c("MEP", "Mono", "Normal like", "Normal"))+
  ylab("Cluster 2 Fraction") + scale_fill_manual(values = c('#80C93B', '#31CEFF', "#D33032", "black"))+
  ggtitle("Cluster 2 fraction and bias") +
  theme(plot.title = element_text(face = "bold", size = 20), 
        axis.title.y = element_text(face = "bold", size = 17),
        axis.title.x = element_blank(),
        axis.text = element_text(face = "bold", size = 14),
        legend.text = element_text(face = "bold", size = 17),
        legend.title = element_text(face = "bold", size = 17))
dev.off()
 

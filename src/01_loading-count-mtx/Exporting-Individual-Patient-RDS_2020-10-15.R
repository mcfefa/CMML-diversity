## Run interactively on the Cluster
## qsub -I -q bigmemQ -l nodes=1:ppn=16 -l walltime=120:00:00 
## module load python/3.7.2
## module load R/3.6.0
## export R_LIBS_USER=$HOME/apps/R:$R_LIB_USER     #<---- if installing in home directory --- didn't do
## R

library(future)
library(bigmemory)
library(Seurat)
sessionInfo()

plan("multiprocess", workers=16)
setwd("/share/lab_padron/Meghan/scRNAseq/CMML/")

### already have an interactive job on 3-5, so should be delegated to one of the 512GB notes
options(future.globals.maxSize=536870912000) ###  512GB: 536870912000
# options(future.globals.maxSize=402653184000) ###  384GB: 536870912000
# #options(future.globals.maxSize=171966464000) ### for 164GB * 1024*1024
# #options(future.globals.maxSize=134217728000) ### for 128GB * 1024*1024

date <- "_2020-10-15"
savedir <- "/share/lab_altrock/MeghanCluster/BrianCluster/individual-CMML-rds/"

### CREATING INDIVIDUAL RDS FILES FOR EACH PATIENT FOR BRIAN TO INTEGRATE WITH PALANTIR EXAMPLE

# CMML SAMPLE IDS
# First 8: LTB3966, LTB4121, LTB5109 (6169), 4J003, 4K001, 4Q001, 5E001, 5H001
# Second 8: SF100109106293, SF100109111451, SF100109110236, SF14040100158,
#           SF14060200025, SF12062800475, SF14072200012, SF13061200056
# Third 8: 4S001, 2V001, SF14101000049, SF16112900158, 6AE001, 6AC001, 
#          6AD001, SF100109101914
# Fourth 8: SF12042500035, SF12091900043, SF12092600014, SF14031800065, 
#           SF14050700419, SF16026800045, SF16072200003, SF16112300029
# Fifth 8: SF13032800016, SF14110400108, SF1411140033, SF14092500135,
#          SF14061300036, SF14080400065, SF15010200008, SF13070900171 

CMML1.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML1_LTB3966/outs/filtered_feature_bc_matrix/")
CMML1 <- CreateSeuratObject(counts = CMML1.data, project = "LTB3966")
CMML1 <- RenameCells(object=CMML1, add.cell.id="LTB3966")
saveRDS(CMML1,paste(savedir,"CMML1_LTB3966",date,".rds",sep=""))

CMML2.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML2_LTB4121/outs/filtered_feature_bc_matrix/")
CMML2 <- CreateSeuratObject(counts = CMML2.data, project = "LTB4121")
CMML2 <- RenameCells(object=CMML2, add.cell.id="LTB4121")
saveRDS(CMML2,paste(savedir,"CMML2_LTB4121",date,".rds",sep=""))

CMML3.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML3_LTB6169/outs/filtered_feature_bc_matrix/")
CMML3 <- CreateSeuratObject(counts = CMML3.data, project = "LTB5109")
CMML3 <- RenameCells(object=CMML3, add.cell.id="LTB5109")
saveRDS(CMML3,paste(savedir,"CMML3_LTB5109",date,".rds",sep=""))

CMML4.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML4_4J003v2/outs/filtered_feature_bc_matrix/")
CMML4 <- CreateSeuratObject(counts = CMML4.data, project = "4J003")
CMML4 <- RenameCells(object=CMML4, add.cell.id="4J003")
saveRDS(CMML4,paste(savedir,"CMML4_4-J-003_postBMT",date,".rds",sep=""))

CMML5.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML5_4K001v2/outs/filtered_feature_bc_matrix/")
CMML5 <- CreateSeuratObject(counts = CMML5.data, project = "4K001")
CMML5 <- RenameCells(object=CMML5, add.cell.id="4K001")
saveRDS(CMML5,paste(savedir,"CMML5_4-K-001_HMA",date,".rds",sep=""))

CMML6.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML6_4Q001/outs/filtered_feature_bc_matrix/")
CMML6 <- CreateSeuratObject(counts = CMML6.data, project = "4Q001")
CMML6 <- RenameCells(object=CMML6, add.cell.id="4Q001")
saveRDS(CMML6,paste(savedir,"CMML6_4-Q-001_HMA",date,".rds",sep=""))

CMML7.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML7_5E001v2/outs/filtered_feature_bc_matrix/")
CMML7 <- CreateSeuratObject(counts = CMML7.data, project = "5E001")
CMML7 <- RenameCells(object=CMML7, add.cell.id="5E001")
saveRDS(CMML7,paste(savedir,"CMML7_5-E-001_HMA",date,".rds",sep=""))

CMML8.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML8_5H001v2/outs/filtered_feature_bc_matrix/")
CMML8 <- CreateSeuratObject(counts = CMML8.data, project = "5H001")
CMML8 <- RenameCells(object=CMML8, add.cell.id="5H001") 
saveRDS(CMML8,paste(savedir,"CMML8_5-H-001",date,".rds",sep=""))

CMML9.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML9_SF100109106293/outs/filtered_feature_bc_matrix/")
CMML9 <- CreateSeuratObject(counts = CMML9.data, project = "SF100109106293")
CMML9 <- RenameCells(object=CMML9, add.cell.id="SF100109106293")
saveRDS(CMML9,paste(savedir,"CMML9_SF-100109-106293",date,".rds",sep=""))

CMML10.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML10_SF100109111451/outs/filtered_feature_bc_matrix/")
CMML10 <- CreateSeuratObject(counts = CMML10.data, project = "SF100109111451")
CMML10 <- RenameCells(object=CMML10, add.cell.id="SF100109111451")
saveRDS(CMML10,paste(savedir,"CMML10_SF-100109-111451",date,".rds",sep=""))

CMML11.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML11_SF100109110236/outs/filtered_feature_bc_matrix/")
CMML11 <- CreateSeuratObject(counts = CMML11.data, project = "SF100109110236")
CMML11 <- RenameCells(object=CMML11, add.cell.id="SF100109110236")
saveRDS(CMML11,paste(savedir,"CMML11_SF-100109-110236",date,".rds",sep=""))

CMML12.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML12_SF14040100158/outs/filtered_feature_bc_matrix/")
CMML12 <- CreateSeuratObject(counts = CMML12.data, project = "SF14040100158")
CMML12 <- RenameCells(object=CMML12, add.cell.id="SF14040100158")
saveRDS(CMML12,paste(savedir,"CMML12_SF-140401-00158",date,".rds",sep=""))

CMML13.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML13_SF14060200025/outs/filtered_feature_bc_matrix/")
CMML13 <- CreateSeuratObject(counts = CMML13.data, project = "SF14060200025")
CMML13 <- RenameCells(object=CMML13, add.cell.id="SF14060200025")
saveRDS(CMML13,paste(savedir,"CMML13_SF-140602-00025",date,".rds",sep=""))

CMML14.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML14_SF12062800475/outs/filtered_feature_bc_matrix/")
CMML14 <- CreateSeuratObject(counts = CMML14.data, project = "SF12062800475")
CMML14 <- RenameCells(object=CMML14, add.cell.id="SF12062800475")
saveRDS(CMML14,paste(savedir,"CMML14_SF-120628-00475",date,".rds",sep=""))

CMML15.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML15_SF14072200012/outs/filtered_feature_bc_matrix/")
CMML15 <- CreateSeuratObject(counts = CMML15.data, project = "SF14072200012")
CMML15 <- RenameCells(object=CMML15, add.cell.id="SF14072200012")
saveRDS(CMML15,paste(savedir,"CMML15_SF-140722-00012",date,".rds",sep=""))

CMML16.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML16_SF13061200056/outs/filtered_feature_bc_matrix/")
CMML16 <- CreateSeuratObject(counts = CMML16.data, project = "SF13061200056")
CMML16 <- RenameCells(object=CMML16, add.cell.id="SF13061200056") 
saveRDS(CMML16,paste(savedir,"CMML16_SF-130612-00056",date,".rds",sep=""))

CMML17.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML17_4S001/outs/filtered_feature_bc_matrix/")
CMML17 <- CreateSeuratObject(counts = CMML17.data, project = "4S001")
CMML17 <- RenameCells(object=CMML17, add.cell.id="4S001")
saveRDS(CMML17,paste(savedir,"CMML17_4-S-001_HMA",date,".rds",sep=""))

CMML18.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML18_2V001/outs/filtered_feature_bc_matrix/")
CMML18 <- CreateSeuratObject(counts = CMML18.data, project = "2V001")
CMML18 <- RenameCells(object=CMML18, add.cell.id="2V001")
saveRDS(CMML18,paste(savedir,"CMML18_2-V-001",date,".rds",sep=""))

CMML19.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML19_SF14101000049/outs/filtered_feature_bc_matrix/")
CMML19 <- CreateSeuratObject(counts = CMML19.data, project = "SF14101000049")
CMML19 <- RenameCells(object=CMML19, add.cell.id="SF14101000049")
saveRDS(CMML19,paste(savedir,"CMML19_SF-141010-00049",date,".rds",sep=""))

CMML20.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML20_SF16112900158/outs/filtered_feature_bc_matrix/")
CMML20 <- CreateSeuratObject(counts = CMML20.data, project = "SF16112900158")
CMML20 <- RenameCells(object=CMML20, add.cell.id="SF16112900158")
saveRDS(CMML20,paste(savedir,"CMML20_SF-161129-00158",date,".rds",sep=""))

CMML21.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML21_6AE001/outs/filtered_feature_bc_matrix/")
CMML21 <- CreateSeuratObject(counts = CMML21.data, project = "6AE001")
CMML21 <- RenameCells(object=CMML21, add.cell.id="6AE001")
saveRDS(CMML21,paste(savedir,"CMML21_6-AE-001",date,".rds",sep=""))

CMML22.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML22_6AC001/outs/filtered_feature_bc_matrix/")
CMML22 <- CreateSeuratObject(counts = CMML22.data, project = "6AC001")
CMML22 <- RenameCells(object=CMML22, add.cell.id="6AC001")
saveRDS(CMML22,paste(savedir,"CMML22_6-AC-001",date,".rds",sep=""))

CMML23.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML23_6AD001/outs/filtered_feature_bc_matrix/")
CMML23 <- CreateSeuratObject(counts = CMML23.data, project = "6AD001")
CMML23 <- RenameCells(object=CMML23, add.cell.id="6AD001")

CMML24.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML24_SF100109101914/outs/filtered_feature_bc_matrix/")
CMML24 <- CreateSeuratObject(counts = CMML24.data, project = "SF100109101914")
CMML24 <- RenameCells(object=CMML24, add.cell.id="SF100109101914")
saveRDS(CMML24,paste(savedir,"CMML24_SF-100109-101914",date,".rds",sep=""))

CMML25.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML25_SF12042500035/outs/filtered_feature_bc_matrix/")
CMML25 <- CreateSeuratObject(counts=CMML25.data, project="SF12042500035")
CMML25 <- RenameCells(object=CMML25, add.cell.id="SF12042500035")
saveRDS(CMML25,paste(savedir,"CMML25_SF-120425-00035",date,".rds",sep=""))

## note: this patient SF-120919-00043 is actually the same patient at the same exact time point as 6-AD-001
CMML26.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML26_SF12091900043/outs/filtered_feature_bc_matrix/")
CMML26 <- CreateSeuratObject(counts=CMML26.data, project="SF12091900043")
CMML26 <- RenameCells(object=CMML26, add.cell.id="6AD001")

CMMLpt_6AD001 <- merge(x=CMML23, y=CMML26, merge.data=TRUE, project = "6AD001Merged")
saveRDS(CMMLpt_6AD001,paste(savedir,"CMML23+26_6-AD-001",date,".rds",sep=""))

CMML27.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML27_SF12092600014/outs/filtered_feature_bc_matrix/")
CMML27 <- CreateSeuratObject(counts=CMML27.data, project="SF12092600014")
CMML27 <- RenameCells(object=CMML27, add.cell.id="SF12092600014")
saveRDS(CMML27,paste(savedir,"CMML27_SF-120926-00014",date,".rds",sep=""))

CMML28.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML28_SF14031800065/outs/filtered_feature_bc_matrix/")
CMML28 <- CreateSeuratObject(counts=CMML28.data, project="SF14031800065")
CMML28 <- RenameCells(object=CMML28, add.cell.id="SF14031800065")
saveRDS(CMML28,paste(savedir,"CMML28_SF-140318-00065",date,".rds",sep=""))

CMML29.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML29_SF14050700419/outs/filtered_feature_bc_matrix/")
CMML29 <- CreateSeuratObject(counts=CMML29.data, project="SF14050700419")
CMML29 <- RenameCells(object=CMML29, add.cell.id="SF14050700419")
saveRDS(CMML29,paste(savedir,"CMML29_SF-140507-00419",date,".rds",sep=""))

CMML30.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML30_SF16026800045/outs/filtered_feature_bc_matrix/")
CMML30 <- CreateSeuratObject(counts=CMML30.data, project="SF16026800045")
CMML30 <- RenameCells(object=CMML30, add.cell.id="SF16026800045")
saveRDS(CMML30,paste(savedir,"CMML30_SF-160268-00045",date,".rds",sep=""))

CMML31.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML31_SF16072200003/outs/filtered_feature_bc_matrix/")
CMML31 <- CreateSeuratObject(counts=CMML31.data, project="SF16072200003")
CMML31 <- RenameCells(object=CMML31, add.cell.id="SF16072200003")
saveRDS(CMML31,paste(savedir,"CMML31_SF-160722-00003",date,".rds",sep=""))

CMML32.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML32_SF16112300029/outs/filtered_feature_bc_matrix/")
CMML32 <- CreateSeuratObject(counts=CMML32.data, project="SF16112300029")
CMML32 <- RenameCells(object=CMML32, add.cell.id="SF16112300029")
saveRDS(CMML32,paste(savedir,"CMML32_SF-161123-00029",date,".rds",sep=""))

#### TREATED SAMPLES

CMML33.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML33_SF13032800016/outs/filtered_feature_bc_matrix/")
CMML33 <- CreateSeuratObject(counts=CMML33.data, project="SF13032800016")
CMML33 <- RenameCells(object=CMML33, add.cell.id="SF13032800016")
saveRDS(CMML33,paste(savedir,"CMML33_SF-130328-00016",date,".rds",sep=""))

CMML34.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML34_SF14110400108/outs/filtered_feature_bc_matrix/")
CMML34 <- CreateSeuratObject(counts=CMML34.data, project="SF14110400108")
CMML34 <- RenameCells(object=CMML34, add.cell.id="SF14110400108")
saveRDS(CMML34,paste(savedir,"CMML34_SF-141104-00108",date,".rds",sep=""))

CMML35.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML35_SF14111400033/outs/filtered_feature_bc_matrix/")
CMML35 <- CreateSeuratObject(counts=CMML35.data, project="SF14111400033")
CMML35 <- RenameCells(object=CMML35, add.cell.id="SF14111400033")
saveRDS(CMML35,paste(savedir,"CMML35_SF-141114-00033",date,".rds",sep=""))

CMML36.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML36_SF14092500135/outs/filtered_feature_bc_matrix/")
CMML36 <- CreateSeuratObject(counts=CMML36.data, project="SF14092500135")
CMML36 <- RenameCells(object=CMML36, add.cell.id="SF14092500135")
saveRDS(CMML36,paste(savedir,"CMML36_SF-140925-00135",date,".rds",sep=""))

CMML37.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML37_SF14061300036/outs/filtered_feature_bc_matrix/")
CMML37 <- CreateSeuratObject(counts=CMML37.data, project="SF14061300036")
CMML37 <- RenameCells(object=CMML37, add.cell.id="SF14061300036")
saveRDS(CMML37,paste(savedir,"CMML37_SF-140613-00036",date,".rds",sep=""))

CMML38.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML38_SF14080400065/outs/filtered_feature_bc_matrix/")
CMML38 <- CreateSeuratObject(counts=CMML38.data, project="SF14080400065")
CMML38 <- RenameCells(object=CMML38, add.cell.id="SF14080400065")
saveRDS(CMML38,paste(savedir,"CMML38_SF-140804-00065",date,".rds",sep=""))

CMML39.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML39_SF15010200008/outs/filtered_feature_bc_matrix/")
CMML39 <- CreateSeuratObject(counts=CMML39.data, project="SF15010200008")
CMML39 <- RenameCells(object=CMML39, add.cell.id="SF15010200008")
saveRDS(CMML39,paste(savedir,"CMML39_SF-150102-00008",date,".rds",sep=""))

CMML40.data <- Read10X(data.dir = "/share/lab_padron/Meghan/scRNAseq/CMML/CMML40_SF13070900171/outs/filtered_feature_bc_matrix/")
CMML40 <- CreateSeuratObject(counts=CMML40.data, project="SF13070900171")
CMML40 <- RenameCells(object=CMML40, add.cell.id="SF13070900171")
saveRDS(CMML40,paste(savedir,"CMML40_SF-130709-00171",date,".rds",sep=""))

#<--------

## merged without cutoffs applied

savedir <- "/share/lab_altrock/MeghanCluster/BrianCluster/"

CMML_A <- merge(x=CMML1, y=CMML2, merge.data=TRUE, project = "EarlyCohort")
CMML_B <- merge(x=CMML_A, y=CMML3, merge.data=TRUE, project = "EarlyCohort")
CMML_C <- merge(x=CMML_B, y=CMML4, merge.data=TRUE, project = "EarlyCohort")
CMML_D <- merge(x=CMML_C, y=CMML5, merge.data=TRUE, project = "EarlyCohort")
CMML_E <- merge(x=CMML_D, y=CMML6, merge.data=TRUE, project = "EarlyCohort")
CMML_F <- merge(x=CMML_E, y=CMML7, merge.data=TRUE, project = "EarlyCohort")
CMML_G <- merge(x=CMML_F, y=CMML8, merge.data=TRUE, project = "EarlyCohort")
CMML_H <- merge(x=CMML_G, y=CMML9, merge.data=TRUE, project = "EarlyCohort")
CMML_I <- merge(x=CMML_H, y=CMML10, merge.data=TRUE, project = "EarlyCohort")
CMML_J <- merge(x=CMML_I, y=CMML11, merge.data=TRUE, project = "EarlyCohort")
CMML_K <- merge(x=CMML_J, y=CMML12, merge.data=TRUE, project = "EarlyCohort")
CMML_L <- merge(x=CMML_K, y=CMML13, merge.data=TRUE, project = "EarlyCohort")
CMML_M <- merge(x=CMML_L, y=CMML14, merge.data=TRUE, project = "EarlyCohort")
CMML_N <- merge(x=CMML_M, y=CMML15, merge.data=TRUE, project = "EarlyCohort")
CMML_O <- merge(x=CMML_N, y=CMML16, merge.data=TRUE, project = "EarlyCohort")
CMML_P <- merge(x=CMML_O, y=CMML17, merge.data=TRUE, project = "EarlyCohort")
CMML_Q <- merge(x=CMML_P, y=CMML18, merge.data=TRUE, project = "EarlyCohort")
CMML_R <- merge(x=CMML_Q, y=CMML19, merge.data=TRUE, project = "EarlyCohort")
CMML_S <- merge(x=CMML_R, y=CMML21, merge.data=TRUE, project = "EarlyCohort")
CMML_T <- merge(x=CMML_S, y=CMML22, merge.data=TRUE, project = "EarlyCohort")
CMML_U <- merge(x=CMML_T, y=CMMLpt_6AD001, merge.data=TRUE, project = "EarlyCohort")
CMML_V <- merge(x=CMML_U, y=CMML24, merge.data=TRUE, project = "EarlyCohort")
CMML_W <- merge(x=CMML_V, y=CMML25, merge.data=TRUE, project = "EarlyCohort")
CMML_X <- merge(x=CMML_W, y=CMML27, merge.data=TRUE, project = "EarlyCohort")
CMML_Y <- merge(x=CMML_X, y=CMML28, merge.data=TRUE, project = "EarlyCohort")
CMML_Z <- merge(x=CMML_Y, y=CMML29, merge.data=TRUE, project = "EarlyCohort")
CMML_AA <- merge(x=CMML_Z, y=CMML30, merge.data=TRUE, project = "EarlyCohort")
CMML_AB <- merge(x=CMML_AA, y=CMML31, merge.data=TRUE, project = "EarlyCohort")
CMML_AC <- merge(x=CMML_AB, y=CMML32, merge.data=TRUE, project = "EarlyCohort")

saveRDS(CMML_AC,paste(savedir,"CMML1-32_Early-Cohort_Raw-No-Cutoffs-Applied",date,".rds",sep=""))

CMML_Tx1 <- merge(x=CMML33, y=CMML34, merge.data=TRUE, project = "TxCohort")
CMML_Tx2 <- merge(x=CMML_Tx1, y=CMML35, merge.data=TRUE, project = "TxCohort")
CMML_Tx3 <- merge(x=CMML_Tx2, y=CMML36, merge.data=TRUE, project = "TxCohort")
CMML_Tx4 <- merge(x=CMML_Tx3, y=CMML37, merge.data=TRUE, project = "TxCohort")
CMML_Tx5 <- merge(x=CMML_Tx4, y=CMML38, merge.data=TRUE, project = "TxCohort")
CMML_Tx6 <- merge(x=CMML_Tx5, y=CMML39, merge.data=TRUE, project = "TxCohort")
CMML_Tx7 <- merge(x=CMML_Tx6, y=CMML40, merge.data=TRUE, project = "TxCohort")

saveRDS(CMML_Tx7,paste(savedir,"CMML33-40_Treated-Cohort_Raw-No-Cutoffs-Applied",date,".rds",sep=""))

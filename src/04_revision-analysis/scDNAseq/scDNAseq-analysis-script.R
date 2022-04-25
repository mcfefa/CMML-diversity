################################################################################
###
### scDNAseq analysis pipeline followed along with Miles et al Nature (2020)
###       https://bowmanr.github.io/scDNA_myeloid/
###       https://github.com/bowmanr/scDNA_myeloid
###
################################################################################

library(rhdf5)
library(dplyr)
library(tidyr)

pt1file <- "/Users/ferrallm/Dropbox (UFL)/papers-in-progress/CMML-scRNAseq-Paper/analysis/scDNAseq/Tapestri_Output_Files/5-M-001.dna+protein.h5"
pt2file <- "/Users/ferrallm/Dropbox (UFL)/papers-in-progress/CMML-scRNAseq-Paper/analysis/scDNAseq/Tapestri_Output_Files/4-I-001.dna+protein.h5"


###############################
### HDF5 and Loom Input
###############################
# Patient 1: 5-M-001 (high CRD, high CD120b)
protein_mat_5M<-h5read(file=pt1file,name="/assays/protein_read_counts/layers/read_counts")
rownames(protein_mat_5M) <- h5read(file=pt1file,name="/assays/protein_read_counts/ca/id")
colnames(protein_mat_5M)<-  h5read(file=pt1file,name="/assays/protein_read_counts/ra/barcode")

print(rownames(protein_mat_5M))
print((protein_mat_5M)[1:5,1:5])

# Patient 2: 4-I-001 (low CRD, low CD120b)
protein_mat_4I<-h5read(file=pt2file,name="/assays/protein_read_counts/layers/read_counts")
rownames(protein_mat_4I) <- h5read(file=pt2file,name="/assays/protein_read_counts/ca/id")
colnames(protein_mat_4I)<-  h5read(file=pt2file,name="/assays/protein_read_counts/ra/barcode")

print(rownames(protein_mat_4I))
print((protein_mat_4I)[1:5,1:5])

###############################
### Extracting variants 
###############################

# VAF cutoff using is that variant must be present in at least 1% of cell (per Miles et al)
VAF_cutoff <- 0.01

# Patient 1: 5-M-001 (high CRD, high CD120b)
NGT_5M<-h5read(file=pt1file,name="/assays/dna_variants/layers/NGT")
NGT_5M[NGT_5M==3]<-NA
VAF_select_5M<-which((rowSums(NGT_5M,na.rm=TRUE)/(ncol(NGT_5M)*2))>VAF_cutoff)
AF_5M<-h5read(file=pt1file,name="/assays/dna_variants/layers/AF",index=list(VAF_select_5M,NULL))
DP_5M<-h5read(file=pt1file,name="/assays/dna_variants/layers/DP",index=list(VAF_select_5M,NULL))
GQ_5M<-h5read(file=pt1file,name="/assays/dna_variants/layers/GQ",index=list(VAF_select_5M,NULL))
NGTlim_5M<-h5read(file=pt1file,name="/assays/dna_variants/layers/NGT",index=list(VAF_select_5M,NULL))
NGTlim_5M[NGTlim_5M==3]<-NA
variants_5M<-h5read(file=pt1file,name="/assays/dna_variants/ca/id",index=list(VAF_select_5M))
cell_barcodes_5M <-h5read(file=pt1file,name="/assays/dna_variants/ra/barcode")
colnames(NGTlim_5M) <-cell_barcodes_5M
rownames(NGTlim_5M) <- variants_5M

print(NGTlim_5M[1:5,1:5])
dim(NGTlim_5M) # [1] 607 393

# Patient 2: 4-I-001 (low CRD, low CD120b)
NGT_4I<-h5read(file=pt2file,name="/assays/dna_variants/layers/NGT")
NGT_4I[NGT_4I==3]<-NA
VAF_select_4I<-which((rowSums(NGT_4I,na.rm=TRUE)/(ncol(NGT_4I)*2))>VAF_cutoff)
AF_4I<-h5read(file=pt2file,name="/assays/dna_variants/layers/AF",index=list(VAF_select_4I,NULL))
DP_4I<-h5read(file=pt2file,name="/assays/dna_variants/layers/DP",index=list(VAF_select_4I,NULL))
GQ_4I<-h5read(file=pt2file,name="/assays/dna_variants/layers/GQ",index=list(VAF_select_4I,NULL))
NGTlim_4I<-h5read(file=pt2file,name="/assays/dna_variants/layers/NGT",index=list(VAF_select_4I,NULL))
NGTlim_4I[NGTlim_4I==3]<-NA
variants_4I<-h5read(file=pt2file,name="/assays/dna_variants/ca/id",index=list(VAF_select_4I))
cell_barcodes_4I <-h5read(file=pt2file,name="/assays/dna_variants/ra/barcode")
colnames(NGTlim_4I) <-cell_barcodes_4I
rownames(NGTlim_4I) <- variants_4I

print(NGTlim_4I[1:5,1:5])
dim(NGTlim_4I) # [1]  261 7889


###############################
### Filtering variants 
###############################

## set filtering criteria
DP_cut=10 #read depth
AF_cut=25 #allele frequency cutoff
GQ_cut=30 #genotype quality cutoff
variant_presence_cutoff=50 #variant must be genotyped in greater than 50% of cells
cell_genotype_cutoff=50 # cell must possess genotype information for at least 50% of the variants of interest

## bind together long form AF, DP, GQ and NGT data

# Patient 1: 5-M-001 (high CRD, high CD120b)
filtered_long_5M<-data.frame(setNames(
  # produce long form allele frequency data
  data.frame(AF_5M, 
             "variants"=variants_5M), 
            c(all_of(cell_barcodes_5M),"variants")) %>%
    pivot_longer(cols=!c(variants), names_to="Cell", values_to="AF"),
  #produce long form allele depth data
  data.frame(DP_5M)%>%
    pivot_longer(cols=everything(),
                 names_to="Cell",
                 values_to="DP")%>%dplyr::select(DP),
  #produce long form genotype quality data
  data.frame(GQ_5M)%>%
    pivot_longer(cols=everything(),
                 names_to="Cell",
                 values_to="GQ")%>%dplyr::select(GQ),
  #produce long form genotype call data
  data.frame(NGTlim_5M)%>%
    pivot_longer(cols=everything(), names_to="Cell",
                 values_to="NGT")%>%dplyr::select(NGT)) %>%
  #filter DP and GQ
  filter(DP>DP_cut&
           GQ>GQ_cut)%>%
  #filter AF for each genotype call
  mutate(pass=case_when(
    NGT==1&(AF>AF_cut)&(AF<(100-AF_cut)) ~ "include",
    NGT==1&((AF<=AF_cut)|(AF>=(100-AF_cut))) ~ "exclude",
    NGT==2&AF>=(100-AF_cut) ~ "include",
    NGT==2&AF<(100-AF_cut) ~ "exclude",
    NGT==0&AF<=AF_cut ~ "include",
    NGT==0&AF>AF_cut ~ "exclude",
    TRUE ~"other"
  ))%>%
  filter(pass=="include")

dim(filtered_long_5M) # [1] 173158      7

# Patient 2: 4-I-001 (low CRD, low CD120b)
filtered_long_4I<-data.frame(setNames(
  # produce long form allele frequency data
  data.frame(AF_4I, 
             "variants"=variants_4I), 
  c(all_of(cell_barcodes_4I),"variants")) %>%
    pivot_longer(cols=!c(variants), names_to="Cell", values_to="AF"),
  #produce long form allele depth data
  data.frame(DP_4I)%>%
    pivot_longer(cols=everything(),
                 names_to="Cell",
                 values_to="DP")%>%dplyr::select(DP),
  #produce long form genotype quality data
  data.frame(GQ_4I)%>%
    pivot_longer(cols=everything(),
                 names_to="Cell",
                 values_to="GQ")%>%dplyr::select(GQ),
  #produce long form genotype call data
  data.frame(NGTlim_4I)%>%
    pivot_longer(cols=everything(), names_to="Cell",
                 values_to="NGT")%>%dplyr::select(NGT)) %>%
  #filter DP and GQ
  filter(DP>DP_cut&
           GQ>GQ_cut)%>%
  #filter AF for each genotype call
  mutate(pass=case_when(
    NGT==1&(AF>AF_cut)&(AF<(100-AF_cut)) ~ "include",
    NGT==1&((AF<=AF_cut)|(AF>=(100-AF_cut))) ~ "exclude",
    NGT==2&AF>=(100-AF_cut) ~ "include",
    NGT==2&AF<(100-AF_cut) ~ "exclude",
    NGT==0&AF<=AF_cut ~ "include",
    NGT==0&AF>AF_cut ~ "exclude",
    TRUE ~"other"
  ))%>%
  filter(pass=="include")

dim(filtered_long_4I) # dim(filtered_long_4I)

## here we check to make sure a variant is not NA (genotype 3) in >50% of cells
## we can also filter out variants that are likely SNPs and only show up as 
##    WT, het or homozygous  This assumes SNPs would never undergo allele 
##    dropout, which is pretty unlikely, so this filter is largely unused

# Patient 1: 5-M-001 (high CRD, high CD120b)
final_variants_5M <-filtered_long_5M%>%
  group_by(variants)%>%
  summarize(diversity=sum(c(0,1,2)%in%NGT),
            gt.mv=(length(NGT)/length(all_of(cell_barcodes_5M)))*100)%>%
  filter(#diversity>1&
    gt.mv>variant_presence_cutoff)%>%
  pull(variants)

length(final_variants_5M) # 512

# Patient 2: 4-I-001 (low CRD, low CD120b)
final_variants_4I <-filtered_long_4I%>%
  group_by(variants)%>%
  summarize(diversity=sum(c(0,1,2)%in%NGT),
            gt.mv=(length(NGT)/length(all_of(cell_barcodes_4I)))*100)%>%
  filter(#diversity>1&
    gt.mv>variant_presence_cutoff)%>%
  pull(variants)

length(final_variants_4I) # 150


## here we filter for cells that now contain genotype information for 
##     at least 50% of our curated set of variants

# Patient 1: 5-M-001 (high CRD, high CD120b)
final_cells_5M<-  filtered_long_5M%>%
  filter(variants%in%final_variants_5M)%>%
  group_by(Cell)%>%
  summarize(gt.mc=(length(NGT)/length(all_of(final_variants_5M)))*100)%>%
  filter(gt.mc>=cell_genotype_cutoff)%>%
  pull(Cell)

length(final_cells_5M) # 393

# Patient 2: 4-I-001 (low CRD, low CD120b)
final_cells_4I<-  filtered_long_4I%>%
  filter(variants%in%final_variants_4I)%>%
  group_by(Cell)%>%
  summarize(gt.mc=(length(NGT)/length(all_of(final_variants_4I)))*100)%>%
  filter(gt.mc>=cell_genotype_cutoff)%>%
  pull(Cell)

length(final_cells_4I) # 7889


## lastly we reconstruct a new NGT matrix of cell-genotype pairs that passed 
## the above filters

# Patient 1: 5-M-001 (high CRD, high CD120b)
final_NGT_5M<- filtered_long_5M%>%
  filter(Cell%in%final_cells_5M&variants%in%final_variants_5M)%>%
  pivot_wider(id_cols=Cell,names_from=variants,values_from=NGT)

dim(final_NGT_5M) # [1] 393 513

# Patient 2: 4-I-001 (low CRD, low CD120b)
final_NGT_4I<- filtered_long_4I%>%
  filter(Cell%in%final_cells_4I&variants%in%final_variants_4I)%>%
  pivot_wider(id_cols=Cell,names_from=variants,values_from=NGT)

dim(final_NGT_4I) # [1] 7889  151


###############################
### Annotating variants 
###############################
# BiocManager::install("VariantAnnotation")
library(VariantAnnotation)
library(GenomicRanges)
library(magrittr)
# install.packages("RMariaDB")
require(RMariaDB)
# BiocManager::install("plyranges")
library(plyranges)

# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
# library(BSgenome.Hsapiens.UCSC.hg38)
# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)
# BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

annotation_key <-read.csv("./src/04_revision-analysis/scDNAseq/Miles_annotation_key.csv")
hg19refseq_txdb <- makeTxDbFromUCSC(genome="hg19",
                                    transcript_ids=annotation_key$ccds_id,
                                    tablename="ccdsGene")

annotation_key%<>%inner_join(select(hg19refseq_txdb,
                                    keys=annotation_key$ccds_id,
                                    columns=c("TXID","TXNAME"),
                                    keytype = "TXNAME"),
                             by=c("ccds_id"="TXNAME"))%>%
  mutate(TXID=as.character(TXID))


# Patient 1: 5-M-001 (high CRD, high CD120b)
# banned <-read.csv("./data/banned_list.csv")
SNV_mat_5M<-data.frame(do.call(cbind,
                            h5read(file=pt1file,name="/assays/dna_variants/ca/",
                                   index=list(VAF_select_5M)))) %>%
  filter(id%in%final_variants_5M)%>%
  mutate(ALT=ifelse(ALT=="*","N",ALT))%>%
  mutate(CHROM=paste0("chr",CHROM))

# Patient 2: 4-I-001 (low CRD, low CD120b)
SNV_mat_4I<-data.frame(do.call(cbind,
                               h5read(file=pt2file,name="/assays/dna_variants/ca/",
                                      index=list(VAF_select_4I)))) %>%
  filter(id%in%final_variants_4I)%>%
  mutate(ALT=ifelse(ALT=="*","N",ALT))%>%
  mutate(CHROM=paste0("chr",CHROM))

#necessary for meaningful GRangess
SNV_mat_5M$REF<-as(SNV_mat_5M$REF, "DNAStringSet")
SNV_mat_4I$REF<-as(SNV_mat_4I$REF, "DNAStringSet")
SNV_mat_5M$ALT<-as(SNV_mat_5M$ALT, "DNAStringSet")
SNV_mat_4I$ALT<-as(SNV_mat_4I$ALT, "DNAStringSet")

variant_gRange_5M<-makeGRangesFromDataFrame(SNV_mat_5M,
                                         seqnames.field = "CHROM",
                                         start.field="POS",
                                         end.field="POS",
                                         keep.extra.columns=TRUE)

variant_gRange_4I<-makeGRangesFromDataFrame(SNV_mat_4I,
                                         seqnames.field = "CHROM",
                                         start.field="POS",
                                         end.field="POS",
                                         keep.extra.columns=TRUE)

#necessary for downstream joining of
variant_gRange_5M$QUERYID<-1:length(variant_gRange_5M)
variant_gRange_4I$QUERYID<-1:length(variant_gRange_4I)

## identify and isolate non coding variants
# Patient 1: 5-M-001 (high CRD, high CD120b)
non_coding_variants_5M <- locateVariants(variant_gRange_5M, 
                                      hg19refseq_txdb,
                                      AllVariants())%>%
  data.frame()%>%
  filter(as.character(LOCATION)!="coding")%>%
  inner_join(variant_gRange_5M,by="QUERYID",copy=TRUE)

# Warning message:
#   In valid.GenomicRanges.seqinfo(x, suggest.trim = TRUE) :
#   GRanges object contains 106 out-of-bound ranges located on sequences 3, 4, 7, 17, 19, 22, 24, 26, and 30. Note that ranges located on
# a sequence whose length is unknown (NA) or on a circular sequence are not considered out-of-bound (use seqlengths() and isCircular()
#                                                                                                    to get the lengths and circularity flags of the underlying sequences). You can use trim() to trim these ranges. See
# ?`trim,GenomicRanges-method` for more information.

# Patient 2: 4-I-001 (low CRD, low CD120b)
non_coding_variants_4I <- locateVariants(variant_gRange_4I, 
                                         hg19refseq_txdb,
                                         AllVariants())%>%
  data.frame()%>%
  filter(as.character(LOCATION)!="coding")%>%
  inner_join(variant_gRange_4I,by="QUERYID",copy=TRUE)

# Warning message:
#   In valid.GenomicRanges.seqinfo(x, suggest.trim = TRUE) :
#   GRanges object contains 20 out-of-bound ranges located on sequences 3, 19, 24, and 26. Note that ranges located on a sequence whose
# length is unknown (NA) or on a circular sequence are not considered out-of-bound (use seqlengths() and isCircular() to get the lengths
#                                                                                   and circularity flags of the underlying sequences). You can use trim() to trim these ranges. See ?`trim,GenomicRanges-method` for more
# information.

## identify and isolate  coding variants
# Patient 1: 5-M-001 (high CRD, high CD120b)
coding_variants_5M  <-  predictCoding(variant_gRange_5M, 
                                   hg19refseq_txdb, 
                                   seqSource=Hsapiens,
                                   varAllele=variant_gRange_5M$ALT)%>%
  data.frame()

# Warning messages:
#   1: In valid.GenomicRanges.seqinfo(x, suggest.trim = TRUE) :
#   GRanges object contains 106 out-of-bound ranges located on sequences 3, 4, 7, 17, 19, 22, 24, 26, and 30. Note that ranges located on
# a sequence whose length is unknown (NA) or on a circular sequence are not considered out-of-bound (use seqlengths() and isCircular()
#                                                                                                    to get the lengths and circularity flags of the underlying sequences). You can use trim() to trim these ranges. See
# ?`trim,GenomicRanges-method` for more information.
# 2: In .predictCodingGRangesList(query, cache[["cdsbytx"]], seqSource,  :
#                                   'varAllele' values with 'N', '.', '+' or '-' were not translated


# Patient 2: 4-I-001 (low CRD, low CD120b)
coding_variants_4I  <-  predictCoding(variant_gRange_4I, 
                                   hg19refseq_txdb, 
                                   seqSource=Hsapiens,
                                   varAllele=variant_gRange_4I$ALT)%>%
  data.frame()

# Warning messages:
#   1: In valid.GenomicRanges.seqinfo(x, suggest.trim = TRUE) :
#   GRanges object contains 20 out-of-bound ranges located on sequences 3, 19, 24, and 26. Note that ranges located on a sequence whose
# length is unknown (NA) or on a circular sequence are not considered out-of-bound (use seqlengths() and isCircular() to get the lengths
#                                                                                   and circularity flags of the underlying sequences). You can use trim() to trim these ranges. See ?`trim,GenomicRanges-method` for more
# information.
# 2: In .predictCodingGRangesList(query, cache[["cdsbytx"]], seqSource,  :
#                                   'varAllele' values with 'N', '.', '+' or '-' were not translated


## Bind it all together into one big table.          
# Patient 1: 5-M-001 (high CRD, high CD120b)
out_5M <- bind_rows(non_coding_variants_5M,coding_variants_5M) %>%
  inner_join(annotation_key)%>%
  mutate(AA=ifelse(!is.na(REFAA),
                   paste0(gene_name,".",REFAA,PROTEINLOC,VARAA),
                   paste0(gene_name,".intronic")))%>%
  dplyr::select(id,AA)

# Patient 2: 4-I-001 (low CRD, low CD120b)
out_4I <- bind_rows(non_coding_variants_4I,coding_variants_4I) %>%
  inner_join(annotation_key)%>%
  mutate(AA=ifelse(!is.na(REFAA),
                   paste0(gene_name,".",REFAA,PROTEINLOC,VARAA),
                   paste0(gene_name,".intronic")))%>%
  dplyr::select(id,AA)


## append Bulk VAF for reference in future cutoffs and allele selection
# Patient 1: 5-M-001 (high CRD, high CD120b)
final_mutation_info_5M<-data.frame(out_5M,
                                "Bulk_VAF"=colSums(final_NGT_5M[,out_5M$id], na.rm=TRUE)/
                                  (nrow(final_NGT_5M)*2)*100) 

print(head(final_mutation_info_5M))

# Patient 2: 4-I-001 (low CRD, low CD120b)
final_mutation_info_4I<-data.frame(out_4I,
                                   "Bulk_VAF"=colSums(final_NGT_4I[,out_4I$id], na.rm=TRUE)/
                                     (nrow(final_NGT_4I)*2)*100) 

print(head(final_mutation_info_4I))

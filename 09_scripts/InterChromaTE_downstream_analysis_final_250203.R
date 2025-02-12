### 03/02/2025
### Ewan Harney

#################################################################
# All downstream analyses using sequence data for InterChromaTE
#
# 1. DEseq2 for RNA 
# 2. DEseq2 for ATAC 
# 3. Combining results and creating complete data table   
# 4. Euler plots to visualise set analysis
# 5. F3 expression versus accessibility (rna-seq vs atac-seq) 
# 6. Enrichment Analysis of genes that were both DE and DA 
# 7. Compare gene expression across F3 & F6 (transgenerational)
# 8. Associations between TEs and expression & accessibility
#
###################################################################

# A list of all packages loaded is provided here
#
# Packages are also loaded in the first section that they are used
# but they may be required for later sections.
# Some of these packages also have dependencies (not listed here).
library(tximport)
library(DESeq2)
library(org.Dm.eg.db)
library(EnhancedVolcano)
library(edgeR)
library(limma)
library(dplyr)
library(eulerr)
library(ggplot2)
library(ggblend)
library(ggrepel)
library(clusterProfiler)
library(rrvgo)
library(chisq.posthoc.test)

################################################################
################################################################
# 1. DEseq2 for RNA 
################################################################
################################################################

# clear environment
rm(list=ls(all=TRUE))

# Load required packages
library(tximport)
library(DESeq2)
library(org.Dm.eg.db)
library(EnhancedVolcano)

# Set working directory for whole script: SuppMat 
workingDir = "C:/../../SuppMat/";
setwd(workingDir); 
getwd();

# read in transcript & gene annotation information
tx2gene<-read.table("interchromate_meta/fbtn_fbgn.tsv",header=TRUE, sep = '\t')
# filenames of Kallisto output
Kfiles1<-read.table("interchromate_meta/Kallisto_filenames.txt",header=FALSE, sep = '\t')
# meta data (which is the same for RNAseq and ATACseq)
coldata<-read.table("interchromate_meta/Sample_info.txt",header=TRUE,row.names=1)

################################################################
# 1.1 Import kallisto results directly and transform them

# Need to provide precise directory location of Kallisto output
dir1="C:/Users/ewanh/Dropbox/Barcelona_IBE/WP1/SuppMat/interchromate_kallisto"

# prep
files1 <- file.path(dir1, Kfiles1$V1, "abundance.tsv")
names(files1) <- sort(rownames(coldata))
# 1. Control v Heatshock F3 Akaa
txi.k1.tsv <- tximport(files1[1:6], type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
# 2. Control v Heatshock F3 Manz
txi.k2.tsv <- tximport(files1[13:18], type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
# 3. Control v Heatshock F6 Akaa
txi.k3.tsv <- tximport(files1[7:12], type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
# 4. Control v Heatshock F6 Manz
txi.k4.tsv <- tximport(files1[19:24], type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
# 5. F3 Population controls
txi.k5.tsv <- tximport(files1[c(1:3,13:15)], type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
# 6. F6 Population controls
txi.k6.tsv <- tximport(files1[c(7:9,19:21)], type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
# 7. Akaa Generational controls
txi.k7.tsv <- tximport(files1[c(1:3,7:9)], type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
# 8. Manz Generational controls
txi.k8.tsv <- tximport(files1[c(13:15,19:21)], type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)

# Define Sample data
coldata$pop <- factor(coldata$pop)
coldata$gen <- factor(coldata$gen)
coldata$trt <- factor(coldata$trt)
coldata$pop_gen <- factor(paste(coldata$pop, coldata$gen, sep = "_"))
coldata$pop_trt <- factor(paste(coldata$pop, coldata$trt, sep = "_"))
coldata$gen_trt <- factor(paste(coldata$gen, coldata$trt, sep = "_"))
coldata$pop_gen_trt <- factor(paste(coldata$pop, coldata$gen, coldata$trt, sep = "_"))

################################################################
# 1.2 Run DEseq2 for different comparisons
#     Generate tables of results for sig DEGs
#     Volcano plots for Heat shock vs Ctrl

##### Akaa F3
# Specifying which sample data
F3Acoldata<-coldata[c(1:6),]
# import from txi object
dds_R_F3A <- DESeqDataSetFromTximport(txi.k1.tsv, colData = F3Acoldata, design = ~ trt)
# exclude genes with low counts
keep_R_F3A <- rowSums(counts(dds_R_F3A) >= 10) >= 3
dds_R_F3Ak <- dds_R_F3A[keep_R_F3A,]
# Run DESeq2
dds_R_F3Ak <- DESeq(dds_R_F3Ak)
# Apply lfcshrink transformation
R_F3A_res <- lfcShrink(dds=dds_R_F3Ak, coef=2, type="apeglm")
# summarise the DESeq2 object
summary(R_F3A_res)
# Select on;y those genes with a significant difference
R_F3A_res_sig<- subset(R_F3A_res, padj < 0.05)
R_F3A_T<-cbind(rownames(R_F3A_res_sig),data.frame(R_F3A_res_sig, row.names = NULL))
# Get gene symbols 
R_F3A_res$name = mapIds(org.Dm.eg.db,keys=rownames(R_F3A_res), column="SYMBOL", keytype="FLYBASE",multiVals="first")
# Enhance volcano plot
EnhancedVolcano(R_F3A_res, lab = R_F3A_res$name, x = 'log2FoldChange', y = 'pvalue', drawConnectors = TRUE, title = "Akaa F3", subtitle = "DE: HS vs Ctrl",legendPosition = "none")

##### Manz F3
F3Mcoldata<-coldata[c(13:18),]
dds_R_F3M <- DESeqDataSetFromTximport(txi.k2.tsv, colData = F3Mcoldata, design = ~ trt)
keep_R_F3M <- rowSums(counts(dds_R_F3M) >= 10) >= 3
dds_R_F3Mk <- dds_R_F3M[keep_R_F3M,]
dds_R_F3Mk <- DESeq(dds_R_F3Mk)
R_F3M_res <- lfcShrink(dds=dds_R_F3Mk, coef=2, type="apeglm")
summary(R_F3M_res)
R_F3M_res_sig<- subset(R_F3M_res, padj < 0.05)
R_F3M_T<-cbind(rownames(R_F3M_res_sig),data.frame(R_F3M_res_sig, row.names = NULL))
R_F3M_res$name = mapIds(org.Dm.eg.db,keys=rownames(R_F3M_res), column="SYMBOL", keytype="FLYBASE",multiVals="first")
EnhancedVolcano(R_F3M_res, lab = R_F3M_res$name, x = 'log2FoldChange', y = 'pvalue', drawConnectors = TRUE, title = "Manz F3", subtitle = "DE: HS vs Ctrl",legendPosition = "none")

##### Akaa F6
F6Acoldata<-coldata[c(7:12),]
dds_R_F6A <- DESeqDataSetFromTximport(txi.k3.tsv, colData = F6Acoldata, design = ~ trt)
keep_R_F6A <- rowSums(counts(dds_R_F6A) >= 10) >= 3
dds_R_F6Ak <- dds_R_F6A[keep_R_F6A,]
dds_R_F6Ak <- DESeq(dds_R_F6Ak)
R_F6A_res <- lfcShrink(dds=dds_R_F6Ak, coef=2, type="apeglm")
summary(R_F6A_res)
R_F6A_res_sig<- subset(R_F6A_res, padj < 0.05)
R_F6A_T<-cbind(rownames(R_F6A_res_sig),data.frame(R_F6A_res_sig, row.names = NULL))
R_F6A_res$name = mapIds(org.Dm.eg.db,keys=rownames(R_F6A_res), column="SYMBOL", keytype="FLYBASE",multiVals="first")
EnhancedVolcano(R_F6A_res, lab = R_F6A_res$name, x = 'log2FoldChange', y = 'pvalue', title = "Akaa F6", subtitle = "DE: HS vs Ctrl",legendPosition = "none")

##### Manz F6
F6Mcoldata<-coldata[c(19:24),]
dds_R_F6M <- DESeqDataSetFromTximport(txi.k4.tsv, colData = F6Mcoldata, design = ~ trt)
keep_R_F6M <- rowSums(counts(dds_R_F6M) >= 10) >= 3
dds_R_F6Mk <- dds_R_F6M[keep_R_F6M,]
dds_R_F6Mk <- DESeq(dds_R_F6Mk)
R_F6M_res <- lfcShrink(dds=dds_R_F6Mk, coef=2, type="apeglm")
summary(R_F6M_res)
R_F6M_res_sig<- subset(R_F6M_res, padj < 0.05)
R_F6M_T<-cbind(rownames(R_F6M_res_sig),data.frame(R_F6M_res_sig, row.names = NULL))
R_F6M_res$name = mapIds(org.Dm.eg.db,keys=rownames(R_F6M_res), column="SYMBOL", keytype="FLYBASE",multiVals="first")
EnhancedVolcano(R_F6M_res, lab = R_F6M_res$name, x = 'log2FoldChange', y = 'pvalue', title = "Manz F6", subtitle = "DE: HS vs Ctrl",legendPosition = "none")

##### Population comparison F3
F3Popcoldata<-coldata[c(1:3,13:15),]
dds_R_F3Pop <- DESeqDataSetFromTximport(txi.k5.tsv, colData = F3Popcoldata, design = ~ pop)
keep_R_F3Pop <- rowSums(counts(dds_R_F3Pop) >= 10) >= 3
dds_R_F3Popk <- dds_R_F3Pop[keep_R_F3Pop,]
dds_R_F3Popk <- DESeq(dds_R_F3Popk)
R_F3Pop_res <- lfcShrink(dds=dds_R_F3Popk, coef=2, type="apeglm")
summary(R_F3Pop_res)
R_F3Pop_res_sig<- subset(R_F3Pop_res, padj < 0.05)
R_F3Pop_C<-cbind(rownames(R_F3Pop_res),data.frame(R_F3Pop_res, row.names = NULL))

##### Population comparison F6
F6Popcoldata<-coldata[c(7:9,19:21),]
dds_R_F6Pop <- DESeqDataSetFromTximport(txi.k6.tsv, colData = F6Popcoldata, design = ~ pop)
keep_R_F6Pop <- rowSums(counts(dds_R_F6Pop) >= 10) >= 3
dds_R_F6Popk <- dds_R_F6Pop[keep_R_F6Pop,]
dds_R_F6Popk <- DESeq(dds_R_F6Popk)
R_F6Pop_res <- lfcShrink(dds=dds_R_F6Popk, coef=2, type="apeglm")
summary(R_F6Pop_res)
R_F6Pop_res_sig<- subset(R_F6Pop_res, padj < 0.05)
R_F6Pop_C<-cbind(rownames(R_F6Pop_res),data.frame(R_F6Pop_res, row.names = NULL))

#### Generation Effect Akaa
GenAcoldata<-coldata[c(1:3,7:9),]
dds_R_GenA <- DESeqDataSetFromTximport(txi.k7.tsv, colData = GenAcoldata, design = ~ gen)
keep_R_GenA <- rowSums(counts(dds_R_GenA) >= 10) >= 3
dds_R_GenAk <- dds_R_GenA[keep_R_GenA,]
dds_R_GenAk <- DESeq(dds_R_GenAk)
R_GenA_res <- lfcShrink(dds=dds_R_GenAk, coef=2, type="apeglm")
summary(R_GenA_res)
R_GenA_res_sig<- subset(R_GenA_res, padj < 0.05)
R_AkaaGen_C<-cbind(rownames(R_GenA_res_sig),data.frame(R_GenA_res_sig, row.names = NULL))

#### Generation Effect Manz
GenMcoldata<-coldata[c(13:15,19:21),]
dds_R_GenM <- DESeqDataSetFromTximport(txi.k8.tsv, colData = GenMcoldata, design = ~ gen)
keep_R_GenM <- rowSums(counts(dds_R_GenM) >= 10) >= 3
dds_R_GenMk <- dds_R_GenM[keep_R_GenM,]
dds_R_GenMk <- DESeq(dds_R_GenMk)
dds_R_GenMk_shr <- lfcShrink(dds=dds_R_GenMk, coef=2, type="apeglm")
R_GenM_res <- lfcShrink(dds=dds_R_GenMk, coef=2, type="apeglm")
summary(R_GenM_res)
R_GenM_res_sig<- subset(R_GenM_res, padj < 0.05)
R_ManzGen_C<-cbind(rownames(R_GenM_res_sig),data.frame(R_GenM_res_sig, row.names = NULL))

################################################################
# 1.3 Generate lists of sets of DEGs 

# List of genes from each DESeq2 comparison
rm(LIST_rna1.ids)
LIST_rna1.ids <- list()
LIST_rna1.ids[['F3_Akaa_Treat']] <- as.character(unique(rownames(R_F3A_res_sig)))
LIST_rna1.ids[['F3_Manz_Treat']] <- as.character(unique(rownames(R_F3M_res_sig)))
LIST_rna1.ids[['F6_Akaa_Treat']] <- as.character(unique(rownames(R_F6A_res_sig)))
LIST_rna1.ids[['F6_Manz_Treat']] <- as.character(unique(rownames(R_F6M_res_sig)))
LIST_rna1.ids[['F3_Population']] <- as.character(unique(rownames(R_F3Pop_res_sig)))
LIST_rna1.ids[['F6_Population']] <- as.character(unique(rownames(R_F6Pop_res_sig)))
LIST_rna1.ids[['Generation_Ak']] <- as.character(unique(rownames(R_GenA_res_sig)))
LIST_rna1.ids[['Generation_Ma']] <- as.character(unique(rownames(R_GenM_res_sig)))

# Lists of genes overlapping (or not) between List 1 groups
rm(LIST_rna2.ids)
LIST_rna2.ids <- list()
LIST_rna2.ids[['F3_Treat_overlap']] <- intersect(LIST_rna1.ids[[1]], LIST_rna1.ids[[2]])
LIST_rna2.ids[['F3_Trt_Akaa_uniq']] <- setdiff(LIST_rna1.ids[[1]], LIST_rna1.ids[[2]])
LIST_rna2.ids[['F3_Trt_Manz_uniq']] <- setdiff(LIST_rna1.ids[[2]], LIST_rna1.ids[[1]])
LIST_rna2.ids[['Pop_difs_overlap']] <- intersect(LIST_rna1.ids[[5]], LIST_rna1.ids[[6]])
LIST_rna2.ids[['Pop_difs_F3_uniq']] <- setdiff(LIST_rna1.ids[[5]], LIST_rna1.ids[[6]])
LIST_rna2.ids[['Pop_difs_F6_uniq']] <- setdiff(LIST_rna1.ids[[6]], LIST_rna1.ids[[5]])
LIST_rna2.ids[['F6_Treat_overlap']] <- intersect(LIST_rna1.ids[[3]], LIST_rna1.ids[[4]])
LIST_rna2.ids[['F6_Trt_Akaa_uniq']] <- setdiff(LIST_rna1.ids[[3]], LIST_rna1.ids[[4]])
LIST_rna2.ids[['F6_Trt_Manz_uniq']] <- setdiff(LIST_rna1.ids[[4]], LIST_rna1.ids[[3]])
LIST_rna2.ids[["F3A_F6A_ovlp"]]<- intersect(LIST_rna1.ids[[1]], LIST_rna1.ids[[3]])
LIST_rna2.ids[["F3M_F6M_ovlp"]]<- intersect(LIST_rna1.ids[[2]], LIST_rna1.ids[[4]])
LIST_rna2.ids[["F3A_F6A_ovlp_NoGenA"]]<- setdiff(LIST_rna2.ids[[10]], LIST_rna1.ids[[7]])
LIST_rna2.ids[["F3M_F6M_ovlp_NoGenM"]]<- setdiff(LIST_rna2.ids[[11]], LIST_rna1.ids[[8]])
LIST_rna2.ids[["F6A_uniq_NoGenA"]]<- setdiff(LIST_rna1.ids[[3]], c(LIST_rna1.ids[[1]],LIST_rna1.ids[[7]]))
LIST_rna2.ids[["F6M_uniq_NoGenM"]]<- setdiff(LIST_rna1.ids[[4]], c(LIST_rna1.ids[[2]],LIST_rna1.ids[[8]]))
LIST_rna2.ids[["F3A_uniq_NoGenA"]]<- setdiff(LIST_rna1.ids[[1]], c(LIST_rna1.ids[[3]],LIST_rna1.ids[[7]]))
LIST_rna2.ids[["F3M_uniq_NoGenM"]]<- setdiff(LIST_rna1.ids[[2]], c(LIST_rna1.ids[[4]],LIST_rna1.ids[[8]]))

# Lists of genes overlapping (or not) accounting for generational differences between controls, except
# the final 2 (which are repeated from LIST_rna1.ids) included to help with a eulerr plot
rm(LIST_rna3.ids)
LIST_rna3.ids <- list()
LIST_rna3.ids[['F3_Trt_Akaa_noGenA']] <- setdiff(LIST_rna1.ids[[1]], LIST_rna1.ids[[7]])
LIST_rna3.ids[['F3_Trt_Manz_noGenM']] <- setdiff(LIST_rna1.ids[[2]], LIST_rna1.ids[[8]])
LIST_rna3.ids[['F6_Trt_Akaa_noGenA']] <- setdiff(LIST_rna1.ids[[3]], LIST_rna1.ids[[7]])
LIST_rna3.ids[['F6_Trt_Manz_noGenM']] <- setdiff(LIST_rna1.ids[[4]], LIST_rna1.ids[[8]])
LIST_rna3.ids[['F3_Trt_overlap_noGenA']] <- intersect(LIST_rna3.ids[[1]], LIST_rna3.ids[[2]])
LIST_rna3.ids[['F3_Trt_Akaa_uniq_noGenM']] <- setdiff(LIST_rna3.ids[[1]], LIST_rna3.ids[[2]])
LIST_rna3.ids[['F3_Trt_Manz_uniq_noGenA']] <- setdiff(LIST_rna3.ids[[2]], LIST_rna3.ids[[1]])
LIST_rna3.ids[['F6_Trt_overlap_noGenA']] <- intersect(LIST_rna3.ids[[3]], LIST_rna3.ids[[4]])
LIST_rna3.ids[['F6_Trt_Akaa_uniq_noGenM']] <- setdiff(LIST_rna3.ids[[3]], LIST_rna3.ids[[4]])
LIST_rna3.ids[['F6_Trt_Manz_uniq_noGenA']] <- setdiff(LIST_rna3.ids[[4]], LIST_rna3.ids[[3]])
LIST_rna3.ids[['F3_Akaa_Treat']] <- as.character(unique(rownames(R_F3A_res_sig)))
LIST_rna3.ids[['F3_Manz_Treat']] <- as.character(unique(rownames(R_F3M_res_sig)))

# Optional: save lists to an r object
save(LIST_rna1.ids,LIST_rna2.ids,LIST_rna3.ids,file ="interchromate_downstream/deseq2_rnaseq_lists_250203.Rdata")

################################################################
################################################################
# 2. DEseq2 for ATAC 
################################################################
################################################################

# Load required packages
library(edgeR)
library(limma)

################################################################
# 2.1 Import results from intepareto as input

# Sample data is the same file used for RNA-seq DESEq2 analysis
#coldata<-read.table("interchromate_meta/Sample_info.txt",header=TRUE,row.names=1)

# Define Sample data
coldata$pop <- factor(coldata$pop)
coldata$gen <- factor(coldata$gen)
coldata$trt <- factor(coldata$trt)
coldata$pop_gen <- factor(paste(coldata$pop, coldata$gen, sep = "_"))
coldata$pop_trt <- factor(paste(coldata$pop, coldata$trt, sep = "_"))
coldata$gen_trt <- factor(paste(coldata$gen, coldata$trt, sep = "_"))
coldata$pop_gen_trt <- factor(paste(coldata$pop, coldata$gen, coldata$trt, sep = "_"))

# Read in the matched data output from intepareto for each comparison
# NB: when runinteParteo.R was run for each comparison, matchedData R objects were given the same generic name.
# Be careful when importing data: check the default column names to make sure they are as expected, 
# and make sure old and new match when reassigning column names

# 1. Control v Heatshock F3 Akaa
load("interchromate_intepareto/interchromate_AkaaF3_Ctrl_v_Heat1000.Rdata")
rownames(matchedData) <- matchedData[,1]
F3Acts<-matchedData[,c(11:13,8:10)]
colnames(F3Acts) <- c("Ak_F3_ctl_1", "Ak_F3_ctl_2", "Ak_F3_ctl_3", "Ak_F3_trt_1","Ak_F3_trt_2", "Ak_F3_trt_3")

# 2. Control v Heatshock F3 Manz
load("interchromate_intepareto/interchromate_ManzF3_Ctrl_v_Heat1000.Rdata")
rownames(matchedData) <- matchedData[,1]
F3Mcts<-matchedData[,c(11:13,8:10)]
colnames(F3Mcts) <- c("Ma_F3_ctl_1", "Ma_F3_ctl_2", "Ma_F3_ctl_3", "Ma_F3_trt_1","Ma_F3_trt_2", "Ma_F3_trt_3")

# 3. Control v Heatshock F6 Akaa
load("interchromate_intepareto/interchromate_AkaaF6_Ctrl_v_Heat1000.Rdata")
rownames(matchedData) <- matchedData[,1]
F6Acts<-matchedData[,c(11:13,8:10)]
colnames(F6Acts) <- c("Ak_F6_ctl_1", "Ak_F6_ctl_2", "Ak_F6_ctl_3", "Ak_F6_trt_1","Ak_F6_trt_2", "Ak_F6_trt_3")

# 4. Control v Heatshock F6 Manz
load("interchromate_intepareto/interchromate_ManzF6_Ctrl_v_Heat1000.Rdata")
rownames(matchedData) <- matchedData[,1]
F6Mcts<-matchedData[,c(11:13,8:10)]
colnames(F6Mcts) <- c("Ma_F6_ctl_1", "Ma_F6_ctl_2", "Ma_F6_ctl_3", "Ma_F6_trt_1","Ma_F6_trt_2", "Ma_F6_trt_3")

# 5. F3 Population Controls
load("interchromate_intepareto/interchromate_AkaaF3_Ctrl_v_ManzF3_Ctrl1000.Rdata")
rownames(matchedData) <- matchedData[,1]
F3Popcts<-matchedData[,c(11:13,8:10)]
colnames(F3Popcts) <- c("Ak_F3_ctl_1", "Ak_F3_ctl_2", "Ak_F3_ctl_3", "Ma_F3_ctl_1", "Ma_F3_ctl_2", "Ma_F3_ctl_3")

# 6. F6 Population Controls
load("interchromate_intepareto/interchromate_AkaaF6_Ctrl_v_ManzF6_Ctrl1000.Rdata")
rownames(matchedData) <- matchedData[,1]
F6Popcts<-matchedData[,c(11:13,8:10)]
colnames(F6Popcts) <- c("Ak_F6_ctl_1", "Ak_F6_ctl_2", "Ak_F6_ctl_3", "Ma_F6_ctl_1", "Ma_F6_ctl_2", "Ma_F6_ctl_3")

# 7. Generation Control Akaa
load("interchromate_intepareto/interchromate_AkaaF3_Ctrl_v_F6_Ctrl1000.Rdata")
rownames(matchedData) <- matchedData[,1]
GenActs<-matchedData[,c(11:13,8:10)]
colnames(GenActs) <- c("Ak_F3_ctl_1", "Ak_F3_ctl_2", "Ak_F3_ctl_3", "Ak_F6_ctl_1", "Ak_F6_ctl_2", "Ak_F6_ctl_3")

# 8. Generation Control Manz
load("interchromate_intepareto/interchromate_ManzF3_Ctrl_v_F6_Ctrl1000.Rdata")
rownames(matchedData) <- matchedData[,1]
GenMcts<-matchedData[,c(11:13,8:10)]
colnames(GenMcts) <- c("Ma_F3_ctl_1", "Ma_F3_ctl_2", "Ma_F3_ctl_3", "Ma_F6_ctl_1", "Ma_F6_ctl_2", "Ma_F6_ctl_3")

################################################################
# 2.2 Run DEseq2 and generate table of results for genes 
#     associated with sig DARs

##### Akaa F3
# Specifying which sample data
F3Acoldata<-coldata[c(1:6),]
# import the data
dds_A_F3A <- DESeqDataSetFromMatrix(countData = F3Acts, colData = F3Acoldata, design = ~ trt)
# exclude genes with low counts
keep_A_F3A <- rowSums(counts(dds_A_F3A) >= 0) >= 3
dds_A_F3Ak <- dds_A_F3A[keep_A_F3A,]
# Run DESeq2
dds_A_F3Ak <- DESeq(dds_A_F3Ak)
# Apply lfcshrink transformation
A_F3A_res <- lfcShrink(dds=dds_A_F3Ak, coef=2, type="apeglm")
# summarise the DESeq2 object
summary(A_F3A_res)
# Select on;y those genes with a significant difference
A_F3A_res_sig<- subset(A_F3A_res, padj < 0.05)
A_F3A_T<-cbind(rownames(A_F3A_res_sig),data.frame(A_F3A_res_sig, row.names = NULL))
# Get gene symbols 
A_F3A_res$name = mapIds(org.Dm.eg.db,keys=rownames(A_F3A_res), column="SYMBOL", keytype="FLYBASE",multiVals="first")
# Enhanced volcano plot
EnhancedVolcano(A_F3A_res, lab = A_F3A_res$name, x = 'log2FoldChange', y = 'pvalue', title = "Akaa F3", subtitle = "DA: HS vs Ctrl", pCutoff = 10e-4, FCcutoff = 0.585,legendPosition = "none")

##### Manz F3
F3Mcoldata<-coldata[c(13:18),]
dds_A_F3M <- DESeqDataSetFromMatrix(countData = F3Mcts, colData = F3Mcoldata, design = ~ trt)
keep_A_F3M <- rowSums(counts(dds_A_F3M) >= 0) >= 3
dds_A_F3Mk <- dds_A_F3M[keep_A_F3M,]
dds_A_F3Mk <- DESeq(dds_A_F3Mk)
A_F3M_res <- lfcShrink(dds=dds_A_F3Mk, coef=2, type="apeglm")
summary(A_F3M_res)
A_F3M_res_sig<- subset(A_F3M_res, padj < 0.05)
A_F3M_T<-cbind(rownames(A_F3M_res_sig),data.frame(A_F3M_res_sig, row.names = NULL))
A_F3M_res$name = mapIds(org.Dm.eg.db,keys=rownames(A_F3M_res), column="SYMBOL", keytype="FLYBASE",multiVals="first")
EnhancedVolcano(A_F3M_res, lab = A_F3M_res$name, x = 'log2FoldChange', y = 'pvalue', title = "Manz F3", subtitle = "DA: HS vs Ctrl", pCutoff = 10e-4, FCcutoff = 0.585,legendPosition = "none")

##### Akaa F6
F6Acoldata<-coldata[c(7:12),]
dds_A_F6A <- DESeqDataSetFromMatrix(countData = F6Acts, colData = F6Acoldata, design = ~ trt)
keep_A_F6A <- rowSums(counts(dds_A_F6A) >= 0) >= 3
dds_A_F6Ak <- dds_A_F6A[keep_A_F6A,]
dds_A_F6Ak <- DESeq(dds_A_F6Ak)
A_F6A_res <- lfcShrink(dds=dds_A_F6Ak, coef=2, type="apeglm")
summary(A_F6A_res)
A_F6A_res_sig<- subset(A_F6A_res, padj < 0.05)
A_F6A_T<-cbind(rownames(A_F6A_res_sig),data.frame(A_F6A_res_sig, row.names = NULL))
A_F6A_res$name = mapIds(org.Dm.eg.db,keys=rownames(A_F6A_res), column="SYMBOL", keytype="FLYBASE",multiVals="first")
EnhancedVolcano(A_F6A_res, lab = A_F6A_res$name, x = 'log2FoldChange', y = 'pvalue', title = "Manz F3", subtitle = "DA: HS vs Ctrl", pCutoff = 10e-4, FCcutoff = 0.3219,legendPosition = "none")

##### Manz F6
F6Mcoldata<-coldata[c(19:24),]
dds_A_F6M <- DESeqDataSetFromMatrix(countData = F6Mcts, colData = F6Mcoldata,  design = ~ trt)
keep_A_F6M <- rowSums(counts(dds_A_F6M) >= 0) >= 3
dds_A_F6Mk <- dds_A_F6M[keep_A_F6M,]
dds_A_F6Mk <- DESeq(dds_A_F6Mk)
A_F6M_res <- lfcShrink(dds=dds_A_F6Mk, coef=2, type="apeglm")
summary(A_F6M_res)
A_F6M_res_sig<- subset(A_F6M_res, padj < 0.05)
A_F6M_T<-cbind(rownames(A_F6M_res_sig),data.frame(A_F6M_res_sig, row.names = NULL))
A_F6M_res$name = mapIds(org.Dm.eg.db,keys=rownames(A_F6M_res), column="SYMBOL", keytype="FLYBASE",multiVals="first")
EnhancedVolcano(A_F6M_res, lab = A_F6M_res$name, x = 'log2FoldChange', y = 'pvalue', title = "Manz F6", subtitle = "DA: HS vs Ctrl", pCutoff = 10e-4, FCcutoff = 0.3219,legendPosition = "none")

##### Population comparison F3
F3Popcoldata<-coldata[c(1:3,13:15),]
dds_A_F3Pop <- DESeqDataSetFromMatrix(countData = F3Popcts, colData = F3Popcoldata,  design = ~ pop)
keep_A_F3Pop <- rowSums(counts(dds_A_F3Pop) >= 10) >= 3
dds_A_F3Popk <- dds_A_F3Pop[keep_A_F3Pop,]
dds_A_F3Popk <- DESeq(dds_A_F3Popk)
A_F3Pop_res <- lfcShrink(dds=dds_A_F3Popk, coef=2, type="apeglm")
summary(A_F3Pop_res)
A_F3Pop_res_sig<- subset(A_F3Pop_res, padj < 0.05)
A_F3Pop_C<-cbind(rownames(A_F3Pop_res),data.frame(A_F3Pop_res, row.names = NULL))

##### Population comparison F6
F6Popcoldata<-coldata[c(7:9,19:21),]
dds_A_F6Pop <- DESeqDataSetFromMatrix(countData = F6Popcts, colData = F6Popcoldata, design = ~ pop)
keep_A_F6Pop <- rowSums(counts(dds_A_F6Pop) >= 10) >= 3
dds_A_F6Popk <- dds_A_F6Pop[keep_A_F6Pop,]
dds_A_F6Popk <- DESeq(dds_A_F6Popk)
A_F6Pop_res <- lfcShrink(dds=dds_A_F6Popk, coef=2, type="apeglm")
summary(A_F6Pop_res)
A_F6Pop_res_sig<- subset(A_F6Pop_res, padj < 0.05)
A_F6Pop_C<-cbind(rownames(A_F6Pop_res),data.frame(A_F6Pop_res, row.names = NULL))

### Generation Effect Akaa
GenAcoldata<-coldata[c(1:3,7:9),]
dds_A_GenA <- DESeqDataSetFromMatrix(countData = GenActs, colData = GenAcoldata, design = ~ gen)
keep_A_GenA <- rowSums(counts(dds_A_GenA) >= 10) >= 3
dds_A_GenAk <- dds_A_GenA[keep_A_GenA,]
dds_A_GenAk <- DESeq(dds_A_GenAk)
A_GenA_res <- lfcShrink(dds=dds_A_GenAk, coef=2, type="apeglm")
summary(A_GenA_res)
A_GenA_res_sig<- subset(A_GenA_res, padj < 0.05)
A_AkaaGen_C<-cbind(rownames(A_GenA_res_sig),data.frame(A_GenA_res_sig, row.names = NULL))

### Generation Effect Manz
GenMcoldata<-coldata[c(13:15,19:21),]
dds_A_GenM <- DESeqDataSetFromMatrix(countData = GenMcts, colData = GenMcoldata, design = ~ gen)
keep_A_GenM <- rowSums(counts(dds_A_GenM) >= 10) >= 3
dds_A_GenMk <- dds_A_GenM[keep_A_GenM,]
dds_A_GenMk <- DESeq(dds_A_GenMk)
A_GenM_res <- lfcShrink(dds=dds_A_GenMk, coef=2, type="apeglm")
summary(A_GenM_res)
A_GenM_res_sig<- subset(A_GenM_res, padj < 0.05)
A_ManzGen_C<-cbind(rownames(A_GenM_res_sig),data.frame(A_GenM_res_sig, row.names = NULL))

################################################################
# 2.3 Generate lists of sets of genes associated with DARs 

# List of genes from each DESeq2 comparison
rm(LIST_atac1.ids)
LIST_atac1.ids <- list()
LIST_atac1.ids[['F3_Akaa_Treat']] <- as.character(unique(rownames(A_F3A_res_sig)))
LIST_atac1.ids[['F3_Manz_Treat']] <- as.character(unique(rownames(A_F3M_res_sig)))
LIST_atac1.ids[['F6_Akaa_Treat']] <- as.character(unique(rownames(A_F6A_res_sig)))
LIST_atac1.ids[['F6_Manz_Treat']] <- as.character(unique(rownames(A_F6M_res_sig)))
LIST_atac1.ids[['F3_Population']] <- as.character(unique(rownames(A_F3Pop_res_sig)))
LIST_atac1.ids[['F6_Population']] <- as.character(unique(rownames(A_F6Pop_res_sig)))
LIST_atac1.ids[['Generation_Ak']] <- as.character(unique(rownames(A_GenA_res_sig)))
LIST_atac1.ids[['Generation_Ma']] <- as.character(unique(rownames(A_GenM_res_sig)))

# Lists of genes overlapping (or not) between List 1 groups
rm(LIST_atac2.ids)
LIST_atac2.ids <- list()
LIST_atac2.ids[['F3_Treat_overlap']] <- intersect(LIST_atac1.ids[[1]], LIST_atac1.ids[[2]])
LIST_atac2.ids[['F3_Trt_Akaa_uniq']] <- setdiff(LIST_atac1.ids[[1]], LIST_atac1.ids[[2]])
LIST_atac2.ids[['F3_Trt_Manz_uniq']] <- setdiff(LIST_atac1.ids[[2]], LIST_atac1.ids[[1]])
LIST_atac2.ids[['Pop_difs_overlap']] <- intersect(LIST_atac1.ids[[5]], LIST_atac1.ids[[6]])
LIST_atac2.ids[['Pop_difs_F3_uniq']] <- setdiff(LIST_atac1.ids[[5]], LIST_atac1.ids[[6]])
LIST_atac2.ids[['Pop_difs_F6_uniq']] <- setdiff(LIST_atac1.ids[[6]], LIST_atac1.ids[[5]])
LIST_atac2.ids[['F6_Treat_overlap']] <- intersect(LIST_atac1.ids[[3]], LIST_atac1.ids[[4]])
LIST_atac2.ids[['F6_Trt_Akaa_uniq']] <- setdiff(LIST_atac1.ids[[3]], LIST_atac1.ids[[4]])
LIST_atac2.ids[['F6_Trt_Manz_uniq']] <- setdiff(LIST_atac1.ids[[4]], LIST_atac1.ids[[3]])
LIST_atac2.ids[["F3A_F6A_ovlp"]]<- intersect(LIST_atac1.ids[[1]], LIST_atac1.ids[[3]])
LIST_atac2.ids[["F3M_F6M_ovlp"]]<- intersect(LIST_atac1.ids[[2]], LIST_atac1.ids[[4]])
LIST_atac2.ids[["F3A_F6A_ovlp_NoGenA"]]<- setdiff(LIST_atac2.ids[[10]], LIST_atac1.ids[[7]])
LIST_atac2.ids[["F3M_F6M_ovlp_NoGenM"]]<- setdiff(LIST_atac2.ids[[11]], LIST_atac1.ids[[8]])
LIST_atac2.ids[["F6A_uniq_NoGenA"]]<- setdiff(LIST_atac1.ids[[3]], c(LIST_atac1.ids[[1]],LIST_atac1.ids[[7]]))
LIST_atac2.ids[["F6M_uniq_NoGenM"]]<- setdiff(LIST_atac1.ids[[4]], c(LIST_atac1.ids[[2]],LIST_atac1.ids[[8]]))
LIST_atac2.ids[["F3A_uniq_NoGenA"]]<- setdiff(LIST_atac1.ids[[1]], c(LIST_atac1.ids[[3]],LIST_atac1.ids[[7]]))
LIST_atac2.ids[["F3M_uniq_NoGenM"]]<- setdiff(LIST_atac1.ids[[2]], c(LIST_atac1.ids[[4]],LIST_atac1.ids[[8]]))

# Lists of genes overlapping (or not) accounting for generational differences between controls, except
# the final 2 (which are repeated from LIST_rna1.ids) included to help with a eulerr plotrm(LIST_atac3.ids)
LIST_atac3.ids <- list()
LIST_atac3.ids[['F3_Trt_Akaa_noGenA']] <- setdiff(LIST_atac1.ids[[1]], LIST_atac1.ids[[7]])
LIST_atac3.ids[['F3_Trt_Manz_noGenM']] <- setdiff(LIST_atac1.ids[[2]], LIST_atac1.ids[[8]])
LIST_atac3.ids[['F6_Trt_Akaa_noGenA']] <- setdiff(LIST_atac1.ids[[3]], LIST_atac1.ids[[7]])
LIST_atac3.ids[['F6_Trt_Manz_noGenM']] <- setdiff(LIST_atac1.ids[[4]], LIST_atac1.ids[[8]])
LIST_atac3.ids[['F3_Trt_overlap_noGenA']] <- intersect(LIST_atac3.ids[[1]], LIST_atac3.ids[[2]])
LIST_atac3.ids[['F3_Trt_Akaa_uniq_noGenM']] <- setdiff(LIST_atac3.ids[[1]], LIST_atac3.ids[[2]])
LIST_atac3.ids[['F3_Trt_Manz_uniq_noGenA']] <- setdiff(LIST_atac3.ids[[2]], LIST_atac3.ids[[1]])
LIST_atac3.ids[['F6_Trt_overlap_noGenA']] <- intersect(LIST_atac3.ids[[3]], LIST_atac3.ids[[4]])
LIST_atac3.ids[['F6_Trt_Akaa_uniq_noGenM']] <- setdiff(LIST_atac3.ids[[3]], LIST_atac3.ids[[4]])
LIST_atac3.ids[['F6_Trt_Manz_uniq_noGenA']] <- setdiff(LIST_atac3.ids[[4]], LIST_atac3.ids[[3]])
LIST_atac3.ids[['F3_Akaa_Treat']] <- as.character(unique(rownames(A_F3A_res_sig)))
LIST_atac3.ids[['F3_Manz_Treat']] <- as.character(unique(rownames(A_F3M_res_sig)))

# Optional: save lists to an r object
save(LIST_atac1.ids,LIST_atac2.ids,LIST_atac3.ids,file ="interchromate_downstream/deseq2_atacseq_lists_250203.Rdata")

################################################################
################################################################
# 3. Combining results and creating complete data table   
################################################################
################################################################

# Load required packages
library(dplyr)

# Load lists of DE and DA genes if sections 1 and 2 have already run
#load("interchromate_downstream/deseq2_rnaseq_lists_250203.Rdata")
#load("interchromate_downstream/deseq2_atacseq_lists_250203.Rdata")

################################################################
# 3.1 Generate lists of sets of genes that are both DE 
#     and associated with DARs

# Overlap between DEGs and DARs
rm(LIST_rna_atac1.ids)
LIST_rna_atac1.ids <- list()
LIST_rna_atac1.ids[['F3_Akaa_ovlp']] <- intersect(LIST_atac1.ids[[1]], LIST_rna1.ids[[1]])
LIST_rna_atac1.ids[['F3_Manz_ovlp']] <- intersect(LIST_atac1.ids[[2]], LIST_rna1.ids[[2]])
LIST_rna_atac1.ids[['F3_ovlp_ovlp']] <- intersect(LIST_atac2.ids[[1]], LIST_rna2.ids[[1]])
LIST_rna_atac1.ids[['F3_AkUniq_ovlp']] <- intersect(LIST_atac2.ids[[2]], LIST_rna2.ids[[2]])
LIST_rna_atac1.ids[['F3_MaUniq_ovlp']] <- intersect(LIST_atac2.ids[[3]], LIST_rna2.ids[[3]])

# OVerlap or differences between populations in genes that are both DEGs and DARs
rm(LIST_rna_atac2.ids)
LIST_rna_atac2.ids <- list()
LIST_rna_atac2.ids[['F3_Akaa_rna']] <- LIST_rna1.ids[[1]]
LIST_rna_atac2.ids[['F3_Akaa_atac']] <- LIST_atac1.ids[[1]]
LIST_rna_atac2.ids[['F3_Manz_rna']] <-  LIST_rna1.ids[[2]]
LIST_rna_atac2.ids[['F3_Manz_atac']] <- LIST_atac1.ids[[2]]
LIST_rna_atac2.ids[['core']] <- intersect(LIST_rna_atac1.ids[[1]], LIST_rna_atac1.ids[[2]])
LIST_rna_atac2.ids[['Ak_reg']] <- setdiff(LIST_rna_atac1.ids[[1]], LIST_rna_atac1.ids[[2]])
LIST_rna_atac2.ids[['Ma_reg']] <- setdiff(LIST_rna_atac1.ids[[2]], LIST_rna_atac1.ids[[1]])


################################################################
# 3.2 Read in the integrated data output from intepareto for 
#     each heat shock comparison
#     NB: when runinteParteo.R was run for each comparison, 
#     Integ_zscores R objects were given the same generic name.
#     Be careful when importing data (make sure R.data 
#     objects are distinct and called correctly)

# Control v Heatshock F3 Akaa
load("interchromate_intepareto/interchromate_AkaaF3_Ctrl_v_Heat1000.Rdata")
F3A_integ<-Integ_zscores
# Control v Heatshock F3 Manz
load("interchromate_intepareto/interchromate_ManzF3_Ctrl_v_Heat1000.Rdata")
F3M_integ<-Integ_zscores
# Control v Heatshock F6 Akaa
load("interchromate_intepareto/interchromate_AkaaF6_Ctrl_v_Heat1000.Rdata")
F6A_integ<-Integ_zscores
# Control v Heatshock F6 Manz
load("interchromate_intepareto/interchromate_ManzF6_Ctrl_v_Heat1000.Rdata")
F6M_integ<-Integ_zscores

################################################################
# 3.3 Read in TE data for Pop1 (Akaa) and Pop2 (Manz), 
#     reference (ref) and non-ref (den)

gene_TEs_P1r<-read.table("interchromate_TEconsensus/genes_flybase_TE_F3_P1_strict_promoter_1000_ref.lst", header = F, stringsAsFactors = F)
gene_TEs_P1d<-read.table("interchromate_TEconsensus/genes_flybase_TE_F3_P1_strict_promoter_1000_den.lst", header = F, stringsAsFactors = F)
gene_TEs_P2r<-read.table("interchromate_TEconsensus/genes_flybase_TE_F3_P2_strict_promoter_1000_ref.lst", header = F, stringsAsFactors = F)
gene_TEs_P2d<-read.table("interchromate_TEconsensus/genes_flybase_TE_F3_P2_strict_promoter_1000_den.lst", header = F, stringsAsFactors = F)

################################################################
# 3.4 Assign genes new value based on proximity to TE 
#     (0 vs R vs D | 0 vs 1)
#     R for reference, D for non-reference (denovo)

F3A_integ$TE_pres <- "0"
F3A_integ[F3A_integ$ensembl_gene_id%in%gene_TEs_P1r$V1,]$TE_pres <- "R"
F3A_integ[F3A_integ$ensembl_gene_id%in%gene_TEs_P1d$V1,]$TE_pres <- "D"
F3A_integ$TE_pres<-factor(F3A_integ$TE_pres,levels =c("0","R","D"))

F3A_integ$TE_pres2 <- "0"
F3A_integ[F3A_integ$ensembl_gene_id%in%gene_TEs_P1r$V1,]$TE_pres2 <- "1"
F3A_integ[F3A_integ$ensembl_gene_id%in%gene_TEs_P1d$V1,]$TE_pres2 <- "1"

F3M_integ$TE_pres <- "0"
F3M_integ[F3M_integ$ensembl_gene_id%in%gene_TEs_P2r$V1,]$TE_pres <- "R"
F3M_integ[F3M_integ$ensembl_gene_id%in%gene_TEs_P2d$V1,]$TE_pres <- "D"
F3M_integ$TE_pres<-factor(F3M_integ$TE_pres,levels =c("0","R","D"))

F3M_integ$TE_pres2 <- "0"
F3M_integ[F3M_integ$ensembl_gene_id%in%gene_TEs_P2r$V1,]$TE_pres2 <- "1"
F3M_integ[F3M_integ$ensembl_gene_id%in%gene_TEs_P2d$V1,]$TE_pres2 <- "1"

F6A_integ$TE_pres <- "0"
F6A_integ[F6A_integ$ensembl_gene_id%in%gene_TEs_P1r$V1,]$TE_pres <- "R"
F6A_integ[F6A_integ$ensembl_gene_id%in%gene_TEs_P1d$V1,]$TE_pres <- "D"
F6A_integ$TE_pres<-factor(F6A_integ$TE_pres,levels =c("0","R","D"))

F6A_integ$TE_pres2 <- "0"
F6A_integ[F6A_integ$ensembl_gene_id%in%gene_TEs_P1r$V1,]$TE_pres2 <- "1"
F6A_integ[F6A_integ$ensembl_gene_id%in%gene_TEs_P1d$V1,]$TE_pres2 <- "1"

F6M_integ$TE_pres <- "0"
F6M_integ[F6M_integ$ensembl_gene_id%in%gene_TEs_P2r$V1,]$TE_pres <- "R"
F6M_integ[F6M_integ$ensembl_gene_id%in%gene_TEs_P2d$V1,]$TE_pres <- "D"
F6M_integ$TE_pres<-factor(F6M_integ$TE_pres,levels =c("0","R","D"))

F6M_integ$TE_pres2 <- "0"
F6M_integ[F6M_integ$ensembl_gene_id%in%gene_TEs_P2r$V1,]$TE_pres2 <- "1"
F6M_integ[F6M_integ$ensembl_gene_id%in%gene_TEs_P2d$V1,]$TE_pres2 <- "1"

################################################################
# 3.5  Assign information about DE status (not DE, shared by 
#      both pops, unique) bases on LIST_rna objects

# F3 Akaa - DE in F3?
F3A_integ$F3A_DE<-"0"
F3A_integ[F3A_integ$ensembl_gene_id %in% LIST_rna2.ids[['F3_Treat_overlap']], ]$F3A_DE<-"Shared"
F3A_integ[F3A_integ$ensembl_gene_id %in% LIST_rna2.ids[['F3_Trt_Akaa_uniq']], ]$F3A_DE<-"Unique"
F3A_integ[!F3A_integ$ensembl_gene_id %in% LIST_rna1.ids[['F3_Akaa_Treat']], ]$F3A_DE<-"Not_DE"

F3A_DE_df <- F3A_integ %>%
  group_by(F3A_DE) %>%
  dplyr::summarise(count = n(),mean = mean(z.atac, na.rm = TRUE), sd = sd(z.atac, na.rm = TRUE)) %>% as.data.frame()
F3A_DE_df

# F3 Manz - DE in F3?
F3M_integ$F3M_DE<-"0"
F3M_integ[F3M_integ$ensembl_gene_id %in% LIST_rna2.ids[['F3_Treat_overlap']], ]$F3M_DE<-"Shared"
F3M_integ[F3M_integ$ensembl_gene_id %in% LIST_rna2.ids[['F3_Trt_Manz_uniq']], ]$F3M_DE<-"Unique"
F3M_integ[!F3M_integ$ensembl_gene_id %in% LIST_rna1.ids[['F3_Manz_Treat']], ]$F3M_DE<-"Not_DE"

F3M_DE_df <- F3M_integ %>%
  group_by(F3M_DE) %>%
  dplyr::summarise(count = n(),mean = mean(z.atac, na.rm = TRUE), sd = sd(z.atac, na.rm = TRUE)) %>% as.data.frame()
F3M_DE_df

# F6 Akaa - DE in F3 (prev gen)?
F6A_integ$F3A_DE<-"0"
F6A_integ[F6A_integ$ensembl_gene_id %in% LIST_rna2.ids[['F3_Treat_overlap']], ]$F3A_DE<-"Shared"
F6A_integ[F6A_integ$ensembl_gene_id %in% LIST_rna2.ids[['F3_Trt_Akaa_uniq']], ]$F3A_DE<-"Unique"
F6A_integ[!F6A_integ$ensembl_gene_id %in% LIST_rna1.ids[['F3_Akaa_Treat']], ]$F3A_DE<-"Not_DE"

F6A_DEinF3_df <- F6A_integ %>%
  group_by(F3A_DE) %>%
  dplyr::summarise(count = n(),mean = mean(z.atac, na.rm = TRUE), sd = sd(z.atac, na.rm = TRUE)) %>% as.data.frame()
F6A_DEinF3_df

# F6 Manz - DE in F3 (prev gen)?
F6M_integ$F3M_DE<-"0"
F6M_integ[F6M_integ$ensembl_gene_id %in% LIST_rna2.ids[['F3_Treat_overlap']], ]$F3M_DE<-"Shared"
F6M_integ[F6M_integ$ensembl_gene_id %in% LIST_rna2.ids[['F3_Trt_Manz_uniq']], ]$F3M_DE<-"Unique"
F6M_integ[!F6M_integ$ensembl_gene_id %in% LIST_rna1.ids[['F3_Manz_Treat']], ]$F3M_DE<-"Not_DE"

F6M_DEinF3_df <- F6M_integ %>%
  group_by(F3M_DE) %>%
  dplyr::summarise(count = n(),mean = mean(z.atac, na.rm = TRUE), sd = sd(z.atac, na.rm = TRUE)) %>% as.data.frame()
F6M_DEinF3_df

# F6 Akaa - DE in F6 (same gen)?
F6A_integ$F6A_DE<-"0"
F6A_integ[F6A_integ$ensembl_gene_id %in% LIST_rna2.ids[['F6_Treat_overlap']], ]$F6A_DE<-"Shared"
F6A_integ[F6A_integ$ensembl_gene_id %in% LIST_rna2.ids[['F6_Trt_Akaa_uniq']], ]$F6A_DE<-"Unique"
F6A_integ[!F6A_integ$ensembl_gene_id %in% LIST_rna1.ids[['F6_Akaa_Treat']], ]$F6A_DE<-"Not_DE"

F6A_DE_df <- F6A_integ %>%
  group_by(F6A_DE) %>%
  dplyr::summarise(count = n(),mean = mean(z.atac, na.rm = TRUE), sd = sd(z.atac, na.rm = TRUE)) %>% as.data.frame()
F6A_DE_df

# F6 Manz - DE in F6 (same gen)?
F6M_integ$F6M_DE<-"0"
F6M_integ[F6M_integ$ensembl_gene_id %in% LIST_rna2.ids[['F6_Treat_overlap']], ]$F6M_DE<-"Shared"
F6M_integ[F6M_integ$ensembl_gene_id %in% LIST_rna2.ids[['F6_Trt_Manz_uniq']], ]$F6M_DE<-"Unique"
F6M_integ[!F6M_integ$ensembl_gene_id %in% LIST_rna1.ids[['F6_Manz_Treat']], ]$F6M_DE<-"Not_DE"

F6M_DE_df <- F6M_integ %>%
  group_by(F6M_DE) %>%
  dplyr::summarise(count = n(),mean = mean(z.atac, na.rm = TRUE), sd = sd(z.atac, na.rm = TRUE)) %>% as.data.frame()
F6M_DE_df

################################################################
# 3.6  Assign information about DA status (not DA, shared by 
#      both pops, unique) bases on LIST_atac objects

# F3 Akaa - DA in F3?
F3A_integ$F3A_DA<-"0"
F3A_integ[F3A_integ$ensembl_gene_id %in% LIST_atac2.ids[['F3_Treat_overlap']], ]$F3A_DA<-"Shared"
F3A_integ[F3A_integ$ensembl_gene_id %in% LIST_atac2.ids[['F3_Trt_Akaa_uniq']], ]$F3A_DA<-"Unique"
F3A_integ[!F3A_integ$ensembl_gene_id %in% LIST_atac1.ids[['F3_Akaa_Treat']], ]$F3A_DA<-"Not_DA"

F3A_DA_df <- F3A_integ %>%
  group_by(F3A_DA) %>%
  dplyr::summarise(count = n(),mean = mean(z.atac, na.rm = TRUE), sd = sd(z.atac, na.rm = TRUE)) %>% as.data.frame()
F3A_DA_df

# F3 Manz - DA in F3?
F3M_integ$F3M_DA<-"0"
F3M_integ[F3M_integ$ensembl_gene_id %in% LIST_atac2.ids[['F3_Treat_overlap']], ]$F3M_DA<-"Shared"
F3M_integ[F3M_integ$ensembl_gene_id %in% LIST_atac2.ids[['F3_Trt_Manz_uniq']], ]$F3M_DA<-"Unique"
F3M_integ[!F3M_integ$ensembl_gene_id %in% LIST_atac1.ids[['F3_Manz_Treat']], ]$F3M_DA<-"Not_DA"

F3M_DA_df <- F3M_integ %>%
  group_by(F3M_DA) %>%
  dplyr::summarise(count = n(),mean = mean(z.atac, na.rm = TRUE), sd = sd(z.atac, na.rm = TRUE)) %>% as.data.frame()
F3M_DA_df

# F6 Akaa - DA in F3 (prev gen)?
F6A_integ$F3A_DA<-"0"
F6A_integ[F6A_integ$ensembl_gene_id %in% LIST_atac2.ids[['F3_Treat_overlap']], ]$F3A_DA<-"Shared"
F6A_integ[F6A_integ$ensembl_gene_id %in% LIST_atac2.ids[['F3_Trt_Akaa_uniq']], ]$F3A_DA<-"Unique"
F6A_integ[!F6A_integ$ensembl_gene_id %in% LIST_atac1.ids[['F3_Akaa_Treat']], ]$F3A_DA<-"Not_DA"

F6A_DA_df <- F6A_integ %>%
  group_by(F3A_DA) %>%
  dplyr::summarise(count = n(),mean = mean(z.atac, na.rm = TRUE), sd = sd(z.atac, na.rm = TRUE)) %>% as.data.frame()
F6A_DA_df

# F6 Manz - DA in F3 (prev gen)?
F6M_integ$F3M_DA<-"0"
F6M_integ[F6M_integ$ensembl_gene_id %in% LIST_atac2.ids[['F3_Treat_overlap']], ]$F3M_DA<-"Shared"
F6M_integ[F6M_integ$ensembl_gene_id %in% LIST_atac2.ids[['F3_Trt_Manz_uniq']], ]$F3M_DA<-"Unique"
F6M_integ[!F6M_integ$ensembl_gene_id %in% LIST_atac1.ids[['F3_Manz_Treat']], ]$F3M_DA<-"Not_DA"

F6M_DA_df <- F6M_integ %>%
  group_by(F3M_DA) %>%
  dplyr::summarise(count = n(),mean = mean(z.atac, na.rm = TRUE), sd = sd(z.atac, na.rm = TRUE)) %>% as.data.frame()
F6M_DA_df

# Little DA in F6 (3 genes in Manz) so don't add this as a column to the data table
# But here are the genes that were DA in the F6:
F6M_integ[F6M_integ$ensembl_gene_id %in% LIST_atac1.ids[['F6_Manz_Treat']], ]
# FBgn0260962 (piccolo), FBgn0038115 (methanethiol oxidase), FBgn0040071 (tara)
# reduced to only 2 genes when accounting for generational differences
F6M_integ[F6M_integ$ensembl_gene_id %in% LIST_atac3.ids[['F6_Trt_Manz_noGenM']],]
# FBgn0260962 (piccolo), FBgn0038115 (methanethiol oxidase)

################################################################
# 3.7 Combine F3 and F6 data sets

Akaa_intergen <- merge(F3A_integ[,c(1:4,7:10)],F6A_integ[,c(1:4,10)],by = "ensembl_gene_id")
Manz_intergen <- merge(F3M_integ[,c(1:4,7:10)],F6M_integ[,c(1:4,10)],by = "ensembl_gene_id")

################################################################
# 3.8 Add regulation information (when genes are both DE and DA)
#     and other information related to directional change

# State which genes show DE and DA in F3 (candidates for regulated or misregulated)
Akaa_intergen$DEGDAG<-"0"
Akaa_intergen[Akaa_intergen$ensembl_gene_id %in% LIST_rna_atac2.ids[['core']], ]$DEGDAG<-"Shared"
Akaa_intergen[Akaa_intergen$ensembl_gene_id %in% LIST_rna_atac2.ids[['Ak_reg']], ]$DEGDAG<-"Unique"

Manz_intergen$DEGDAG<-"0"
Manz_intergen[Manz_intergen$ensembl_gene_id %in% LIST_rna_atac2.ids[['core']], ]$DEGDAG<-"Shared"
Manz_intergen[Manz_intergen$ensembl_gene_id %in% LIST_rna_atac2.ids[['Ma_reg']], ]$DEGDAG<-"Unique"

# Qualify regulation vs misregulation and the direction
Akaa_intergen$Zdir<-"0"
Akaa_intergen[Akaa_intergen$RNAseq.log2FoldChange.x > 0 & Akaa_intergen$atac.log2FoldChange.x > 0, ]$Zdir<-"RegUp"
Akaa_intergen[Akaa_intergen$RNAseq.log2FoldChange.x < 0 & Akaa_intergen$atac.log2FoldChange.x < 0, ]$Zdir<-"RegDn"
Akaa_intergen[Akaa_intergen$RNAseq.log2FoldChange.x > 0 & Akaa_intergen$atac.log2FoldChange.x < 0, ]$Zdir<-"MisRegExUp"
Akaa_intergen[Akaa_intergen$RNAseq.log2FoldChange.x < 0 & Akaa_intergen$atac.log2FoldChange.x > 0, ]$Zdir<-"MisRegExDn"

Manz_intergen$Zdir<-"0"
Manz_intergen[Manz_intergen$RNAseq.log2FoldChange.x > 0 & Manz_intergen$atac.log2FoldChange.x > 0, ]$Zdir<-"RegUp"
Manz_intergen[Manz_intergen$RNAseq.log2FoldChange.x < 0 & Manz_intergen$atac.log2FoldChange.x < 0, ]$Zdir<-"RegDn"
Manz_intergen[Manz_intergen$RNAseq.log2FoldChange.x > 0 & Manz_intergen$atac.log2FoldChange.x < 0, ]$Zdir<-"MisRegExUp"
Manz_intergen[Manz_intergen$RNAseq.log2FoldChange.x < 0 & Manz_intergen$atac.log2FoldChange.x > 0, ]$Zdir<-"MisRegExDn"

# Qualify the direction of RNA-seq changes 
Akaa_intergen$Rdir<-"0"
Akaa_intergen[Akaa_intergen$RNAseq.log2FoldChange.x > 0, ]$Rdir<-"RNAUp"
Akaa_intergen[Akaa_intergen$RNAseq.log2FoldChange.x < 0, ]$Rdir<-"RNADn"

Manz_intergen$Rdir<-"0"
Manz_intergen[Manz_intergen$RNAseq.log2FoldChange.x > 0, ]$Rdir<-"RNAUp"
Manz_intergen[Manz_intergen$RNAseq.log2FoldChange.x < 0, ]$Rdir<-"RNADn"

# Qualifying the direction of atac-seq changes 
Akaa_intergen$Adir<-"0"
Akaa_intergen[Akaa_intergen$atac.log2FoldChange.x > 0, ]$Adir<-"ATACUp"
Akaa_intergen[Akaa_intergen$atac.log2FoldChange.x < 0, ]$Adir<-"ATACDn"

Manz_intergen$Adir<-"0"
Manz_intergen[Manz_intergen$atac.log2FoldChange.x > 0, ]$Adir<-"ATACUp"
Manz_intergen[Manz_intergen$atac.log2FoldChange.x < 0, ]$Adir<-"ATACDn"

# Make sure this directional information is coded as factors
Akaa_intergen$Zdir <- factor(Akaa_intergen$Zdir, levels = c("MisRegExDn", "MisRegExUp","RegDn", "RegUp"))
Manz_intergen$Zdir <- factor(Manz_intergen$Zdir, levels = c("MisRegExDn", "MisRegExUp","RegDn", "RegUp"))

Akaa_intergen$Rdir <- factor(Akaa_intergen$Rdir, levels = c("RNADn","RNAUp"))
Manz_intergen$Rdir <- factor(Manz_intergen$Rdir, levels = c("RNADn","RNAUp"))

Akaa_intergen$Adir <- factor(Akaa_intergen$Adir, levels = c("ATACDn","ATACUp"))
Manz_intergen$Adir <- factor(Manz_intergen$Adir, levels = c("ATACDn","ATACUp"))

# Summary table with all=TRUE to keep all results

summary_table_Fbgn <- merge(Akaa_intergen[,c(1,5,2,3,4,14,7,8,9,10,11,12)], Manz_intergen[,c(1,5,2,3,4,14,7,8,9,10,11,12)], by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all = TRUE)

## Save this final object combining all information, and also make an R object
write.table(summary_table_Fbgn,file = "interchromate_downstream/summary_table_Fbgn.txt", sep = "\t")
save(Akaa_intergen,Manz_intergen,F3A_integ,F6A_integ,F3M_integ,F6M_integ,summary_table_Fbgn, file = "interchromate_downstream/interchromate_integrated_results_250203.Rdata")

# Can load it later if already complete
#load("interchromate_downstream/interchromate_integrated_results_250203.Rdata")

################################################################
# 3.9 Using the Akaa TEs and Manz TEs data (read in at start of section 3) we
# will create tables of genes that are close to polymorphic TEs and which have
# unique patterns of expression or accessbility

## Going to load an alternative version of the TE data for these output tables that have more 
## information about insertions and their associated genes
Akaa_TEs <- read.table("interchromate_TEconsensus/TE_P1_all_insertions_and_genes.txt", header = FALSE)
Manz_TEs <- read.table("interchromate_TEconsensus/TE_P2_all_insertions_and_genes.txt", header = FALSE)
## Also a table with strand information for genes to say whether insertions are upstream or downstream
fbstrand <- read.table("interchromate_meta/FBgn_strand.txt", header = FALSE)


## Merge the full Akaa TE objects with the significant RNA and ATAC results from F3, then combine with metadata
## to get symbol and strand
AkaaTEtab <- merge(Akaa_TEs, Akaa_intergen[,c(1:3,5,7,8)], by.x = "V8", by.y = "ensembl_gene_id")
AkaaTEtab2 <- merge(AkaaTEtab, Manz_intergen[,c(1:3,5,7,8)], by.x = "V8", by.y = "ensembl_gene_id", all.x = TRUE)
AkaaTEtab3 <- na.omit(AkaaTEtab2[(AkaaTEtab2$F3A_DE !="Not_DE" | AkaaTEtab2$F3A_DA !="Not_DA") & AkaaTEtab2$TE_pres.y == "0", ])
AkaaTEtab4 <- merge(AkaaTEtab3,fbstrand, by.x = "V8", by.y = "V1")

## Merge the full Manz TE objects with the significant RNA and ATAC results from F3
ManzTEtab <- merge(Manz_TEs, Manz_intergen[,c(1:3,5,7,8)], by.x = "V8", by.y = "ensembl_gene_id")
ManzTEtab2 <- merge(ManzTEtab, Akaa_intergen[,c(1:3,5,7,8)], by.x = "V8", by.y = "ensembl_gene_id", all.x = TRUE)
ManzTEtab3 <- na.omit(ManzTEtab2[(ManzTEtab2$F3M_DE !="Not_DE" | ManzTEtab2$F3M_DA !="Not_DA") & ManzTEtab2$TE_pres.y == "0", ])
ManzTEtab4 <- merge(ManzTEtab3,fbstrand, by.x = "V8", by.y = "V1")

## Merge the full Akaa TE objects with the significant RNA and ATAC results from F6
AkaaTEF6tab <- merge(Akaa_TEs, Akaa_intergen[,c(1,5,9,12)], by.x = "V8", by.y = "ensembl_gene_id")
AkaaTEF6tab2 <- merge(AkaaTEF6tab, Manz_intergen[,c(1,5,9,12)], by.x = "V8", by.y = "ensembl_gene_id", all.x = TRUE)
AkaaTEF6tab3 <- na.omit(AkaaTEF6tab2[AkaaTEF6tab2$F6A_DE !="Not_DE" & AkaaTEF6tab2$TE_pres.y == "0", ])
AkaaTEF6tab4 <- merge(AkaaTEF6tab3,fbstrand, by.x = "V8", by.y = "V1")

## Merge the full Manz TE objects with the significant RNA and ATAC results from F6
ManzTEF6tab <- merge(Manz_TEs, Manz_intergen[,c(1,5,9,12)], by.x = "V8", by.y = "ensembl_gene_id")
ManzTEF6tab2 <- merge(ManzTEF6tab, Akaa_intergen[,c(1,5,9,12)], by.x = "V8", by.y = "ensembl_gene_id", all.x = TRUE)
ManzTEF6tab3 <- na.omit(ManzTEF6tab2[ManzTEF6tab2$F6M_DE !="Not_DE" & ManzTEF6tab2$TE_pres.y == "0", ])
ManzTEF6tab4 <- merge(ManzTEF6tab3,fbstrand, by.x = "V8", by.y = "V1")

## Class TE as upstream/downstream (seperate or overlapping or intronic)
TE_sumlist <- list(AkaaTEtab4,ManzTEtab4,AkaaTEF6tab4,ManzTEF6tab4)
inserttype <- function(x) 
{
  x$insertion <- 
    ifelse((x$V3.x < x$V10 & x$V3.y == "+"), "upstream_sprt", 
       ifelse((x$V3.x < x$V10 & x$V3.y == "-"),"downstream_sprt", 
          ifelse((x$V2.x > x$V11 & x$V3.y == "+"),"downstream_sprt", 
             ifelse((x$V2.x > x$V11 & x$V3.y == "-"),"upstream_sprt",
                ifelse((x$V2.x < x$V10 & x$V3.x > x$V10 & x$V3.y == "+"), "upstream_ovlp", 
                   ifelse((x$V2.x < x$V10 & x$V3.x > x$V10 & x$V3.y == "-"),"downstream_ovlp", 
                      ifelse((x$V2.x < x$V11 & x$V3.x > x$V11 & x$V3.y == "+"),"downstream_ovlp", 
                         ifelse((x$V2.x < x$V11 & x$V3.x > x$V11 & x$V3.y == "-"),"upstream_ovlp",
                            ifelse((x$V2.x > x$V10 & x$V3.x < x$V11),"intronic", "unclear")))))))))
return(x)
}

TE_sumlist_updated <- lapply(TE_sumlist, inserttype)

AkaaTEtab5 <- TE_sumlist_updated[[1]][,c(7,2:5,1,22,24,23,12,15,13,16,17,20,18,21,6)]
ManzTEtab5 <- TE_sumlist_updated[[2]][,c(7,2:5,1,22,24,23,12,15,13,16,17,20,18,21,6)]
AkaaTEF6tab5 <- TE_sumlist_updated[[3]][,c(7,2:5,1,18,20,19,13,14,16,17)]
ManzTEF6tab5 <- TE_sumlist_updated[[4]][,c(7,2:5,1,18,20,19,13,14,16,17)]

# Save these output tables
write.table(AkaaTEtab5,file = "interchromate_downstream/Akaa_TE_F3_tab.txt", sep = "\t")
write.table(ManzTEtab5,file = "interchromate_downstream/Manz_TE_F3_tab.txt", sep = "\t")
write.table(AkaaTEF6tab5,file = "interchromate_downstream/Akaa_TE_F6_tab.txt", sep = "\t")
write.table(ManzTEF6tab5,file = "interchromate_downstream/Manz_TE_F6_tab.txt", sep = "\t")


################################################################
################################################################
# 4. Euler plots to visualise set analysis
################################################################
################################################################

# load required packages
library(eulerr)

########################################################
# 4.1 Overlap between heat shock DEGs   

##### F3 pop comp
plot(euler(LIST_rna1.ids[c(1,2)], shape= "circle"),quantities = TRUE, fill =c("#6295cd", "#cb4f42"),lty =c(1,1), labels = c("Ak-F3", "Ma-F3"))
##### F6 pop comp
plot(euler(LIST_rna3.ids[c(3,4)], shape= "ellipse"),quantities = TRUE, fill =c("#BED8F4", "#FFA9A0"),lty =c(1,1), labels = c("Ak-F6", "Ma-F6"))
##### Akaa gen comp
plot(euler(LIST_rna3.ids[c(11,3)], shape= "ellipse"),quantities = TRUE, fill =c("#6295cd", "#BED8F4"),lty =c(1,1), labels = c("Ak-F3", "Ak-F6"))
##### Manz gen comp
plot(euler(LIST_rna3.ids[c(12,4)], shape= "ellipse"),quantities = TRUE, fill =c("#cb4f42", "#FFA9A0"),lty =c(1,1), labels = c("Ma-F3", "Ma-F6"))

########################################################
# 4.2 Overlap between heat shock DARs   

##### F3 pop comp
plot(euler(LIST_atac1.ids[c(1,2)], shape= "circle"),quantities = TRUE, fill =c("#6295cd", "#cb4f42"),lty =c(1,1), labels = c("Ak-F3", "Ma-F3"))
##### F6 pop comp
plot(euler(LIST_atac3.ids[c(3,4)], shape= "ellipse"),quantities = TRUE, fill =c("#BED8F4", "#FFA9A0"),lty =c(1,1), labels = c("Ak-F6", "Ma-F6"))
##### Akaa gen comp
plot(euler(LIST_atac3.ids[c(11,3)], shape= "ellipse"),quantities = TRUE, fill =c("#6295cd", "#BED8F4"),lty =c(1,1), labels = c("Ak-F3", "Ak-F6"))
##### Manz gen comp
plot(euler(LIST_atac3.ids[c(12,4)], shape= "ellipse"),quantities = TRUE, fill =c("#cb4f42", "#FFA9A0"),lty =c(1,1), labels = c("Ma-F3", "Ma-F6"))

########################################################
# 4.3 overlap between DEGs and DARs

##### F3 Akaa HS DE/DA intersection
plot(euler(LIST_rna_atac2.ids[c(1,2)], shape= "ellipse"),quantities = TRUE, fill =c("#3762CF", "#2C9AC9"),lty =c(1,1), labels = c("Ak F3 DE RNA", "Ak F3 DA ATAC"))
# overlap of 449
##### F3 Manz HS DE/DA intersection
plot(euler(LIST_rna_atac2.ids[c(3,4)], shape= "ellipse"),quantities = TRUE, fill =c("#DD657D", "#ED956D"),lty =c(1,1), labels = c("Ma F3 DE RNA", "Ma F3 DA ATAC"))
# overlap of 426
##### Overlap between populations (all F3 HS DE/DA)
plot(euler(LIST_rna_atac1.ids[c(1,2)], shape= "ellipse"),quantities = TRUE, fill =c("#6295cd", "#cb4f42"),lty =c(1,1), labels = c("Ak F3 DE+DA", "Ma F3 DE+DA"))

################################################################
################################################################
# 5. F3 expression versus accessibility (rna-seq vs atac-seq) 
################################################################
################################################################

# load required packages
library(ggplot2)
library(ggblend)
library(ggrepel)

### This function allows clearer plotting of dense data (e.g. overlapping values)
trace(grDevices:::png, quote({
  if (missing(type) && missing(antialias)) {
    type <- "cairo-png"
    antialias <- "subpixel"
  }
}), print = FALSE)

################################################################
# 5.1 model showing correlation between rna and atac, with DE class 
#     (not DE, shared, unique) as an interaction.
Akaa_lm<- lm(RNAseq.log2FoldChange.x ~ atac.log2FoldChange.x*F3A_DE, data = Akaa_intergen)
summary(Akaa_lm)
anova(Akaa_lm)

Manz_lm<- lm(RNAseq.log2FoldChange.x ~ atac.log2FoldChange.x*F3M_DE, data = Manz_intergen)
summary(Manz_lm)
anova(Manz_lm)

################################################################
# 5.2 plots of rna-seq lfc vs atac lfc

# gene symbols to be included in figure
Akaa_intergen$name = mapIds(org.Dm.eg.db, keys=Akaa_intergen$ensembl_gene_id, column="SYMBOL", keytype="FLYBASE",multiVals="first")
Manz_intergen$name = mapIds(org.Dm.eg.db, keys=Manz_intergen$ensembl_gene_id, column="SYMBOL",keytype="FLYBASE",multiVals="first")

# ggplot 2 theme to be used in figure
plot_theme1 <- theme(axis.title.x = element_text(size=14),
                    axis.title.y = element_text(size=14),
                    axis.text.x = element_text(size=12),
                    axis.text.y = element_text(size=12),
                    legend.title = element_text(size=14),
                    legend.text = element_text(size=12),
                    legend.position="bottom",
                    legend.direction="vertical")

# Akaa rna vs atac
ggplot(Akaa_intergen, aes(x=atac.log2FoldChange.x, y=RNAseq.log2FoldChange.x, color=as.factor(F3A_DE), fill=as.factor(F3A_DE), shape=as.factor(F3A_DE), label = name)) + geom_point(size = 2, alpha = 0.3) * (blend("lighten") + blend("multiply", alpha = 0.5)) +
  geom_smooth(method=lm, alpha = 0.3,aes(color=as.factor(F3A_DE))) +
  geom_text_repel(aes(label=ifelse((z.atac.x>17.996 & F3A_DE != "Not_DE"),as.character(name),'')),  point.padding = 0.3, box.padding = 0.3)+
  #geom_text(aes(label=ifelse(z.atac.x>17.996,as.character(name),'')),hjust=0,vjust=0)+
  scale_color_manual(values = c("grey70","#C73EE4", "#FFA639"), labels = c("Not DEG", "Shared", "Unique"))+
  scale_fill_manual(values = c("grey70","#C73EE4", "#FFA639"), labels = c("Not DEG", "Shared", "Unique"))+
  scale_shape_manual(values = c(1,15,16), labels = c("Not DEG", "Shared", "Unique"))+
  geom_hline(yintercept=0)+ geom_vline(xintercept=0)+
  coord_cartesian(xlim=c(-1.1,2.4), ylim=c(-7.6,14))+
  labs(fill = "DE Class", color = "DE Class", shape = "DE Class", y = "RNA LFC", x = "ATAC LFC")+
  theme_bw()+
  plot_theme1
# warning messages associated with blending (overlapping points), but plots appear correct

# Manz rna vs atac
ggplot(Manz_intergen, aes(x=atac.log2FoldChange.x, y=RNAseq.log2FoldChange.x, color=as.factor(F3M_DE), fill=as.factor(F3M_DE), shape=as.factor(F3M_DE), label = name)) + geom_point(size = 2, alpha = 0.3) * (blend("lighten") + blend("multiply", alpha = 0.5)) +
  geom_smooth(method=lm, alpha = 0.3,aes(color=as.factor(F3M_DE))) +
  geom_text_repel(aes(label=ifelse((z.atac.x>12 & F3M_DE != "Not_DE"),as.character(name),'')),  point.padding = 0.3, box.padding = 0.3)+
  scale_color_manual(values = c("grey70","#C73EE4", "#FFA639"), labels = c("Not DEG", "Shared", "Unique"))+
  scale_fill_manual(values = c("grey70","#C73EE4", "#FFA639"), labels = c("Not DEG", "Shared", "Unique"))+
  scale_shape_manual(values = c(1,15,16), labels = c("Not DEG", "Shared", "Unique"))+
  geom_hline(yintercept=0)+ geom_vline(xintercept=0)+
  coord_cartesian(xlim=c(-1.1,2.4), ylim=c(-7.6,14))+
  labs(fill = "DE Class", color = "DE Class", shape = "DE Class", y = "RNA LFC", x = "ATAC LFC")+
  theme_bw()+
  plot_theme1
# warning messages associated with blending (overlapping points), but plots appear correct

################################################################
################################################################
# 6. Enrichment Analysis of genes that were both DE and DA 
################################################################
################################################################

# load required  packages
library(clusterProfiler)
library(rrvgo)

################################################################
# 6.1 First check percentage of genes with  regulation/
#      misregulation and their direction
#      DEDA = Differential Expression, Differential Accessibility
#      AS = Akaa Shared, AU = Akaa Unique
#      MS = Manz Shared, MU = Manz Unique,  

DEDA_AS<- Akaa_intergen[Akaa_intergen$ensembl_gene_id%in%LIST_rna_atac2.ids[['core']], c(2,3)]
(nrow(DEDA_AS[DEDA_AS$RNAseq.log2FoldChange.x>0 & DEDA_AS$atac.log2FoldChange.x>0,])/nrow(DEDA_AS))*100
(nrow(DEDA_AS[DEDA_AS$RNAseq.log2FoldChange.x<0 & DEDA_AS$atac.log2FoldChange.x<0,])/nrow(DEDA_AS))*100
(nrow(DEDA_AS[DEDA_AS$RNAseq.log2FoldChange.x>0 & DEDA_AS$atac.log2FoldChange.x<0,])/nrow(DEDA_AS))*100
(nrow(DEDA_AS[DEDA_AS$RNAseq.log2FoldChange.x<0 & DEDA_AS$atac.log2FoldChange.x>0,])/nrow(DEDA_AS))*100

DEDA_AU<- Akaa_intergen[Akaa_intergen$ensembl_gene_id%in%LIST_rna_atac2.ids[['Ak_reg']], c(2,3)]
(nrow(DEDA_AU[DEDA_AU$RNAseq.log2FoldChange.x>0 & DEDA_AU$atac.log2FoldChange.x>0,])/nrow(DEDA_AU))*100
(nrow(DEDA_AU[DEDA_AU$RNAseq.log2FoldChange.x<0 & DEDA_AU$atac.log2FoldChange.x<0,])/nrow(DEDA_AU))*100
(nrow(DEDA_AU[DEDA_AU$RNAseq.log2FoldChange.x>0 & DEDA_AU$atac.log2FoldChange.x<0,])/nrow(DEDA_AU))*100
(nrow(DEDA_AU[DEDA_AU$RNAseq.log2FoldChange.x<0 & DEDA_AU$atac.log2FoldChange.x>0,])/nrow(DEDA_AU))*100

DEDA_MS<- Manz_intergen[Manz_intergen$ensembl_gene_id%in%LIST_rna_atac2.ids[['core']], c(2,3)]
(nrow(DEDA_MS[DEDA_MS$RNAseq.log2FoldChange.x>0 & DEDA_MS$atac.log2FoldChange.x>0,])/nrow(DEDA_MS))*100
(nrow(DEDA_MS[DEDA_MS$RNAseq.log2FoldChange.x<0 & DEDA_MS$atac.log2FoldChange.x<0,])/nrow(DEDA_MS))*100
(nrow(DEDA_MS[DEDA_MS$RNAseq.log2FoldChange.x>0 & DEDA_MS$atac.log2FoldChange.x<0,])/nrow(DEDA_MS))*100
(nrow(DEDA_MS[DEDA_MS$RNAseq.log2FoldChange.x<0 & DEDA_MS$atac.log2FoldChange.x>0,])/nrow(DEDA_MS))*100

DEDA_MU<- Manz_intergen[Manz_intergen$ensembl_gene_id%in%LIST_rna_atac2.ids[['Ma_reg']], c(2,3)]
(nrow(DEDA_MU[DEDA_MU$RNAseq.log2FoldChange.x>0 & DEDA_MU$atac.log2FoldChange.x>0,])/nrow(DEDA_MU))*100
(nrow(DEDA_MU[DEDA_MU$RNAseq.log2FoldChange.x<0 & DEDA_MU$atac.log2FoldChange.x<0,])/nrow(DEDA_MU))*100
(nrow(DEDA_MU[DEDA_MU$RNAseq.log2FoldChange.x>0 & DEDA_MU$atac.log2FoldChange.x<0,])/nrow(DEDA_MU))*100
(nrow(DEDA_MU[DEDA_MU$RNAseq.log2FoldChange.x<0 & DEDA_MU$atac.log2FoldChange.x>0,])/nrow(DEDA_MU))*100
# Shared genes show identical patterns of regulation
# Results table = Table 1 of manuscript

################################################################
# 6.2 GO enrichment (Biological Processes) for the DE/DA genes 
#     by reg/misreg direction

# select shared and unique DEDA genes from the two populations
#     AS = Akaa Shared, AU = Akaa Unique, 
#     MS = Manz Shared, MU = Manz Unique
AS_DEDA <- Akaa_intergen[Akaa_intergen$ensembl_gene_id%in%LIST_rna_atac2.ids[['core']], c(1:3)]
AU_DEDA <- Akaa_intergen[Akaa_intergen$ensembl_gene_id%in%LIST_rna_atac2.ids[['Ak_reg']], c(1:3)]
MS_DEDA <- Manz_intergen[Manz_intergen$ensembl_gene_id%in%LIST_rna_atac2.ids[['core']], c(1:3)]
MU_DEDA <- Manz_intergen[Manz_intergen$ensembl_gene_id%in%LIST_rna_atac2.ids[['Ma_reg']], c(1:3)]
# AS and MS DEDA should be identical for Manz and Akaa (given results of code chunk 6.1), but check

## Get list of all gene ids that appear in our analysis. This 
# will be 'universe' for GO enrichment analysis
DmelBG_A <- Akaa_intergen$ensembl_gene_id
DmelBG_M <- Manz_intergen$ensembl_gene_id
DmelBG_All<-unique(sort(c(DmelBG_A,DmelBG_M)))

# GO enrichment naming conventions:
# CU = Chromatin (Accessibility) Up, CD = Chromatin (Accessibility) Down
# EU = Expression Up, ED = Expression Down

#################################
# 6.2.1 Akaa Shared GO enrichment
# CUEU
GO_AS_CUEU <- enrichGO(gene = AS_DEDA[AS_DEDA$RNAseq.log2FoldChange.x>0 & AS_DEDA$atac.log2FoldChange.x>0,1] ,
                       OrgDb = org.Dm.eg.db,keyType = 'FLYBASE', readable = T, ont = "BP",
                       universe = DmelBG_All, pvalueCutoff = 0.05, qvalueCutoff = 0.05)
# How many GO terms?
sum(GO_AS_CUEU@result$p.adjust < 0.05)
# [1] 22
# extract adjusted p values from results for revigo
GO_AS_CUEU_scores <- setNames(-log10(GO_AS_CUEU$p.adjust),GO_AS_CUEU$ID)
# use revigo to calculate similarity between significant GO terms
GO_AS_CUEU_simtrx <- calculateSimMatrix(GO_AS_CUEU$ID, orgdb="org.Dm.eg.db",ont="BP",method="Rel")
# Cluster GO terms into smaller list based on similarity and q values
GO_AS_CUEU_redtrm <- reduceSimMatrix(GO_AS_CUEU_simtrx,GO_AS_CUEU_scores,threshold=0.7,orgdb="org.Dm.eg.db")
# plot grouped GO terms following revigo summarising
treemapPlot(GO_AS_CUEU_redtrm)
# Write table t of simplified GO terms 
GO_AS_CUEU_merge<-merge(GO_AS_CUEU@result[GO_AS_CUEU@result$p.adjust < 0.05,], GO_AS_CUEU_redtrm, by.x = "ID", by.y = "go")
GO_AS_CUEU_ord<-GO_AS_CUEU_merge[with(GO_AS_CUEU_merge, order(GO_AS_CUEU_merge$cluster, -GO_AS_CUEU_merge$parentSimScore)),]
GO_AS_CUEU_out = GO_AS_CUEU_ord[!duplicated(GO_AS_CUEU_ord$cluster),]
write.table(GO_AS_CUEU_out,file = "interchromate_GOenrich/GO_AS_CUEU_out.txt", sep = "\t")

# CDED
GO_AS_CDED <- enrichGO(gene = AS_DEDA[AS_DEDA$RNAseq.log2FoldChange.x<0 & AS_DEDA$atac.log2FoldChange.x<0,1] ,
                       OrgDb = org.Dm.eg.db,keyType = 'FLYBASE', readable = T, ont = "BP",
                       universe = DmelBG_All, pvalueCutoff = 0.05, qvalueCutoff = 0.05)
sum(GO_AS_CDED@result$p.adjust < 0.05)
#[1] 0

# CUED
GO_AS_CUED <- enrichGO(gene = AS_DEDA[AS_DEDA$RNAseq.log2FoldChange.x<0 & AS_DEDA$atac.log2FoldChange.x>0,1] ,
                       OrgDb = org.Dm.eg.db,keyType = 'FLYBASE', readable = T, ont = "BP",
                       universe = DmelBG_All, pvalueCutoff = 0.05, qvalueCutoff = 0.05)
sum(GO_AS_CUED@result$p.adjust < 0.05)
#[1] 0

# CDEU
GO_AS_CDEU <- enrichGO(gene = AS_DEDA[AS_DEDA$RNAseq.log2FoldChange.x>0 & AS_DEDA$atac.log2FoldChange.x<0,1] ,
                       OrgDb = org.Dm.eg.db,keyType = 'FLYBASE', readable = T, ont = "BP",
                       universe = DmelBG_All, pvalueCutoff = 0.05, qvalueCutoff = 0.05)
sum(GO_AS_CDEU@result$p.adjust < 0.05)
#[1] 1
write.table(GO_AS_CDEU,file = "interchromate_GOenrich/GO_AS_CDEU_out.txt", sep = "\t")

#################################
# 6.2.2 Manz Shared GO enrichment should be identical
# CUEU
GO_MS_CUEU <- enrichGO(gene = MS_DEDA[MS_DEDA$RNAseq.log2FoldChange.x>0 & MS_DEDA$atac.log2FoldChange.x>0,1] ,
                       OrgDb = org.Dm.eg.db,keyType = 'FLYBASE', readable = T, ont = "BP",
                       universe = DmelBG_All, pvalueCutoff = 0.05, qvalueCutoff = 0.05)
sum(GO_MS_CUEU@result$p.adjust < 0.05)
#[1] 22
# Are Manz and Akaa results identical?
identical(GO_MS_CUEU,GO_AS_CUEU)
# [1] TRUE
GO_MS_CDED <- enrichGO(gene = MS_DEDA[MS_DEDA$RNAseq.log2FoldChange.x<0 & MS_DEDA$atac.log2FoldChange.x<0,1] ,
                       OrgDb = org.Dm.eg.db,keyType = 'FLYBASE', readable = T, ont = "BP",
                       universe = DmelBG_All, pvalueCutoff = 0.05, qvalueCutoff = 0.05)
sum(GO_MS_CDED@result$p.adjust < 0.05)
#[1] 0
GO_MS_CUED <- enrichGO(gene = MS_DEDA[MS_DEDA$RNAseq.log2FoldChange.x<0 & MS_DEDA$atac.log2FoldChange.x>0,1] ,
                       OrgDb = org.Dm.eg.db,keyType = 'FLYBASE', readable = T, ont = "BP",
                       universe = DmelBG_All, pvalueCutoff = 0.05, qvalueCutoff = 0.05)
sum(GO_MS_CUED@result$p.adjust < 0.05)
#[1] 0
GO_MS_CDEU <- enrichGO(gene = MS_DEDA[MS_DEDA$RNAseq.log2FoldChange.x>0 & MS_DEDA$atac.log2FoldChange.x<0,1] ,
                       OrgDb = org.Dm.eg.db,keyType = 'FLYBASE', readable = T, ont = "BP",
                       universe = DmelBG_All, pvalueCutoff = 0.05, qvalueCutoff = 0.05)
sum(GO_MS_CDEU@result$p.adjust < 0.05)
# Are Manz and Akaa results identical?
identical(GO_MS_CDEU,GO_AS_CDEU)
# [1] TRUE

#################################
# 6.2.3 Akaa Unique
# CUEU
GO_AU_CUEU <- enrichGO(gene = AU_DEDA[AU_DEDA$RNAseq.log2FoldChange.x>0 & AU_DEDA$atac.log2FoldChange.x>0,1] ,
                       OrgDb = org.Dm.eg.db,keyType = 'FLYBASE', readable = T, ont = "BP",
                       universe = DmelBG_All, pvalueCutoff = 0.05, qvalueCutoff = 0.05)
sum(GO_AU_CUEU@result$p.adjust < 0.05)
# 0

# CDED
GO_AU_CDED <- enrichGO(gene = AU_DEDA[AU_DEDA$RNAseq.log2FoldChange.x<0 & AU_DEDA$atac.log2FoldChange.x<0,1] ,
                       OrgDb = org.Dm.eg.db,keyType = 'FLYBASE', readable = T, ont = "BP",
                       universe = DmelBG_All, pvalueCutoff = 0.05, qvalueCutoff = 0.05)
sum(GO_AU_CDED@result$p.adjust < 0.05)
# 19
GO_AU_CDED_scores <- setNames(-log10(GO_AU_CDED$p.adjust),GO_AU_CDED$ID)
GO_AU_CDED_simtrx <- calculateSimMatrix(GO_AU_CDED$ID, orgdb="org.Dm.eg.db",ont="BP",method="Rel")
GO_AU_CDED_redtrm <- reduceSimMatrix(GO_AU_CDED_simtrx,GO_AU_CDED_scores,threshold=0.7,orgdb="org.Dm.eg.db")
treemapPlot(GO_AU_CDED_redtrm)
GO_AU_CDED_merge<-merge(GO_AU_CDED@result[GO_AU_CDED@result$p.adjust < 0.05,], GO_AU_CDED_redtrm, by.x = "ID", by.y = "go")
GO_AU_CDED_ord<-GO_AU_CDED_merge[with(GO_AU_CDED_merge, order(GO_AU_CDED_merge$cluster, -GO_AU_CDED_merge$parentSimScore)),]
GO_AU_CDED_out = GO_AU_CDED_ord[!duplicated(GO_AU_CDED_ord$cluster),]
write.table(GO_AU_CDED_out,file = "interchromate_GOenrich/GO_AU_CDED_out.txt", sep = "\t")

# CUED
GO_AU_CUED <- enrichGO(gene = AU_DEDA[AU_DEDA$RNAseq.log2FoldChange.x<0 & AU_DEDA$atac.log2FoldChange.x>0,1] ,
                       OrgDb = org.Dm.eg.db,keyType = 'FLYBASE', readable = T, ont = "BP",
                       universe = DmelBG_All, pvalueCutoff = 0.05, qvalueCutoff = 0.05)
sum(GO_AU_CUED@result$p.adjust < 0.05)
# 0

#CDEU
GO_AU_CDEU <- enrichGO(gene = AU_DEDA[AU_DEDA$RNAseq.log2FoldChange.x>0 & AU_DEDA$atac.log2FoldChange.x<0,1] ,
                       OrgDb = org.Dm.eg.db,keyType = 'FLYBASE', readable = T, ont = "BP",
                       universe = DmelBG_All, pvalueCutoff = 0.05, qvalueCutoff = 0.05)
sum(GO_AU_CDEU@result$p.adjust < 0.05)
# 28
GO_AU_CDEU_scores <- setNames(-log10(GO_AU_CDEU$p.adjust),GO_AU_CDEU$ID)
GO_AU_CDEU_simtrx <- calculateSimMatrix(GO_AU_CDEU$ID, orgdb="org.Dm.eg.db",ont="BP",method="Rel")
GO_AU_CDEU_redtrm <- reduceSimMatrix(GO_AU_CDEU_simtrx,GO_AU_CDEU_scores,threshold=0.7,orgdb="org.Dm.eg.db")
treemapPlot(GO_AU_CDEU_redtrm)
GO_AU_CDEU_merge<-merge(GO_AU_CDEU@result[GO_AU_CDEU@result$p.adjust < 0.05,], GO_AU_CDEU_redtrm, by.x = "ID", by.y = "go")
GO_AU_CDEU_ord<-GO_AU_CDEU_merge[with(GO_AU_CDEU_merge, order(GO_AU_CDEU_merge$cluster, -GO_AU_CDEU_merge$parentSimScore)),]
GO_AU_CDEU_out = GO_AU_CDEU_ord[!duplicated(GO_AU_CDEU_ord$cluster),]
write.table(GO_AU_CDEU_out,file = "interchromate_GOenrich/GO_AU_CDEU_out.txt", sep = "\t")

#################################
# 6.2.4 Manz Unique
# CUEU
GO_MU_CUEU <- enrichGO(gene = MU_DEDA[MU_DEDA$RNAseq.log2FoldChange.x>0 & MU_DEDA$atac.log2FoldChange.x>0,1] ,
                       OrgDb = org.Dm.eg.db,keyType = 'FLYBASE', readable = T, ont = "BP",
                       universe = DmelBG_All, pvalueCutoff = 0.05, qvalueCutoff = 0.05)
sum(GO_MU_CUEU@result$p.adjust < 0.05)
# 0

# CDED
GO_MU_CDED <- enrichGO(gene = MU_DEDA[MU_DEDA$RNAseq.log2FoldChange.x<0 & MU_DEDA$atac.log2FoldChange.x<0,1] ,
                       OrgDb = org.Dm.eg.db,keyType = 'FLYBASE', readable = T, ont = "BP",
                       universe = DmelBG_All, pvalueCutoff = 0.05, qvalueCutoff = 0.05)
sum(GO_MU_CDED@result$p.adjust < 0.05)
# 0

# CUED
GO_MU_CUED <- enrichGO(gene = MU_DEDA[MU_DEDA$RNAseq.log2FoldChange.x<0 & MU_DEDA$atac.log2FoldChange.x>0,1] ,
                       OrgDb = org.Dm.eg.db,keyType = 'FLYBASE', readable = T, ont = "BP",
                       universe = DmelBG_All, pvalueCutoff = 0.05, qvalueCutoff = 0.05)
sum(GO_MU_CUED@result$p.adjust < 0.05)
# 0

#CDEU
GO_MU_CDEU <- enrichGO(gene = MU_DEDA[MU_DEDA$RNAseq.log2FoldChange.x>0 & MU_DEDA$atac.log2FoldChange.x<0,1] ,
                       OrgDb = org.Dm.eg.db,keyType = 'FLYBASE', readable = T, ont = "BP",
                       universe = DmelBG_All, pvalueCutoff = 0.05, qvalueCutoff = 0.05)
sum(GO_MU_CDEU@result$p.adjust < 0.05)
# 9
GO_MU_CDEU_scores <- setNames(-log10(GO_MU_CDEU$p.adjust),GO_MU_CDEU$ID)
GO_MU_CDEU_simtrx <- calculateSimMatrix(GO_MU_CDEU$ID, orgdb="org.Dm.eg.db",ont="BP",method="Rel")
GO_MU_CDEU_redtrm <- reduceSimMatrix(GO_MU_CDEU_simtrx,GO_MU_CDEU_scores,threshold=0.7,orgdb="org.Dm.eg.db")
treemapPlot(GO_MU_CDEU_redtrm)
GO_MU_CDEU_merge<-merge(GO_MU_CDEU@result[GO_MU_CDEU@result$p.adjust < 0.05,], GO_MU_CDEU_redtrm, by.x = "ID", by.y = "go")
GO_MU_CDEU_ord<-GO_MU_CDEU_merge[with(GO_MU_CDEU_merge, order(GO_MU_CDEU_merge$cluster, -GO_MU_CDEU_merge$parentSimScore)),]
GO_MU_CDEU_out = GO_MU_CDEU_ord[!duplicated(GO_MU_CDEU_ord$cluster),]
write.table(GO_MU_CDEU_out,file = "interchromate_GOenrich/GO_MU_CDEU_out.txt", sep = "\t")


################################################################
################################################################
# 7. Compare gene expression across F3 & F6 (transgenerational)
################################################################
################################################################

# No addtional packages required

################################################################
# 7.1 Select genes from rna list that show overlap between 
#     F3 and F6 in Akaa and Manz (transgenerational)

# Akaa 
ovlpA_intergen <- Akaa_intergen[Akaa_intergen$ensembl_gene_id %in% LIST_rna2.ids[['F3A_F6A_ovlp_NoGenA']], ]
nrow(ovlpA_intergen)
# [1] 5
# Too few transgeneerational genes in Akaa to do  test
# # Gene IDs that are transgen in Manz, but lfc values taken from Akaa
ovlpAinM_intergen <- Manz_intergen[Manz_intergen$ensembl_gene_id %in% LIST_rna2.ids[['F3A_F6A_ovlp_NoGenA']], ]

# Focus on Manz
ovlpM_intergen <- Manz_intergen[Manz_intergen$ensembl_gene_id %in% LIST_rna2.ids[['F3M_F6M_ovlp_NoGenM']], ]
# Gene IDs that are transgen in Manz, but lfc values taken from Akaa
ovlpMinA_intergen <- Akaa_intergen[Akaa_intergen$ensembl_gene_id %in% LIST_rna2.ids[['F3M_F6M_ovlp_NoGenM']], ]

# Generate list of genes DE in both F3 and F6 Akaa, and associated with TEs
TGA_Akaa_TE <- merge(ovlpA_intergen, Akaa_TEs, by.x = "ensembl_gene_id", by.y = "V8", all.x = TRUE)
TGA_Akaa_Ma_TE <- merge(ovlpAinM_intergen, Manz_TEs, by.x = "ensembl_gene_id", by.y = "V8", all.x = TRUE)
TGA_Akaa_TE2 <- TGA_Akaa_TE[TGA_Akaa_TE$TE_pres !="0", ]
TGA_Akaa_TE3 <- merge(TGA_Akaa_TE2[,c(1:5,7:10,12,17:22)],TGA_Akaa_Ma_TE[,c(1:5,7:10,12,17:22)], by = "ensembl_gene_id", all.x = TRUE )
TGA_Akaa_TE4 <- merge(TGA_Akaa_TE3,fbstrand, by.x = "ensembl_gene_id", by.y = "V1")

# Generate list of genes DE in both F3 and F6 Manz, and associated with TEs
TGA_Manz_TE <- merge(ovlpM_intergen, Manz_TEs, by.x = "ensembl_gene_id", by.y = "V8", all.x = TRUE)
TGA_Manz_Ak_TE <- merge(ovlpMinA_intergen, Akaa_TEs, by.x = "ensembl_gene_id", by.y = "V8", all.x = TRUE)
TGA_Manz_TE2 <- TGA_Manz_TE[TGA_Manz_TE$TE_pres !="0", ]
TGA_Manz_TE3 <- merge(TGA_Manz_TE2[,c(1:5,7:10,12,17:22)],TGA_Manz_Ak_TE[,c(1:5,7:10,12,17:22)], by = "ensembl_gene_id", all.x = TRUE )
TGA_Manz_TE4 <- merge(TGA_Manz_TE3,fbstrand, by.x = "ensembl_gene_id", by.y = "V1")

# Optional: save these objects 
write.table(TGA_Akaa_TE4,file = "interchromate_downstream/Akaa_TEs_F3F6transgen.txt", sep = "\t")
write.table(TGA_Manz_TE4,file = "interchromate_downstream/Manz_TEs_F3F6transgen.txt", sep = "\t")

################################################################
# 7.2 linear models: regression of F3 expr. LFC vs F6 expr. LFC

# Manzanares F3 vs F6
TG_mM<- lm(RNAseq.log2FoldChange.y ~ RNAseq.log2FoldChange.x*F3M_DE, data = ovlpM_intergen)
anova(TG_mM)
summary(TG_mM)
TG_mM_a<- lm(RNAseq.log2FoldChange.y ~ RNAseq.log2FoldChange.x, data = ovlpM_intergen)
anova(TG_mM_a)
summary(TG_mM_a)
#                         Estimate   Std. Error  t value Pr(>|t|)    
# (Intercept)             0.0009858  0.0219475   0.045    0.964    
# RNAseq.log2FoldChange.x 0.3463169  0.0349393   9.912   <2e-16 ***

# Gene IDs from Manz, LFC values from Akaa (MinA)
TG_mMinA<- lm(RNAseq.log2FoldChange.y ~ RNAseq.log2FoldChange.x*F3A_DE, data = ovlpMinA_intergen)
anova(TG_mMinA)
summary(TG_mMinA)
TG_mMinA_a<- lm(RNAseq.log2FoldChange.y ~ RNAseq.log2FoldChange.x, data = ovlpMinA_intergen)
anova(TG_mMinA_a)
summary(TG_mMinA_a)
#                         Estimate  Std. Error t value Pr(>|t|)    
# (Intercept)             -0.01417    0.00803  -1.764   0.0801 .  
# RNAseq.log2FoldChange.x -0.14460    0.01461  -9.898   <2e-16 ***

################################################################
# 7.3 linear models: regression of F3 expr. LFC vs F6 expr. LFC

# Colours for figures:
# Blue for Akaa: "#6295cd"
# Red for Manz:  "#cb4f42"

# ggplot theme for these figures
plot_theme2 <- theme(axis.title.x = element_text(size = 16),
                    axis.text.x = element_text(size = 14),
                    axis.text.y = element_text(size = 14),
                    axis.title.y = element_text(size = 16),
                    title = element_text(size = 20),
                    plot.title = element_text(hjust = 0.5),
                    legend.position = "none")

# plot of Manzanares transgenerational gene IDs, Akaa LFC values
plot_TG_MinA <- ggplot(ovlpMinA_intergen, aes(x=RNAseq.log2FoldChange.x, y=RNAseq.log2FoldChange.y)) + geom_point(size = 2, alpha = 0.3, color = "#6295cd") * (blend("lighten") + blend("multiply", alpha = 0.5)) +
  geom_smooth(method=lm, alpha = 0.3,aes(fill = "#6295cd"), color = "#6295cd") +
  scale_color_manual(values = c("#6295cd"))+
  scale_fill_manual(values = c("#6295cd"))+
  scale_shape_manual(values = c(16))+
  geom_hline(yintercept=0)+ geom_vline(xintercept=0)+
  coord_cartesian(xlim=c(-1.3,1.2), ylim=c(-1.15,1.25))+
  labs(fill = "Response type", color = "Response type", shape = "Response type",
       y = "RNA LFC in F6", x = "RNA LFC in F3", title = "Akaa")+
  theme_bw()+   
  plot_theme2

# plot of Manzanares transgenerational genes
plot_TG_M <-ggplot(ovlpM_intergen, aes(x=RNAseq.log2FoldChange.x, y=RNAseq.log2FoldChange.y)) + geom_point(size = 2, alpha = 0.3, color = "#cb4f42") * (blend("lighten") + blend("multiply", alpha = 0.5)) +
  geom_smooth(method=lm, alpha = 0.3,aes(fill = "#cb4f42"), color = "#cb4f42") +
  scale_color_manual(values = c("#cb4f42"))+
  scale_fill_manual(values = c("#cb4f42"))+
  scale_shape_manual(values = c(16))+
  geom_hline(yintercept=0)+ geom_vline(xintercept=0)+
  coord_cartesian(xlim=c(-1.3,1.3), ylim=c(-1.15,1.25))+
  labs(fill = "Response type", color = "Response type", shape = "Response type",
       y = "RNA LFC in F6", x = "RNA LFC in F3", title = "Manz")+
  theme_bw()+
  plot_theme2

################################################################
# 7.4 Gene enrichement for the Manzanares transgen genes
# 

# Subset Manz data by directions of F3 and F6 LFC
M_intergen_3up6up <- subset(Manz_intergen, (RNAseq.log2FoldChange.x > 0 & RNAseq.log2FoldChange.y > 0))
M_intergen_3up6dn <- subset(Manz_intergen, (RNAseq.log2FoldChange.x > 0 & RNAseq.log2FoldChange.y < 0))
M_intergen_3dn6up <- subset(Manz_intergen, (RNAseq.log2FoldChange.x < 0 & RNAseq.log2FoldChange.y > 0))
M_intergen_3dn6dn <- subset(Manz_intergen, (RNAseq.log2FoldChange.x < 0 & RNAseq.log2FoldChange.y < 0))

# Select gene IDs that were significantly DE in both F3 qand F6 in Manz
M_3up6up_sig <- M_intergen_3up6up[M_intergen_3up6up$ensembl_gene_id %in% LIST_rna2.ids[['F3M_F6M_ovlp_NoGenM']],]$ensembl_gene_id
M_3up6dn_sig <- M_intergen_3up6dn[M_intergen_3up6dn$ensembl_gene_id %in% LIST_rna2.ids[['F3M_F6M_ovlp_NoGenM']],]$ensembl_gene_id
M_3dn6up_sig <- M_intergen_3dn6up[M_intergen_3dn6up$ensembl_gene_id %in% LIST_rna2.ids[['F3M_F6M_ovlp_NoGenM']],]$ensembl_gene_id
M_3dn6dn_sig <- M_intergen_3dn6dn[M_intergen_3dn6dn$ensembl_gene_id %in% LIST_rna2.ids[['F3M_F6M_ovlp_NoGenM']],]$ensembl_gene_id

##########################
# 7.4.1 Manz F3 up, F6 up
M_3up6up_enrich <- enrichGO(gene = M_3up6up_sig ,
                            OrgDb = org.Dm.eg.db,keyType = 'FLYBASE', readable = T,
                            ont = "BP",
                            universe = DmelBG_All,
                            pvalueCutoff = 0.05, 
                            qvalueCutoff = 0.05)

sum(M_3up6up_enrich@result$p.adjust<0.05)
# 146
M_3up6up_enrich_scores <- setNames(-log10(M_3up6up_enrich$p.adjust), M_3up6up_enrich$ID)
M_3up6up_enrich_simtrx <- calculateSimMatrix(M_3up6up_enrich$ID, orgdb="org.Dm.eg.db",ont="BP",method="Rel")
M_3up6up_enrich_redtrm <- reduceSimMatrix(M_3up6up_enrich_simtrx,M_3up6up_enrich_scores,threshold=0.7,orgdb="org.Dm.eg.db")
treemapPlot(M_3up6up_enrich_redtrm)
M_3up6up_merge<-merge(M_3up6up_enrich@result[M_3up6up_enrich@result$p.adjust < 0.05,], M_3up6up_enrich_redtrm, by.x = "ID", by.y = "go")
M_3up6up_ord<-M_3up6up_merge[with(M_3up6up_merge, order(M_3up6up_merge$cluster, -M_3up6up_merge$parentSimScore)),]
M_3up6up_out = M_3up6up_ord[!duplicated(M_3up6up_ord$cluster),]
write.table(M_3up6up_out,file = "interchromate_GOenrich/M_3up6up_out.txt", sep = "\t")

##########################
# 7.4.2 Manz F3 up, F6 down
M_3up6dn_enrich <- enrichGO(gene = M_3up6dn_sig ,
                            OrgDb = org.Dm.eg.db,keyType = 'FLYBASE', readable = T,
                            ont = "BP",
                            universe = DmelBG_All,
                            pvalueCutoff = 0.05, 
                            qvalueCutoff = 0.05)
sum(M_3up6dn_enrich@result$p.adjust < 0.05)
# 0

##########################
#7.4.3 Manz F3 dn, F6 up
M_3dn6up_enrich <- enrichGO(gene = M_3dn6up_sig ,
                            OrgDb = org.Dm.eg.db,keyType = 'FLYBASE', readable = T,
                            ont = "BP",
                            universe = DmelBG_All,
                            pvalueCutoff = 0.05, 
                            qvalueCutoff = 0.05)
sum(M_3dn6up_enrich@result$p.adjust < 0.05)
# 13
M_3dn6up_enrich_scores <- setNames(-log10(M_3dn6up_enrich$p.adjust), M_3dn6up_enrich$ID)
M_3dn6up_enrich_simtrx <- calculateSimMatrix(M_3dn6up_enrich$ID, orgdb="org.Dm.eg.db",ont="BP",method="Rel")
M_3dn6up_enrich_redtrm <- reduceSimMatrix(M_3dn6up_enrich_simtrx,M_3dn6up_enrich_scores,threshold=0.7,orgdb="org.Dm.eg.db")
treemapPlot(M_3dn6up_enrich_redtrm)
M_3dn6up_merge<-merge(M_3dn6up_enrich@result[M_3dn6up_enrich@result$p.adjust < 0.05,], M_3dn6up_enrich_redtrm, by.x = "ID", by.y = "go")
M_3dn6up_ord<-M_3dn6up_merge[with(M_3dn6up_merge, order(M_3dn6up_merge$cluster, -M_3dn6up_merge$parentSimScore)),]
M_3dn6up_out = M_3dn6up_ord[!duplicated(M_3dn6up_ord$cluster),]
write.table(M_3dn6up_out,file = "interchromate_GOenrich/M_3dn6up_out.txt", sep = "\t")

##########################
# 7.4.4 Manz F3 dn, F6 dn
M_3dn6dn_enrich <- enrichGO(gene = M_3dn6dn_sig ,
                            OrgDb = org.Dm.eg.db,keyType = 'FLYBASE', readable = T,
                            ont = "BP",
                            universe = DmelBG_All,
                            pvalueCutoff = 0.05, 
                            qvalueCutoff = 0.05)
sum(M_3dn6dn_enrich@result$p.adjust < 0.05)
# 0


################################################################
################################################################
# 8. Associations between TEs and expression & accessibility
################################################################
################################################################

#load required packages
library(chisq.posthoc.test)

################################################################
# 8.1 Chi square tests for association between TEs and gene 
#    DE or DA status (not DE or DA, shared or unique)

# Akaa DE
chi_A_DE<-chisq.test(Akaa_intergen$TE_pres, Akaa_intergen$F3A_DE, correct=FALSE)
chi_A_DE
chi_A_DE$expected
chi_A_DE$observed
tab_A_DE<-table(Akaa_intergen$TE_pres, Akaa_intergen$F3A_DE)
chisq.posthoc.test(tab_A_DE,method = "BH") 
# X-squared = 16.659, df = 4, p-value = 0.002251
# Check text file for posthoc test results

# Manz DE
chi_M_DE<-chisq.test(Manz_intergen$TE_pres, Manz_intergen$F3M_DE, correct=FALSE)
chi_M_DE
chi_M_DE$expected
chi_M_DE$observed
tab_M_DE<-table(Manz_intergen$TE_pres, Manz_intergen$F3M_DE)
chisq.posthoc.test(tab_M_DE,method = "BH")
# X-squared = 19.202, df = 4, p-value = 0.0007174
# Check text file for posthoc test results

# Akaa DA
chi_A_DA<-chisq.test(Akaa_intergen$TE_pres, Akaa_intergen$F3A_DA, correct=FALSE)
chi_A_DA
chi_A_DA$expected
chi_A_DA$observed
tab_A_DA<-table(Akaa_intergen$TE_pres, Akaa_intergen$F3A_DA)
chisq.posthoc.test(tab_A_DA,method = "BH") 
# X-squared = 16.069, df = 4, p-value = 0.002928
# Check text file for posthoc test results

# Manz DA
chi_M_DA<-chisq.test(Manz_intergen$TE_pres, Manz_intergen$F3M_DA, correct=FALSE)
chi_M_DA
chi_M_DA$expected
chi_M_DA$observed
tab_M_DA<-table(Manz_intergen$TE_pres, Manz_intergen$F3M_DA)
chisq.posthoc.test(tab_M_DA,method = "BH") 
# X-squared = 6.0575, df = 4, p-value = 0.1949
# Check text file for posthoc test results

################################################################
# 8.2 Chi square tests for association between TEs and genes 
#     DEDA status (not DEDA, shared DEDA, unique DEDA)

# Akaa DEDA
chi_A_DEDA<-chisq.test(Akaa_intergen$TE_pres, Akaa_intergen$DEGDAG, correct=FALSE)
chi_A_DEDA$observed
chi_A_DEDA$expected
tab_A_DEDA<-table(Akaa_intergen$TE_pres, Akaa_intergen$DEGDAG)
chisq.posthoc.test(tab_A_DEDA,method = "BH") 
# X-squared = 9.8074, df = 4, p-value = 0.0438

# Manz DEDA
chi_M_DEDA<-chisq.test(Manz_intergen$TE_pres, Manz_intergen$DEGDAG, correct=FALSE)
chi_M_DEDA
chi_M_DEDA$expected
tab_M_DEDA<-table(Manz_intergen$TE_pres, Manz_intergen$DEGDAG)
chisq.posthoc.test(tab_M_DEDA,method = "BH") 
# X-squared = 3.9996, df = 4, p-value = 0.4061

# Just about significant effect for Akaa, but it's very slight, and represents an association
# with genes without TEs


################################################################
# 8.3 Chi square tests for association between TEs and direction
#     of change for DE or DA genes (non-DE or non-DA genes 
#     excluded)

# Akaa differential expression
ch_A_DE_dir<-chisq.test(Akaa_intergen[Akaa_intergen$F3A_DE!="Not_DE",]$TE_pres, Akaa_intergen[Akaa_intergen$F3A_DE!="Not_DE",]$Rdir, correct=FALSE)
ch_A_DE_dir
ch_A_DE_dir$expected
ch_A_DE_dir$observed
tab_A_DE_dir<-table(Akaa_intergen[Akaa_intergen$F3A_DE!="Not_DE",]$TE_pres, Akaa_intergen[Akaa_intergen$F3A_DE!="Not_DE",]$Rdir)
chisq.posthoc.test(tab_A_DE_dir,method = "BH") 
# X-squared = 18.07, df = 2, p-value = 0.0001191
# Check text file for posthoc test results

# Akaa differential accessibility
ch_A_DA_dir<-chisq.test(Akaa_intergen[Akaa_intergen$F3A_DA!="Not_DA",]$TE_pres, Akaa_intergen[Akaa_intergen$F3A_DA!="Not_DA",]$Adir, correct=FALSE)
ch_A_DA_dir
ch_A_DA_dir$expected
tab_A_DA_dir<-table(Akaa_intergen[Akaa_intergen$F3A_DA!="Not_DA",]$TE_pres, Akaa_intergen[Akaa_intergen$F3A_DA!="Not_DA",]$Adir)
chisq.posthoc.test(tab_A_DA_dir,method = "BH") 
# X-squared = 1.784, df = 2, p-value = 0.4098

# Manz differential expression
ch_M_DE_dir<-chisq.test(Manz_intergen[Manz_intergen$F3M_DE!="Not_DE",]$TE_pres, Manz_intergen[Manz_intergen$F3M_DE!="Not_DE",]$Rdir, correct=FALSE)
ch_M_DE_dir
ch_M_DE_dir$expected
tab_M_DE_dir<-table(Manz_intergen[Manz_intergen$F3M_DE!="Not_DE",]$TE_pres, Manz_intergen[Manz_intergen$F3M_DE!="Not_DE",]$Rdir)
chisq.posthoc.test(tab_M_DE_dir,method = "BH") 
# X-squared = 2.0498, df = 2, p-value = 0.3588


# Manz differential accessbility
ch_M_DA_dir<-chisq.test(Manz_intergen[Manz_intergen$F3M_DA!="Not_DA",]$TE_pres, Manz_intergen[Manz_intergen$F3M_DA!="Not_DA",]$Adir, correct=FALSE)
ch_M_DA_dir
ch_M_DA_dir$expected
ch_M_DA_dir$observed
tab_M_DA_dir<-table(Manz_intergen[Manz_intergen$F3M_DA!="Not_DA",]$TE_pres, Manz_intergen[Manz_intergen$F3M_DA!="Not_DA",]$Adir)
chisq.posthoc.test(tab_M_DA_dir,method = "BH") 
# X-squared = 38.539, df = 2, p-value = 4.28e-09
# Check text file for posthoc test results

################################################################
# 8.4 Chi square tests for association between TEs and direction
#     of change for DEDA genes (non DEDA genes excluded). 

# Subset to only include DEDA genes
Akaa_intergen_reg<- subset(Akaa_intergen, (DEGDAG !="0" ))
Manz_intergen_reg<- subset(Manz_intergen, (DEGDAG !="0" ))

ch_A_DEDA_dir<-chisq.test(Akaa_intergen_reg$TE_pres, Akaa_intergen_reg$Zdir, correct=FALSE)
ch_A_DEDA_dir
ch_A_DEDA_dir$expected
tab_A_DEDA_dir<-table(Akaa_intergen_reg$TE_pres, Akaa_intergen_reg$Zdir)
chisq.posthoc.test(tab_A_DEDA_dir,method = "BH") 
# X-squared = 6.4509, df = 6, p-value = 0.3746

ch_M_DEDA_dir<-chisq.test(Manz_intergen_reg$TE_pres, Manz_intergen_reg$Zdir, correct=FALSE)
ch_M_DEDA_dir
ch_M_DEDA_dir$expected
ch_M_DEDA_dir$observed
tab_M_DEDA_dir<-table(Manz_intergen_reg$TE_pres, Manz_intergen_reg$Zdir)
chisq.posthoc.test(tab_M_DEDA_dir,method = "BH") 
# X-squared = 20.637, df = 6, p-value = 0.002131

################################################################
# 8.5 Chi square tests for association between TEs and direction
#     of change for expression and accessibility considering
#     everything (not just sig DEDA)

ch_A_ALL_dir<-chisq.test(Akaa_intergen$TE_pres, Akaa_intergen$Zdir, correct=FALSE)
ch_A_ALL_dir
ch_A_ALL_dir$expected
ch_A_ALL_dir$observed
tab_A_ALL_dir<-table(Akaa_intergen$TE_pres, Akaa_intergen$Zdir)
chisq.posthoc.test(tab_A_ALL_dir,method = "BH") 
# X-squared = 32.853, df = 6, p-value = 1.119e-05
# Check text file for posthoc test results

ch_M_ALL_dir<-chisq.test(Manz_intergen$TE_pres, Manz_intergen$Zdir, correct=FALSE)
ch_M_ALL_dir
ch_M_ALL_dir$expected
ch_M_ALL_dir$observed
tab_M_ALL_dir<-table(Manz_intergen$TE_pres, Manz_intergen$Zdir)
chisq.posthoc.test(tab_M_ALL_dir,method = "BH") 
# X-squared = 44.929, df = 6, p-value = 4.834e-08
# Check text file for posthoc test results

################################################################
# 8.6 ggplots for chi-square tests from 8.3 and 8.5

# ggplot theme for these figures
plot_theme3 <- theme(axis.title.x = element_text(size=14),
                    axis.title.y = element_text(size=14),
                    axis.text.x = element_text(size=12),
                    axis.text.y = element_text(size=12),
                    legend.title = element_text(size=14),
                    legend.text = element_text(size=12),
                    legend.position="none")

# Akaa TE associations with DE direction
ggplot(Akaa_intergen[Akaa_intergen$F3A_DE!="Not_DE",], aes(x = TE_pres, fill = Rdir)) + 
  geom_bar(position = "fill")+
  ggtitle("Akaa")+
  geom_text(stat = "count", aes(label = after_stat(count), 
                                y = after_stat(count)/sum(after_stat(count)), 
                                group=Rdir),
            position=position_fill(vjust = c(0.5,0.5)))+
  labs(y = "Proportion", x = "Type of TE present", fill = "RNA direction")+
  scale_fill_manual(values = c("grey90", "grey60"), labels = c("Exp Dn","Exp Up"))+
  scale_x_discrete(labels = c("No TE","Reference","Non Ref"))+
  theme_classic()+
  plot_theme3

# Akaa TE associations with DA direction
ggplot(Akaa_intergen[Akaa_intergen$F3A_DA!="Not_DA",], aes(x = TE_pres, fill = Adir)) + 
  geom_bar(position = "fill")+
  ggtitle("Akaa")+
  geom_text(stat = "count", aes(label = after_stat(count), 
                                y = after_stat(count)/sum(after_stat(count)), 
                                group=Adir),
            position=position_fill(vjust = c(0.5,0.5)))+
  labs(y = "Proportion", x = "Type of TE present", fill = "ATAC direction")+
  scale_fill_manual(values = c("grey90", "grey60"), labels = c("Acc Dn","Acc Up"))+
  scale_x_discrete(labels = c("No TE","Reference","Non Ref"))+
  theme_classic()+
  plot_theme3

## Akaa TE associations with all genes expression & accessibility direction
ggplot(Akaa_intergen, aes(x = TE_pres, fill = Zdir)) + 
  geom_bar(position = "fill")+
  ggtitle("Akaa")+
  geom_text(stat = "count", aes(label = after_stat(count), 
                                y = after_stat(count)/sum(after_stat(count)), 
                                group=Zdir),
            position=position_fill(vjust = c(0.5,0.5,0.5,0.5)))+
  labs(y = "Proportion", x = "Type of TE present", fill = "Correspondence between\nExpression and Accessibility")+
  scale_fill_manual(values = c("#e7d4e8","#9970ab", "#d9f0d3", "#5aae61"), labels = c("Discord Exp Dn Acc Up", "Discord Exp Up Acc Dn", "Concord Exp Dn Acc Dn", "Concord Exp Up Acc Up"))+
  scale_x_discrete(labels = c("No TE","Reference","Non Ref"))+
  theme_classic()+
  plot_theme3

# Manz TE associations with DE direction
ggplot(Manz_intergen[Manz_intergen$F3M_DE!="Not_DE",], aes(x = TE_pres, fill = Rdir)) + 
  geom_bar(position = "fill")+
  ggtitle("Manz")+
  geom_text(stat = "count", aes(label = after_stat(count), 
                                y = after_stat(count)/sum(after_stat(count)), 
                                group=Rdir),
            position=position_fill(vjust = c(0.5,0.5)))+
  labs(y = "Proportion", x = "Type of TE present", fill = "Expression direction")+
  scale_fill_manual(values = c("grey90", "grey60"), labels = c("Exp Dn","Exp Up"))+
  scale_x_discrete(labels = c("No TE","Reference","Non Ref"))+
  theme_classic()+
  plot_theme3

# Manz TE associations with DE direction
ggplot(Manz_intergen[Manz_intergen$F3M_DA!="Not_DA",], aes(x = TE_pres, fill = Adir)) + 
  geom_bar(position = "fill")+
  ggtitle("Manz")+
  geom_text(stat = "count", aes(label = after_stat(count), 
                                y = after_stat(count)/sum(after_stat(count)), 
                                group=Adir),
            position=position_fill(vjust = c(0.5,0.5)))+
  labs(y = "Proportion", x = "Type of TE present", fill = "Accessibility direction")+
  scale_fill_manual(values = c("grey90", "grey60"), labels = c("Acc Dn","Acc Up"))+
  scale_x_discrete(labels = c("No TE","Reference","Non Ref"))+
  theme_classic()+
  plot_theme3

# Manz TE associations with all genes expression & accessibility direction
ggplot(Manz_intergen, aes(x = TE_pres, fill = Zdir)) + 
  geom_bar(position = "fill")+
  ggtitle("Manz")+
  geom_text(stat = "count", aes(label = after_stat(count), 
                                y = after_stat(count)/sum(after_stat(count)), 
                                group=Zdir),
            position=position_fill(vjust = c(0.5,0.5,0.5,0.5)))+
  labs(y = "Proportion", x = "Type of TE present", fill = "Correspondence between\nRNA and ATAC")+
  scale_fill_manual(values = c("#e7d4e8","#9970ab", "#d9f0d3", "#5aae61"), labels = c("Discord Exp Dn Acc Up", "Discord Exp Up Acc Dn", "Concord Exp Dn Acc Dn", "Concord Exp Up Acc Up"))+
  scale_x_discrete(labels = c("No TE","Reference","Non Ref"))+
  theme_classic()+
  plot_theme3

################################################################
# 8.7 enrichment analysis based on genes associated with TEs
#     and showing more than expected association with 
#     expression or accessibility 

# Subset the gene ids based on the interesting association between TE and expression or accessibility
Ak_TED_ED<-subset(Akaa_intergen, ( TE_pres == "D" & F3A_DE != "Not_DE" & RNAseq.log2FoldChange.x < 0))$ensembl_gene_id
Ma_TED_CU<-subset(Manz_intergen, ( TE_pres == "D" & F3M_DA != "Not_DA" & atac.log2FoldChange.x > 0 ))$ensembl_gene_id
Ma_TER_CU<-subset(Manz_intergen, ( TE_pres == "R" & F3M_DA != "Not_DA" & atac.log2FoldChange.x > 0 ))$ensembl_gene_id

## Get list of all gene ids that appear in our analysis. This 
# will be 'universe' for GO enrichment analysis
DmelBG_A <- Akaa_intergen$ensembl_gene_id
DmelBG_M <- Manz_intergen$ensembl_gene_id
DmelBG_All<-unique(sort(c(DmelBG_A,DmelBG_M)))

# Akaa non-ref TE associated with reduced F3 expression
Ak_TED_ED_enrich <- enrichGO(gene = Ak_TED_ED ,
                         OrgDb = org.Dm.eg.db,keyType = 'FLYBASE', readable = T,
                         ont = "BP",
                         universe = DmelBG_All,
                         pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.05)
sum(Ak_TED_ED_enrich@result$p.adjust < 0.05)
#  108
Ak_TED_ED_enrich_scores <- setNames(-log10(Ak_TED_ED_enrich$p.adjust), Ak_TED_ED_enrich$ID)
Ak_TED_ED_enrich_simtrx <- calculateSimMatrix(Ak_TED_ED_enrich$ID, orgdb="org.Dm.eg.db",ont="BP",method="Rel")
Ak_TED_ED_enrich_redtrm <- reduceSimMatrix(Ak_TED_ED_enrich_simtrx,Ak_TED_ED_enrich_scores,threshold=0.7,orgdb="org.Dm.eg.db")
treemapPlot(Ak_TED_ED_enrich_redtrm)
Ak_TED_ED_merge<-merge(Ak_TED_ED_enrich@result[Ak_TED_ED_enrich@result$p.adjust < 0.05,], Ak_TED_ED_enrich_redtrm, by.x = "ID", by.y = "go")
Ak_TED_ED_ord<-Ak_TED_ED_merge[with(Ak_TED_ED_merge, order(Ak_TED_ED_merge$cluster, -Ak_TED_ED_merge$parentSimScore)),]
Ak_TED_ED_out = Ak_TED_ED_ord[!duplicated(Ak_TED_ED_ord$cluster),]
write.table(Ak_TED_ED_out,file = "interchromate_GOenrich/Ak_TED_ED_out.txt", sep = "\t")

# Manz non-ref TE associated with increased F3 expression
Ma_TED_CU_enrich <- enrichGO(gene = Ma_TED_CU ,
                         OrgDb = org.Dm.eg.db,keyType = 'FLYBASE', readable = T,
                         ont = "BP",
                         universe = DmelBG_All,
                         pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.05)
sum(Ma_TED_CU_enrich@result$p.adjust < 0.05)
# 74
Ma_TED_CU_enrich_scores <- setNames(-log10(Ma_TED_CU_enrich$p.adjust), Ma_TED_CU_enrich$ID)
Ma_TED_CU_enrich_simtrx <- calculateSimMatrix(Ma_TED_CU_enrich$ID, orgdb="org.Dm.eg.db",ont="BP",method="Rel")
Ma_TED_CU_enrich_redtrm <- reduceSimMatrix(Ma_TED_CU_enrich_simtrx,Ma_TED_CU_enrich_scores,threshold=0.7,orgdb="org.Dm.eg.db")
treemapPlot(Ma_TED_CU_enrich_redtrm)
Ma_TED_CU_merge<-merge(Ma_TED_CU_enrich@result[Ma_TED_CU_enrich@result$p.adjust < 0.05,], Ma_TED_CU_enrich_redtrm, by.x = "ID", by.y = "go")
Ma_TED_CU_ord<-Ma_TED_CU_merge[with(Ma_TED_CU_merge, order(Ma_TED_CU_merge$cluster, -Ma_TED_CU_merge$parentSimScore)),]
Ma_TED_CU_out = Ma_TED_CU_ord[!duplicated(Ma_TED_CU_ord$cluster),]
write.table(Ma_TED_CU_out,file = "interchromate_GOenrich/Ma_TED_CU_out.txt", sep = "\t")

# Manz ref TE associated with increased F3 expression
Ma_TER_CU_enrich <- enrichGO(gene = Ma_TER_CU,
                         OrgDb = org.Dm.eg.db,keyType = 'FLYBASE', readable = T,
                         ont = "BP",
                         universe = DmelBG_All,
                         pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.05)
sum(Ma_TER_CU_enrich@result$p.adjust < 0.05)
# 3
Ma_TER_CU_enrich_scores <- setNames(-log10(Ma_TER_CU_enrich$p.adjust), Ma_TER_CU_enrich$ID)
Ma_TER_CU_enrich_simtrx <- calculateSimMatrix(Ma_TER_CU_enrich$ID, orgdb="org.Dm.eg.db",ont="BP",method="Rel")
Ma_TER_CU_enrich_redtrm <- reduceSimMatrix(Ma_TER_CU_enrich_simtrx,Ma_TER_CU_enrich_scores,threshold=0.7,orgdb="org.Dm.eg.db")
treemapPlot(Ma_TER_CU_enrich_redtrm)
Ma_TER_CU_merge<-merge(Ma_TER_CU_enrich@result[Ma_TER_CU_enrich@result$p.adjust < 0.05,], Ma_TER_CU_enrich_redtrm, by.x = "ID", by.y = "go")
Ma_TER_CU_ord<-Ma_TER_CU_merge[with(Ma_TER_CU_merge, order(Ma_TER_CU_merge$cluster, -Ma_TER_CU_merge$parentSimScore)),]
Ma_TER_CU_out = Ma_TER_CU_ord[!duplicated(Ma_TER_CU_ord$cluster),]
write.table(Ma_TER_CU_out,file = "interchromate_GOenrich/Ma_TER_CU_out.txt", sep = "\t")

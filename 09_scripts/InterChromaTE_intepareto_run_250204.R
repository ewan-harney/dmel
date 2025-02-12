library(intePareto)
library(data.table)
library(apeglm)

rm(list=ls(all=TRUE))

## Set working directory and load customised doMatch function
setwd("C:/Users/ewanh/Dropbox/Bioinformatics_UPF/ATAC_DS_analysis")
source("doMatch_4_function.R")

## Load annotation data (see InterChromaTE_intepareto_prepTranscriptTable_250204.sh)
annotation_data <- fread("annotation_dm_r6_v46.tsv", header = F, stringsAsFactors = F)
colnames(annotation_data) <- c("ensembl_transcript_id", "ensembl_gene_id",
                               "external_gene_name", "chromosome_name", "start_position",
                               "end_position", "transcript_start", "transcript_end",
                               "transcription_start_site", "transcript_length", "strand")

## What comparison are we doing?
## These strings will be added to the object name at the end
comparison <-"AkaaF3_Heat_v_AkaaF3_Ctrl"
promoterLength <- as.numeric(1000)

## Where to read the kallisto *abundance.tsv files from
rna_data <- read.table("rna_meta_v2.txt", header = F, stringsAsFactors = F)
colnames(rna_data) <- c("condition", "files")
rna_data <- rna_data[rna_data$condition=="AkaaF3Heat" | rna_data$condition=="AkaaF3Ctrl",]

## Where to read the nf-core/atacseq *mLb.clN.sorted.bam files from
atac_data <- read.table("atac_meta_v2.txt", header = F, stringsAsFactors = F)
colnames(atac_data) <- c("mark", "condition", "files")
atac_data <- atac_data[atac_data$condition=="AkaaF3Heat" | atac_data$condition=="AkaaF3Ctrl",]

## Use the doMatch_4 function (a slightly altered version of the base doMatch function in intepareto) to match rna and atac reads based on postions (assuming atac data is for promoters within 1kb of expressed genes)
res <- doMatch_4(rnaMeta = rna_data,
                 chipMeta = atac_data,
                 region = "promoter",
                 promoter = 1000,
                 fragLength = 50,
                 method = "weighted.mean",
                 annotation = annotation_data,
                 conditions = c("AkaaF3Heat","AkaaF3Ctrl"))

## Matched data object will be used for ATAC-seq DESeq2 analysis
matchedData <- res[["matched.data"]]

## Calculate a Z score with the doIntegration function
df_final <- doIntegration(res = res, # result list from "doMatch" function
                          type = "apeglm", # shrinkage estimator, default is "apeglm"
                          ref = "AkaaF3Ctrl", # specifying the reference level
                          apeAdapt = FALSE)

## Extract just the matched ensemble gene IDs, RNA and ATAC l2fc values and the z score
Integ_zscores<- na.omit(merge(matchedData,df_final, by="row.names",all.x=TRUE)[,c(2,15:17)])

## Save the two objects to an Rdata object
save(matchedData,Integ_zscores, file =paste0("interchromate_intepareto/interchromate_",comparison,promoterLength,".Rdata"))


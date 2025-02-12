#!/bin/bash

#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ewan.harney@gmail.com
#SBATCH --cpus-per-task=12
#SBATCH --partition=normal
#SBATCH --job-name=createTranscriptsTable
#SBATCH -o logs/createTranscriptsTable_%a.out # Standard output 
#SBATCH -e logs/createTranscriptsTable_%a.err # Standard error  

DIR="/homes/users/eharney/scratch/intepareto"

# Creat the file to be populated with transcripts
>  ${DIR}/annotation_dm_r6_v46.tsv

# Select only features annotated as mRNA from the .gtf file
transcripts=$(awk '$3 == "mRNA"' ${DIR}/Dmel_FB_r6_v46_chr.gtf |cut -f9 | cut -f 6 -d ' ' | sort -u |tr -d '"' |tr -d ';')

# Function to gather information on transcripts from different fields
while read transcript
do
gene=$(grep "$transcript" ${DIR}/Dmel_FB_r6_v46_chr.gtf | awk ' $3 == "mRNA" '  | cut -f 9 | cut -f2 -d' ' |tr -d '"' |tr -d ';' )
symbol=$(grep "$transcript" ${DIR}/Dmel_FB_r6_v46_chr.gtf | awk ' $3 == "mRNA" ' | cut -f 9 | cut -f4 -d' ' |tr -d '"' |tr -d ';' )
chr=$(grep "$transcript" ${DIR}/Dmel_FB_r6_v46_chr.gtf | awk ' $3 == "mRNA" ' | cut -f 1 )
geneStart=$(grep "$gene" ${DIR}/Dmel_FB_r6_v46_chr.gtf | awk ' $3 == "gene" ' |cut -f4)
geneEnd=$(grep "$gene" ${DIR}/Dmel_FB_r6_v46_chr.gtf | awk ' $3 == "gene" ' |cut -f5)
transcriptStart=$(grep "$transcript" ${DIR}/Dmel_FB_r6_v46_chr.gtf | awk ' $3 == "mRNA" ' | cut -f 4 )
transcriptEnd=$(grep "$transcript" ${DIR}/Dmel_FB_r6_v46_chr.gtf | awk ' $3 == "mRNA" ' | cut -f 5 )
strand=$(grep "$transcript" ${DIR}/Dmel_FB_r6_v46_chr.gtf | awk ' $3 == "mRNA" ' | cut -f 7 )
if [[ $strand == "+" ]]; then
#statements
strand="1"
TSSstart=$transcriptStart
elif [[ $strand == "-" ]]; then
strand="-1"
TSSstart=$transcriptEnd
fi
transcriptLength=$(grep "$transcript" ${DIR}/Dmel_FB_r6_v46_chr.gtf | awk ' $3 == "exon" ' | awk 'BEGIN { FS=OFS="\t" } {print $5-$4+1}' |  awk '{ sum += $0} END {print sum}' )
echo -e "$transcript\t$gene\t$symbol\t$chr\t$geneStart\t$geneEnd\t$transcriptStart\t$transcriptEnd\t$TSSstart\t$transcriptLength\t$strand" >>  ${DIR}/annotation_dm_r6_v46.tsv
done <<< "$transcripts"


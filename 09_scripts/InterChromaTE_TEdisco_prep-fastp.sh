#!/bin/bash

#SBATCH --partition=normal
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --mail-type=START,END,FAIL
#SBATCH -o logs/intePareto_%a.out # Standard output
#SBATCH -e logs/intePareto_%a.err # Standard error
#SBATCH --mail-user=ewan.harney@gmail.com

module load fastp

line=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat DNA_fastp_filelist.txt))
INDIR=/gpfs42/projects/lab_jgonzalez/gonzalez_lab/raw_data/InterchromaTE/WP1/220329_WP1_DNA
DIR=/gpfs42/robbyfs/scratch/lab_jgonzalez/eharney/fastp
sample=$(echo "$line" | cut -f1)
output=$(echo "$line" | cut -f2)
fastp -w 12 -h ${DIR}/${output}fastp.html \
-i ${INDIR}/${sample}1_001.fastq.gz -I ${INDIR}/${sample}2_001.fastq.gz \
-o ${DIR}/${output}1_fastp.fastq \
-O ${DIR}/${output}2_fastp.fastq

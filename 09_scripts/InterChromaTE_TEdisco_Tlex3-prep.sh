#!/bin/bash

#SBATCH --partition=normal
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --mail-type=START,END,FAIL
#SBATCH -o Tlex/logs/tlex_prep_jan23_%a.out # Standard output
#SBATCH -e Tlex/logs/tlex_prep_jan23_%a.err # Standard error
#SBATCH --mail-user=ewan.harney@gmail.com

DIR="/gpfs42/robbyfs/scratch/lab_jgonzalez/eharney"

## Make the directories that tlex will read
mkdir -p ${DIR}/Tlex/F3-P11/F3-P11
mkdir -p ${DIR}/Tlex/F3-P12/F3-P12
mkdir -p ${DIR}/Tlex/F3-P13/F3-P13
mkdir -p ${DIR}/Tlex/F3-P21/F3-P21
mkdir -p ${DIR}/Tlex/F3-P22/F3-P22
mkdir -p ${DIR}/Tlex/F3-P23/F3-P23

## Stay in the same directory. Copy the files from fastp folder and rename to be accepted by tlex
# Copy and rename R1 files
cp ${DIR}/fastp/F3_P11_R1_fastp.fastq.gz ${DIR}/Tlex/F3-P11/F3-P11/F3-P11_1.fastq.gz
cp ${DIR}/fastp/F3_P12_R1_fastp.fastq.gz ${DIR}/Tlex/F3-P12/F3-P12/F3-P12_1.fastq.gz
cp ${DIR}/fastp/F3_P13_R1_fastp.fastq.gz ${DIR}/Tlex/F3-P13/F3-P13/F3-P13_1.fastq.gz
cp ${DIR}/fastp/F3_P21_R1_fastp.fastq.gz ${DIR}/Tlex/F3-P21/F3-P21/F3-P21_1.fastq.gz
cp ${DIR}/fastp/F3_P22_R1_fastp.fastq.gz ${DIR}/Tlex/F3-P22/F3-P22/F3-P22_1.fastq.gz
cp ${DIR}/fastp/F3_P23_R1_fastp.fastq.gz ${DIR}/Tlex/F3-P23/F3-P23/F3-P23_1.fastq.gz
# R2 files
cp ${DIR}/fastp/F3_P11_R2_fastp.fastq.gz ${DIR}/Tlex/F3-P11/F3-P11/F3-P11_2.fastq.gz
cp ${DIR}/fastp/F3_P12_R2_fastp.fastq.gz ${DIR}/Tlex/F3-P12/F3-P12/F3-P12_2.fastq.gz
cp ${DIR}/fastp/F3_P13_R2_fastp.fastq.gz ${DIR}/Tlex/F3-P13/F3-P13/F3-P13_2.fastq.gz
cp ${DIR}/fastp/F3_P21_R2_fastp.fastq.gz ${DIR}/Tlex/F3-P21/F3-P21/F3-P21_2.fastq.gz
cp ${DIR}/fastp/F3_P22_R2_fastp.fastq.gz ${DIR}/Tlex/F3-P22/F3-P22/F3-P22_2.fastq.gz
cp ${DIR}/fastp/F3_P23_R2_fastp.fastq.gz ${DIR}/Tlex/F3-P23/F3-P23/F3-P23_2.fastq.gz
# gunzip R1 files
gunzip ${DIR}/Tlex/F3-P11/F3-P11/F3-P11_1.fastq.gz
gunzip ${DIR}/Tlex/F3-P12/F3-P12/F3-P12_1.fastq.gz
gunzip ${DIR}/Tlex/F3-P13/F3-P13/F3-P13_1.fastq.gz
gunzip ${DIR}/Tlex/F3-P21/F3-P21/F3-P21_1.fastq.gz
gunzip ${DIR}/Tlex/F3-P22/F3-P22/F3-P22_1.fastq.gz
gunzip ${DIR}/Tlex/F3-P23/F3-P23/F3-P23_1.fastq.gz
# R2 files
gunzip ${DIR}/Tlex/F3-P11/F3-P11/F3-P11_2.fastq.gz
gunzip ${DIR}/Tlex/F3-P12/F3-P12/F3-P12_2.fastq.gz
gunzip ${DIR}/Tlex/F3-P13/F3-P13/F3-P13_2.fastq.gz
gunzip ${DIR}/Tlex/F3-P21/F3-P21/F3-P21_2.fastq.gz
gunzip ${DIR}/Tlex/F3-P22/F3-P22/F3-P22_2.fastq.gz
gunzip ${DIR}/Tlex/F3-P23/F3-P23/F3-P23_2.fastq.gz


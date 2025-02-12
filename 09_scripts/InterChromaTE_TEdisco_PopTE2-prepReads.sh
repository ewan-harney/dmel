#!/bin/bash

#SBATCH --partition=normal
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --mail-type=START,END,FAIL
#SBATCH -o logs/preplog_P1P2_redo_%a.out # Standard output
#SBATCH -e logs/preplog_P1P2_redo_%a.err # Standard error
#SBATCH --mail-user=ewan.harney@gmail.com

module load PoPoolationTE2/v1.10.03
module load BWA/0.7.17-GCC-9.3.0

# Run with:
# sbatch --array=1-6 PopTE2_prep_v2_redo_v3_P1P2.sh
DIR="/gpfs42/robbyfs/scratch/lab_jgonzalez/eharney"

# Get the sample
data=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat DNA_bwamem_filelist_P1P2.txt))
sample=$(echo "$data" | cut -f2 )
site=$(echo "$data" | cut -f3 )
genome="${DIR}/PopTE2/input/Dmel_FB_r6_v46/Dmel_FB_r6_v46_MASK.chr.TE.merged.fasta"

#### Here is the pre-processing of the fastq
# First I copy reads
read1=$(echo "${DIR}/fastp/${sample}1_fastp.fastq.gz")
read2=$(echo "${DIR}/fastp/${sample}2_fastp.fastq.gz")

cp ${read1} ${DIR}/PopTE2/fastq/
cp ${read2} ${DIR}/PopTE2/fastq/

# Unzip 
gunzip ${DIR}/PopTE2/fastq/${sample}1_fastp.fastq.gz &
gunzip ${DIR}/PopTE2/fastq/${sample}2_fastp.fastq.gz

# Modify the Illumina header
sed 's/ 1:N:0:\(.*\)$/\#\1\/1/g' ${DIR}/PopTE2/fastq/${sample}1_fastp.fastq > ${DIR}/PopTE2/fastq/${sample}1_id_fastp.fastq
sed 's/ 2:N:0:\(.*\)$/\#\1\/2/g' ${DIR}/PopTE2/fastq/${sample}2_fastp.fastq > ${DIR}/PopTE2/fastq/${sample}2_id_fastp.fastq

# Zip again
gzip ${DIR}/PopTE2/fastq/${sample}1_id_fastp.fastq
gzip ${DIR}/PopTE2/fastq/${sample}2_id_fastp.fastq

# Remove
rm -f ${DIR}/PopTE2/fastq/${sample}1_fastp.fastq
rm -f ${DIR}/PopTE2/fastq/${sample}2_fastp.fastq

####

# Probably you can start from here directly:
echo $read1 $read2 $genome $sample

#  Map reads to the TE-combined-reference
bwa bwasw \
-t 12 \
${genome} \
${DIR}/PopTE2/fastq/${sample}1_id_fastp.fastq.gz > \
${DIR}/PopTE2/map/local/${sample}1.sam

bwa bwasw \
-t 12 \
${genome} \
${DIR}/PopTE2/fastq/${sample}2_id_fastp.fastq.gz > \
${DIR}/PopTE2/map/local/${sample}2.sam

# Restore paired-end information with PoPoolationTE2 se2pe

java  -jar $EBROOTPOPTE/popte2-v1.10.03.jar se2pe \
--fastq1 ${DIR}/PopTE2/fastq/${sample}1_id_fastp.fastq.gz  \
--fastq2 ${DIR}/PopTE2/fastq/${sample}2_id_fastp.fastq.gz \
--bam1 ${DIR}/PopTE2/map/local/${sample}1.sam \
--bam2 ${DIR}/PopTE2/map/local/${sample}2.sam \
--sort \
--detailed-log \
--output ${DIR}/PopTE2/map/local/${sample}.sort.bam

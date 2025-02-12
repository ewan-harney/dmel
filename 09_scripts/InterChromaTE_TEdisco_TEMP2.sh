#!/bin/bash

#SBATCH --partition=normal
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --mail-type=START,END,FAIL
#SBATCH -o logs/TEMP2_flybase_PICTEMP_230125_v2_%a.out # Standard output
#SBATCH -e logs/TEMP2_flybase_PICTEMP_230125_v2_%a.err # Standard error
#SBATCH --mail-user=ewan.harney@gmail.com

module load picard/2.25.5-Java-11
module load TEMP2/0.1.4-foss-2016b

# PRESENCE MODULE
module load BEDTools/2.27.1-foss-2016b

# sbatch --array=1-67 TEMP.sh
DIR="/homes/users/eharney/scratch"

cd $DIR/TEMP2

data=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat DNA_TEMP2_flybase_P12.txt))
sample=$(echo "$data" | cut -f2 )
#site=$(echo "$data" | cut -f3 )
genome="${DIR}/genomes/bwa_220706/Dmel_FB_r6_v46_chr.fasta"
repeatmasker="${DIR}/TEMP2/input/Dmel_FB_r6_v46_chr.fasta.out.bed"
genome2bit="${DIR}/TEMP2/input/Dmel_FB_r6_v46_chr.fasta.2bit"

read1=$(echo "${DIR}/fastp/${sample}1_fastp.fastq.gz")
read2=$(echo "${DIR}/fastp/${sample}2_fastp.fastq.gz")
bam=$(echo ${DIR}/TEMP2/map2/${sample}.flybase.sorted.bam)

echo $read1 $read2 $genome $sample

mkdir -p ${DIR}/TEMP2/output2/${sample}_v2/flybase/insertion_m5

TEMP2 insertion \
-l ${read1} \
-r ${read2} \
-i ${bam} \
-g ${genome} \
-R ${DIR}/TEs/consensuses_curated_v4.fasta \
-t ${repeatmasker} \
-o ${DIR}/TEMP2/output2/${sample}_v2/flybase/insertion_m5 \
-p ${sample} \
-c 12 \
-m 5 

#!/bin/bash

#SBATCH --partition=normal
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --mail-type=START,END,FAIL
#SBATCH -o logs/kquant_%a.out # Standard output
#SBATCH -e logs/kquant_%a.err # Standard error
#SBATCH --mail-user=ewan.harney@gmail.com


### use the command sbatch --array=1-24 kallisto_quant_v2.sh

module load kallisto/0.46.1-foss-2019b

line=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat RNA_data-files_kalisto.txt))
infastq=/scratch/lab_jgonzalez/eharney/RNA-kallisto/fastq
resfolder=/scratch/lab_jgonzalez/eharney/RNA-kallisto/results
index=/scratch/lab_jgonzalez/eharney/RNA-kallisto/Dmel_FB_r6_v46_transcript.idx

output=$(echo "$line" | cut -f1)
reads_1=$(echo "$line" | cut -f2)
reads_2=$(echo "$line" | cut -f3)

echo ${index}
echo ${resfolder}/${output}
echo ${infastq}/${reads_1}
echo ${infastq}/${reads_2}

mkdir ${resfolder}/${output}
kallisto quant -i ${index} -o ${resfolder}/${output} -b 100 ${infastq}/${reads_1} ${infastq}/${reads_2}

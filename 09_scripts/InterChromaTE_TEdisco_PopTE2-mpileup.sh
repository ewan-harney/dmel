#!/bin/bash

#SBATCH --partition=normal
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --mail-type=START,END,FAIL
#SBATCH -o logs/ppileup_230123_%a.out # Standard output
#SBATCH -e logs/ppileup_230123_%a.err # Standard error
#SBATCH --mail-user=ewan.harney@gmail.com

module load PoPoolationTE2/v1.10.03

DIR="/gpfs42/robbyfs/scratch/lab_jgonzalez/eharney"

cd $DIR/PopTE2/

data=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat PopTE2_miltippileup_list_redo2.txt))
mapref=$(echo "$data" | cut -f1 )
sample1=$(echo "$data" | cut -f2 )
sample2=$(echo "$data" | cut -f3 )
sample3=$(echo "$data" | cut -f4 )
outname=$(echo "$data" | cut -f5 )

java -Xmx32g -jar $EBROOTPOPTE/popte2-v1.10.03.jar ppileup \
--bam ${DIR}/PopTE2/map/${mapref}/${sample1} \
--bam ${DIR}/PopTE2/map/${mapref}/${sample2} \
--bam ${DIR}/PopTE2/map/${mapref}/${sample3} \
--map-qual 15 \
--hier ${DIR}/TEs/te-hierarchy.txt \
--output ${DIR}/PopTE2/output/${mapref}/ppileup/${outname}.${mapref}ref.ppileup.gz \
--detailed-log

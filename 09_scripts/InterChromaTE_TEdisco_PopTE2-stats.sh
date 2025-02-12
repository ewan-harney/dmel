#!/bin/bash

#SBATCH --partition=normal
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --mail-type=START,END,FAIL
#SBATCH -o logs/popte2_stats_230214_%a.out # Standard output
#SBATCH -e logs/popte2_stats_230214_%a.err # Standard error
#SBATCH --mail-user=ewan.harney@gmail.com

module load PoPoolationTE2/v1.10.03

# run with array like so
# sbatch --array=1-6 interChromaTE_TEdisco_PopoolationTE2.sh

DIR="/gpfs42/robbyfs/scratch/lab_jgonzalez/eharney"

cd $DIR/PopTE2/

data=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat PopTE2_multi_stats_list_v2.txt))
# refgen=$(echo "$data" | cut -f1 )
simpleref=$(echo "$data" | cut -f2 )
# sample1=$(echo "$data" | cut -f3 )
# sample2=$(echo "$data" | cut -f4 )
# sample3=$(echo "$data" | cut -f5 )
outname=$(echo "$data" | cut -f6 )

genome="${DIR}/PopTE2/input/Dmel_FB_r6_v46/Dmel_FB_r6_v46_MASK.chr.TE.merged.fasta"

# identifySignatures: the important parameter here is min-count, which is the minimum coverage an insertion breakpoint has to be called
java -Xmx64G -jar $EBROOTPOPTE/popte2-v1.10.03.jar identifySignatures \
--ppileup ${DIR}/PopTE2/output/${simpleref}/ppileup/${outname}.${simpleref}ref.ppileup.gz \
--mode separate \
--output ${DIR}/PopTE2/output2/mincount2/${outname}.${simpleref}ref.signatures \
--min-count 2 \
--signature-window minimumSampleMedian \
--min-valley minimumSampleMedian \
--detailed-log

# Estimate the strand - this is optional, and I didn't end up doing it
#java -Xmx64G -jar $EBROOTPOPTE/popte2-v1.10.03.jar updatestrand \
#--bam ${DIR}/PopTE2/map/${simpleref}/${sample1} \
#--bam ${DIR}/PopTE2/map/${simpleref}/${sample2} \
#--bam ${DIR}/PopTE2/map/${simpleref}/${sample3} \
#--signature ${DIR}/PopTE2/output/mincount2/${outname}.${simpleref}ref.signatures \
#--hier ${DIR}/TEs/te-hierarchy.txt \
#--map-qual 15 \
#--max-disagreement 0.5 \
#--output ${DIR}/PopTE2/output/mincount2/${outname}.${simpleref}ref.strand.signatures

# frequency
java -Xmx64G -jar $EBROOTPOPTE/popte2-v1.10.03.jar frequency \
--ppileup ${DIR}/PopTE2/output/${simpleref}/ppileup/${outname}.${simpleref}ref.ppileup.gz \
--signature ${DIR}/PopTE2/output2/mincount2/${outname}.${simpleref}ref.signatures \
--output ${DIR}/PopTE2/output2/mincount2/${outname}.${simpleref}ref.strand.freqsig \
--detailed-log

# pairupSignatures
java -Xmx64G -jar $EBROOTPOPTE/popte2-v1.10.03.jar pairupSignatures \
--signature ${DIR}/PopTE2/output2/mincount2/${outname}.${simpleref}ref.strand.freqsig  \
--ref-genome ${genome} \
--hier ${DIR}/TEs/te-hierarchy.txt \
--min-distance -200 \
--max-distance 300 \
--output ${DIR}/PopTE2/output2/mincount2/${outname}.${simpleref}ref.strand.teinsertions \
--detailed-log \
--output-detail medium

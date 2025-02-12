#!/bin/bash

#SBATCH --partition=normal
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --mail-type=START,END,FAIL
#SBATCH -o logs/rmlog_%a.out # Standard output
#SBATCH -e logs/rmlog_%a.err # Standard error
#SBATCH --mail-user=ewan.harney@gmail.com

## Repeat mask the reference genome with RepeatMasker
module load RepeatMasker/4.1.2-GCC-9.3.0-HMMER

DIR=/gpfs42/robbyfs/scratch/lab_jgonzalez/eharney

RepeatMasker -e wublast -gccalc -s -cutoff 200 -no_is -nolow -norna -gff -u -pa 12 \
-lib $DIR/TEs/consensuses_curated_v4.fasta \
-dir $DIR/PopTE2/input/Dmel_FB_r6_v46 \
$DIR/genomes/bwa_220706/Dmel_FB_r6_v46_chr.fasta

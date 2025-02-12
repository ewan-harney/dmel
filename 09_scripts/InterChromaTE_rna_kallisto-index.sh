#!/bin/bash

#SBATCH --partition=normal
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --mail-type=START,END,FAIL
#SBATCH -o logs/index_%a.out # Standard output
#SBATCH -e logs/index_%a.err # Standard error
#SBATCH --mail-user=ewan.harney@gmail.com

module load kallisto/0.46.1-foss-2019b

kallisto index -i Dmel_FB_r6_v46_transcript.idx dmel-all-transcript-r6.46.fasta.gz
